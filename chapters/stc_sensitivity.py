"""Isolated LTspice irradiance sensitivity and audit of temperature support.

Run: python chapters/stc_sensitivity.py
Original model/input/manuscript outputs are never edited. Native saved-result
reproduction precedes a separate explicit-25-C isothermal sensitivity. The
existing t/TNOM coupling does not support a defensible temperature sweep.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from soilspice import config, experiments, geometry, modules, spice
from soilspice.style import STYLE, ORDER, ORIENTATIONS, module_label, save
from chapters.longi_reverse_bias_diagnostic import read_raw

GS = (200, 400, 600, 800, 1000)
TAUS = (.71, .44, 0.)
COLORS = {.71: '#CF8425', .44: '#6F300F', 0.: '#26130B'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@dataclass(frozen=True)
class State:
    module: str
    shape: str = 'clean'
    orientation: str = 'portrait'
    h: float = 0.
    tau: float = 1.
    G: int = 1000
    track: str = 'explicit_25C'

    @property
    def key(self):
        return f'{self.module}_{self.shape}_{self.orientation}_h{self.h:g}_tau{self.tau:g}_G{self.G}_{self.track}'


def cases():
    return [(name, shape, orientation, h, tau)
            for name in ORDER for shape, orientation, h, taus in
            [(s,o,80.,TAUS) for s in ('rectangle','triangle') for o in ORIENTATIONS]
            + [('rectangle',o,60.,(.44,)) for o in ORIENTATIONS]
            for tau in taus]


def netlist(state, library):
    """Use repository patching, adding only operating controls and trace output.

    legacy: preserve fitted t_ref, implicit native operating temperature and
    original sweep. explicit_25C: t=25 AND .temp 25, Vi=G, dense global sweep.
    No IS/n/Rs/Rp/Bishop/bypass coefficients or topology are changed.
    """
    module=modules.get(state.module)
    original=spice.read_netlist(module.netlist_path).replace('\r\n','\n')
    text=re.sub(r'^\.lib .+$',lambda _:f'.lib "{library.resolve()}"',original,flags=re.M)
    text=re.sub(r'^(V1\s+\S+\s+\S+)\s+V\s*$',r'\1 0',text,flags=re.M)
    if state.shape=='clean':
        cover={ref:0. for ref in module.instances.values()}
    else:
        grid=geometry.fcov_grid(module,state.shape,state.orientation,state.h,experiments.TRIANGLE_APEX)
        cover=geometry.cell_fcov(module,grid,state.orientation)
    text=spice.patch_cells(text,module,cover,state.tau)
    temperature=module.t_ref if state.track=='legacy' else 25.
    text=spice.set_temperature(text,module,temperature)
    text=spice.set_param(text,'Vi',state.G)
    if state.track!='legacy':
        assert not re.search(r'^\.temp\b',text,flags=re.M|re.I)
        text=text.replace('.end\n','.temp 25\n.end\n')
        sweep=re.search(r'^\.dc\s+V1\s+([\d.]+)\s+([\d.]+)\s+([\d.-]+)',text,flags=re.M|re.I)
        lo,hi=map(float,sweep.groups()[:2])
        step=.01 if state.track=='explicit_25C_refined' else .02
        upper=math.ceil(max(lo,hi)*1.15/step)*step
        directive=f'.dc V1 {upper:.2f} 0 -{step:g}' if lo>hi else f'.dc V1 0 {upper:.2f} {step:g}'
        text=text[:sweep.start()]+directive+text[sweep.end():]
    diode_nodes=[re.search(rf'^D{k}\s+(\S+)\s+(\S+)\s',text,flags=re.M).groups() for k in (1,2,3)]
    source_nodes=re.search(r'^V1\s+(\S+)\s+(\S+)\s',text,flags=re.M).groups()
    requests=['I(V1)',*[f'I(D{k})' for k in (1,2,3)]]
    for node in (*source_nodes,*[node for pair in diode_nodes for node in pair]):
        if node!='0': requests.append(f'V({node})')
    text=text.replace('.end\n','.save '+' '.join(dict.fromkeys(requests))+'\n.end\n')
    return text,diode_nodes,source_nodes


def measure(raw_path, diode_nodes, source_nodes, *, require_voc_bracket=True):
    raw=read_raw(raw_path)
    v,i=raw[config.VOLTAGE_TRACE.lower()],raw[config.CURRENT_TRACE.lower()]*config.CURRENT_SIGN
    order=np.argsort(v)
    v,i=v[order],i[order]
    assert np.isfinite(v).all() and np.isfinite(i).all()
    def node(n):
        return np.zeros_like(v) if n=='0' else raw[f'v({n})'.lower()][order]
    assert np.max(np.abs(node(source_nodes[0])-node(source_nodes[1])-v))<5e-5
    keep=(v>=0)&(i>=0)
    indices=np.flatnonzero(keep)
    if len(indices)<10: raise ValueError('Insufficient first-quadrant curve')
    params=spice.iv_parameters(v[keep],i[keep])
    index=indices[int(np.argmax(v[keep]*i[keep]))]
    if not (params['Pmax_W']>0 and params['Isc_A']>0 and indices[0]<index<indices[-1]):
        raise ValueError('Invalid or sweep-edge MPP')
    crossings=np.flatnonzero((i[:-1]>=0)&(i[1:]<0))
    params['Voc_bracketed']=bool(len(crossings))
    if len(crossings):
        j=crossings[-1]
        params['Voc_V']=float(np.interp(0,i[j:j+2][::-1],v[j:j+2][::-1]))
    elif require_voc_bracket:
        raise ValueError('Sweep does not bracket open circuit')
    # Legacy SunPower exports stop below Voc. Preserve that exact archived
    # sweep for reproduction; require a bracket in every new sensitivity run.
    bypass_mask=0
    for k,(pos,neg) in enumerate(diode_nodes,1):
        current=raw[f'i(d{k})'][order]
        vd=node(pos)-node(neg)
        active=bool(current[index]>max(1e-4,.01*i[index]) and vd[index]>0)
        params[f'bypass_{k}_MPP_A']=float(current[index])
        params[f'bypass_{k}_MPP_voltage_V']=float(vd[index])
        params[f'bypass_{k}_MPP_active']=active
        params[f'bypass_{k}_max_first_quadrant_A']=float(np.max(current[keep]))
        if active: bypass_mask |= 1 << (k-1)
    params['bypass_MPP_mask']=bypass_mask
    params['MPP_grid_step_V']=float(np.median(np.diff(v)))
    return params


class Runner:
    def __init__(self,out,library):
        self.out,self.library=out,library
        self.runs=out/'runs'
        self.runs.mkdir(parents=True,exist_ok=True)
        self.results={}

    def one(self,state):
        start=time.monotonic()
        path=self.runs/(state.key+'.net')
        metadata=path.with_suffix('.json')
        text,diodes,nodes=netlist(state,self.library)
        signature=hashlib.sha256(text.encode()).hexdigest()
        try:
            if metadata.exists() and path.with_suffix('.raw').exists():
                cached=json.loads(metadata.read_text())
                if cached.get('netlist_sha256')==signature and cached.get('status')=='ok':
                    return state,cached
            path.write_text(text,encoding='utf-8')
            raw=spice.run_ltspice(path)
            values=measure(raw,diodes,nodes,require_voc_bracket=state.track!='legacy')
            log=path.with_suffix('.log').read_text(errors='replace')
            temperature=re.search(r'^temp\s*=\s*([\d.+-]+)',log,re.M|re.I)
            result=dict(**asdict(state),status='ok',**values,cell_temperature_C=float(temperature[1]),
                nominal_parameter_temperature_C=modules.get(state.module).t_ref if state.track=='legacy' else 25.,
                netlist_sha256=signature,raw_sha256=digest(raw),netlist=str(path.relative_to(self.out)),
                elapsed_seconds=time.monotonic()-start,solver_warning_count=len(re.findall(r'warning',log,re.I)))
            if state.track!='legacy': assert result['cell_temperature_C']==25
        except Exception as error:
            result=dict(**asdict(state),status='failed',error=repr(error),netlist_sha256=signature,
                        netlist=str(path.relative_to(self.out)),elapsed_seconds=time.monotonic()-start)
        metadata.write_text(json.dumps(result,indent=2)+'\n')
        return state,result

    def batch(self,states,workers,label):
        todo=sorted(set(states)-self.results.keys(),key=lambda s:s.key)
        print(f'{label}: {len(todo)} runs',flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(self.one,state) for state in todo]
            for n,future in enumerate(as_completed(futures),1):
                state,result=future.result()
                self.results[state]=result
                if n%20==0 or n==len(todo):
                    print(f'  {label}: {n}/{len(todo)} completed; '
                          f'{sum(r["status"]!="ok" for r in self.results.values())} failed',flush=True)
        self.flush()

    def flush(self):
        frame=pd.DataFrame(self.results.values())
        frame.to_csv(self.out/'run_manifest.csv',index=False)
        frame[frame.status!='ok'].to_csv(self.out/'failed_simulations.csv',index=False)


def validate(runner,workers):
    saved={shape:pd.read_csv(config.CH1_DIR/'data'/f'iv_grid_{shape}.csv')
           for shape in ('rectangle','triangle')}
    jobs=[State(name,track='legacy') for name in ORDER]
    for name in ORDER:
        for shape,table in saved.items():
            for ori in ORIENTATIONS:
                q=table[(table.module==name)&(table.orientation==ori)&(table.source=='SPICE')]
                available=q.loc[q.tau_I<1,'tau_I'].unique()
                # Several landscape flashes lack the 0.442 campaign state.
                # Reproduce only saved states: use the closest available
                # campaign value and record its exact value in the audit.
                chosen={float(min(available,key=lambda value:abs(value-target))) for target in (.706,.442,0.)}
                jobs += [State(name,shape,ori,80,tau,1000,'legacy') for tau in sorted(chosen)]
    runner.batch(jobs,workers,'Preflight: archived reference')
    rows=[]
    for state in jobs:
        result=runner.results[state]
        shape='rectangle' if state.shape=='clean' else state.shape
        q=saved[shape]
        q=q[(q.module==state.module)&(q.orientation==state.orientation)&(q.source=='SPICE')
            &np.isclose(q.tau_I,state.tau)]
        expected=float((q.V*q.I).max())
        actual=result.get('Pmax_W',np.nan)
        delta=actual-expected
        rows.append(dict(**asdict(state),saved_Pmp_W=expected,reproduced_Pmp_W=actual,
                         difference_W=delta,pass_reproduction=bool(abs(delta)<.02),
                         actual_cell_temperature_C=result.get('cell_temperature_C',np.nan)))
    legacy=pd.DataFrame(rows)
    legacy.to_csv(runner.out/'preflight_archived_reference.csv',index=False)
    if not legacy.pass_reproduction.all():
        raise RuntimeError('Archived-reference preflight failed; see preflight CSV and failed_simulations.csv')
    # Test true 25 C separately, against the same archived states; this must
    # not be conflated with reproducing the original implicit-27-C results.
    requested=[State(s.module,s.shape,s.orientation,s.h,s.tau,s.G) for s in jobs]
    runner.batch(requested,workers,'Preflight: explicit 25 C versus archive')
    records=[]
    for old,new in zip(jobs,requested):
        a,b=runner.results[old],runner.results[new]
        records.append(dict(**asdict(new),archived_Pmp_W=a.get('Pmax_W'),
            explicit_25C_Pmp_W=b.get('Pmax_W'),
            Pmp_difference_pct=100*(b.get('Pmax_W',np.nan)/a.get('Pmax_W',np.nan)-1),
            archived_cell_temperature_C=a.get('cell_temperature_C'),
            explicit_cell_temperature_C=b.get('cell_temperature_C')))
    pd.DataFrame(records).to_csv(runner.out/'preflight_explicit_25C_vs_archive.csv',index=False)
    # Transparent-band controls for every geometry at STC, with matching clean.
    controls=[State(name,shape,ori,h,1.) for name,shape,ori,h,tau in cases()]
    runner.batch(controls,workers,'Preflight: tau=1 zero additional loss')
    records=[]
    for state in sorted(set(controls),key=lambda s:s.key):
        clean=runner.results[State(state.module)]
        soiled=runner.results[state]
        loss=100*(1-soiled.get('Pmax_W',np.nan)/clean.get('Pmax_W',np.nan))
        records.append(dict(**asdict(state),loss_percent=loss,pass_zero_loss=bool(abs(loss)<1e-5)))
    control=pd.DataFrame(records)
    control.to_csv(runner.out/'preflight_tau1.csv',index=False)
    if not control.pass_zero_loss.all():
        raise RuntimeError('Transparent-band preflight failed')
    print('Preflight PASSED: archived reproduction and tau=1; true-25-C mismatch quantified separately.',flush=True)


def compile_results(runner):
    records=[]
    for name,shape,ori,h,tau in cases():
        ref=runner.results[State(name,shape,ori,h,tau)]
        clean_ref=runner.results[State(name)]
        stc=100*(1-ref.get('Pmax_W',np.nan)/clean_ref.get('Pmax_W',np.nan))
        for G in GS:
            a=runner.results[State(name,G=G)]
            b=runner.results[State(name,shape,ori,h,tau,G)]
            ok=a['status']=='ok' and b['status']=='ok'
            loss=100*(1-b.get('Pmax_W',np.nan)/a.get('Pmax_W',np.nan))
            row=dict(module=name,orientation=ori,shape=shape,band_height_mm=h,tau_I=tau,
                irradiance_W_m2=G,cell_temperature_C=25.,nominal_parameter_temperature_C=25.,
                clean_Pmp_W=a.get('Pmax_W',np.nan),soiled_Pmp_W=b.get('Pmax_W',np.nan),
                relative_power_loss_percent=loss,L_STC_percent=stc,Delta_L_pp=loss-stc,
                clean_Vmp_V=a.get('Vmp_V',np.nan),soiled_Vmp_V=b.get('Vmp_V',np.nan),
                clean_Imp_A=a.get('Imp_A',np.nan),soiled_Imp_A=b.get('Imp_A',np.nan),
                clean_bypass_MPP_mask=a.get('bypass_MPP_mask',np.nan),
                soiled_bypass_MPP_mask=b.get('bypass_MPP_mask',np.nan),
                STC_bypass_MPP_mask=ref.get('bypass_MPP_mask',np.nan),
                bypass_MPP_changed_from_STC=b.get('bypass_MPP_mask')!=ref.get('bypass_MPP_mask'),
                status='ok' if ok else 'failed')
            for k in (1,2,3):
                row[f'bypass_{k}_MPP_A']=b.get(f'bypass_{k}_MPP_A',np.nan)
                row[f'bypass_{k}_MPP_voltage_V']=b.get(f'bypass_{k}_MPP_voltage_V',np.nan)
            records.append(row)
    frame=pd.DataFrame(records)
    frame.to_csv(runner.out/'simulation_results.csv',index=False)
    groups=['module','orientation','shape','band_height_mm','tau_I']
    summary=[]
    for key,q in frame.groupby(groups):
        valid=q[np.isfinite(q.Delta_L_pp)]
        worst=valid.loc[valid.Delta_L_pp.abs().idxmax()] if len(valid) else None
        summary.append(dict(zip(groups,key),L_STC_percent=q.L_STC_percent.iloc[0],
            max_abs_Delta_L_pp=valid.Delta_L_pp.abs().max(),
            min_Delta_L_pp=valid.Delta_L_pp.min(),max_Delta_L_pp=valid.Delta_L_pp.max(),
            worst_irradiance_W_m2=worst.irradiance_W_m2 if worst is not None else np.nan,
            valid_conditions=len(valid),requested_conditions=len(GS),
            any_bypass_MPP_change=bool(q.bypass_MPP_changed_from_STC.any())))
    summary=pd.DataFrame(summary)
    summary.to_csv(runner.out/'configuration_summary.csv',index=False)
    frame[frame.bypass_MPP_changed_from_STC].to_csv(runner.out/'bypass_transitions.csv',index=False)
    return frame,summary


def figures(frame,out):
    for quantity,filename,ylabel in (
        ('relative_power_loss_percent','fig_loss_vs_irradiance','Relative power loss (%)'),
        ('Delta_L_pp','fig_deviation_from_STC','Difference from STC loss (percentage points)')):
        with plt.rc_context(STYLE):
            fig,axes=plt.subplots(4,4,figsize=(11,8),sharex=True,sharey=quantity=='relative_power_loss_percent')
            for col,name in enumerate(ORDER):
                for row,(shape,ori) in enumerate([(s,o) for s in ('rectangle','triangle') for o in ORIENTATIONS]):
                    ax=axes[row,col]
                    q=frame[(frame.module==name)&(frame['shape']==shape)&(frame.orientation==ori)&(frame.band_height_mm==80)]
                    for tau in TAUS:
                        g=q[np.isclose(q.tau_I,tau)].sort_values('irradiance_W_m2')
                        ax.plot(g.irradiance_W_m2,g[quantity],color=COLORS[tau],marker='o',ms=3,lw=1.1)
                    if row==0: ax.set_title(module_label(name))
                    if col==0: ax.set_ylabel(f'{shape.capitalize()}\n{ori}',fontsize=8)
                    if row==3: ax.set_xlabel('Irradiance (W m$^{-2}$)')
                    ax.set_xticks([200,600,1000])
                    if quantity=='Delta_L_pp':
                        ax.axhline(0,color='.5',lw=.5)
                        for level in (-2,-1,1,2): ax.axhline(level,color='.8',lw=.5,ls=':')
                    ax.text(.03,.96,'abcdefghijklmnop'[row*4+col],transform=ax.transAxes,va='top',fontweight='bold')
            fig.legend(handles=[Line2D([],[],color=COLORS[t],marker='o',ms=3,label=rf'$\tau_I={t:g}$') for t in TAUS],
                       loc='upper center',ncol=3,bbox_to_anchor=(.5,.965))
            fig.suptitle('80 mm band; explicit 25 °C isothermal sensitivity (unchanged fitted electrical coefficients)',y=.995,fontsize=9)
            fig.supylabel(ylabel,x=.005,fontsize=9)
            fig.tight_layout(rect=(.035,0,1,.92))
            save(fig,out/filename)
    with plt.rc_context(STYLE):
        fig,axes=plt.subplots(2,4,figsize=(9.5,4.5),sharex=True,sharey='row')
        for col,name in enumerate(ORDER):
            for ori,ls in (('portrait','-'),('landscape','--')):
                q=frame[(frame.module==name)&(frame.band_height_mm==60)&(frame.orientation==ori)].sort_values('irradiance_W_m2')
                for ax,key in zip(axes[:,col],('relative_power_loss_percent','Delta_L_pp')):
                    ax.plot(q.irradiance_W_m2,q[key],color='#0072B2',ls=ls,marker='o',ms=3)
            axes[0,col].set_title(module_label(name))
            axes[1,col].axhline(0,color='.5',lw=.5)
            axes[1,col].set_xlabel('Irradiance (W m$^{-2}$)')
            axes[1,col].set_xticks([200,600,1000])
        axes[0,0].set_ylabel('Relative power loss (%)')
        axes[1,0].set_ylabel('Difference from STC (pp)')
        fig.legend(handles=[Line2D([],[],color='#0072B2',ls=ls,label=ori) for ori,ls in (('portrait','-'),('landscape','--'))],
                   loc='upper center',ncol=2,bbox_to_anchor=(.5,.955))
        fig.suptitle('Annual reference geometry: rectangular band, h = 60 mm, τI = 0.44; 25 °C',y=.995,fontsize=9)
        fig.tight_layout(rect=(0,0,1,.89))
        save(fig,out/'fig_annual_reference_geometry')


def grid_convergence(runner,frame,workers):
    """Independent refinement of worst loss and bypass-change configurations."""
    chosen=[frame.loc[frame.Delta_L_pp.abs().idxmax()]]
    annual=frame[frame.band_height_mm==60]
    chosen.append(annual.loc[annual.Delta_L_pp.abs().idxmax()])
    chosen.extend(row for _,row in frame[frame.bypass_MPP_changed_from_STC].iterrows())
    states=[]
    for row in chosen:
        for G in (int(row.irradiance_W_m2),1000):
            states += [State(row.module,G=G,track='explicit_25C_refined'),
                       State(row.module,row['shape'],row.orientation,row.band_height_mm,row.tau_I,G,'explicit_25C_refined')]
    runner.batch(states,workers,'MPP grid convergence: 0.01 V')
    records=[]
    for row in chosen:
        losses=[]
        for G in (int(row.irradiance_W_m2),1000):
            clean=runner.results[State(row.module,G=G,track='explicit_25C_refined')]
            shade=runner.results[State(row.module,row['shape'],row.orientation,row.band_height_mm,row.tau_I,G,'explicit_25C_refined')]
            losses.append(100*(1-shade.get('Pmax_W',np.nan)/clean.get('Pmax_W',np.nan)))
        refined=losses[0]-losses[1]
        records.append(dict(module=row.module,shape=row['shape'],orientation=row.orientation,
            band_height_mm=row.band_height_mm,tau_I=row.tau_I,irradiance_W_m2=row.irradiance_W_m2,
            Delta_L_002V_pp=row.Delta_L_pp,Delta_L_001V_pp=refined,
            refinement_difference_pp=refined-row.Delta_L_pp,
            pass_grid_convergence=bool(abs(refined-row.Delta_L_pp)<.02)))
    convergence=pd.DataFrame(records)
    convergence.to_csv(runner.out/'grid_convergence.csv',index=False)
    if not convergence.pass_grid_convergence.all():
        raise RuntimeError('MPP grid convergence failed; inspect grid_convergence.csv')


def temperature_audit(out):
    records=[]
    for name in ORDER:
        m=modules.get(name)
        text=spice.read_netlist(m.netlist_path)
        records.append(dict(module=name,fitted_nominal_parameter_temperature_C=m.t_ref,
            runtime_default_temperature_C=27,temperature_setter=m.temperature,
            temperature_sweep_supported=False,
            reason='t changes both Iph and diode TNOM; no operating .temp. Varying t is not varying a physical cell temperature.',
            bishop_branch_present=bool(re.search(r'^B1\b',text,re.M)),
            unused_voltage_temperature_coefficient=bool(re.search(r'\.param\s+kv\b',text,re.I)),
            explicit_25C_control='.param t=25; .temp 25; diode nominal temperature remains fixed at 25 throughout irradiance sweep',
            fixed_parameters='IS/IS1, n1, Rs, Rp, Bishop a/b/Vb, bypass model unchanged; no refitting'))
    (out/'temperature_model_audit.json').write_text(json.dumps(records,indent=2)+'\n')
    pd.DataFrame([dict(module=n,shape=shape,orientation=ori,band_height_mm=h,tau_I=tau,
        irradiance_W_m2=g,requested_temperature_C=t,status='not_run_unsupported',
        reason='Existing t/TNOM coupling; no validated independent physical temperature control')
        for n,shape,ori,h,tau in cases() for g in (400,1000) for t in (15,45,65)]).to_csv(out/'temperature_conditions_not_run.csv',index=False)


def report(frame,summary,out,protection):
    good=frame[frame.status=='ok']
    worst=good.loc[good.Delta_L_pp.abs().idxmax()]
    nonstc=good[good.irradiance_W_m2!=1000]
    ranges=nonstc.groupby('module').Delta_L_pp.agg(['min','max'])
    ranges['max_abs']=nonstc.groupby('module').Delta_L_pp.apply(lambda s:s.abs().max())
    ranges.to_csv(out/'module_summary.csv')
    legacy=pd.read_csv(out/'preflight_archived_reference.csv')
    explicit=pd.read_csv(out/'preflight_explicit_25C_vs_archive.csv')
    convergence=pd.read_csv(out/'grid_convergence.csv')
    lines=['# Irradiance sensitivity of the STC loss approximation','','## Scope and temperature limitation','',
        f'{len(frame)} requested loss conditions: 56 band configurations × five irradiances. '
        'Four modules; 80 mm rectangular/triangular bands in both orientations with exact '
        'tau_I=0.71, 0.44, 0; plus 60 mm rectangular bands at exact tau_I=0.44 in both orientations. '
        'These requested transmittances are not silently replaced by the campaign values 0.706/0.442.','',
        'The original pipeline writes t into a parameter that controls both photocurrent and diode TNOM. '
        'It does not set .temp. LTspice logs confirm operating temperature 27 °C in the archived-reproduction '
        'track, with module-specific nominal parameters 24.8 °C (Yingli), 25.4 °C (SunPower), '
        '24.2 °C (JA Solar), 25 °C (LONGi). Therefore the original reference is not strictly a '
        '1000 W/m², actual-25-°C simulation. This cannot be represented as an exact STC reproduction.','',
        'First, archived curves are reproduced using the original settings/grid. Separately, the '
        'irradiance sensitivity uses explicit .temp 25 AND t=25 for all G, without changing fitted '
        'IS/IS1, n, Rs, Rp, breakdown or bypass coefficients. This creates a consistently isothermal '
        '25 °C reference, but its absolute powers can differ from the archive; '
        'preflight_explicit_25C_vs_archive.csv quantifies the difference. No archived result is overwritten. '
        'The requested L_STC is recomputed within this explicit-25-C track, never borrowed from another track.','',
        'The multi-temperature sweep is **not performed**: changing t would change the nominal '
        'diode calibration while changing photocurrent, instead of independently moving operating '
        'temperature with saturation-current and thermal-voltage responses referenced to a fixed TNOM. '
        'The native diode has a temperature law, but the current wrapper does not separate these controls. '
        'No arbitrary thermal law, datasheet power correction, or fitted temperature coefficient is added. '
        'Constants kv, q and k are declared but do not implement a separate cell-voltage temperature model.','',
        '## Model implementation','',
        '- Existing soilspice.spice.read_netlist, patch_cells, set_temperature, set_param, '
        'run_ltspice and iv_parameters are reused. A native binary RAW reader from the previous '
        'LONGi diagnostic supplies the saved traces when PyLTSpice is unavailable.',
        '- Cell photocurrent is [Isc + Isc kc (t−25)] (Vi/1000) [1−f_cov(1−tau_I)]. '
        'Vi is set to each requested irradiance; geometric coverage uses the original '
        'geometry.fcov_grid/cell_fcov and the same four-decimal netlist precision.',
        '- SunPower and LONGi contain the explicit Bishop branch '
        'I_B=(Vj/Rp) a pwr(1−Vj/Vb,−b). Yingli and JA Solar do not contain this branch '
        'in their current templates: reverse current is predominantly the ordinary shunt '
        'and any native diode leakage. All use their existing fixed shunt resistance.',
        '- Module current sign is the configured +I(V1), checked against terminal voltage; '
        'MPP is the global maximum of V I over the entire first quadrant, including all local peaks. '
        'The sensitivity sweep is extended 15% and sampled at 0.02 V; original grids are retained '
        'for archived reproduction. Validity requires positive power, interior MPP and a Voc crossing.',
        '- Bypass activity is measured directly from the three module bypass diode currents '
        'and terminal voltages. An MPP diode is active when its forward current exceeds '
        'max(0.1 mA, 1% of module MPP current) and its voltage is forward. The mask identifies '
        'which physical bypass diode is active; parallel branches are already represented '
        'in the original full-module netlist. These are MPP states, not claims about the entire I–V curve.','',
        '## Validation','',
        'preflight_archived_reference.csv checks all four clean modules and 48 campaign '
        'shading states (two shapes, two orientations, three exact available campaign transmittances). '
        'The unavailable 0.442 landscape campaign states use the saved 0.536 state for '
        'reproduction only; the main sensitivity always uses the requested exact 0.44. '
        'The original SunPower sweep can stop short of Voc, although its MPP is interior. '
        'That truncated sweep is preserved only for archived reproduction; every new '
        'sensitivity curve must bracket Voc. Archived powers are checked '
        'with a 0.02 W reproduction tolerance. The transparent-band control checks every '
        'unique geometry at tau_I=1 with loss tolerance 1e-5 percentage points. '
        'Every run and any failure is retained in run_manifest.csv and failed_simulations.csv. '
        'The worst overall, worst annual-reference and bypass-change cases are recomputed '
        'with a 0.01 V grid; grid_convergence.csv reports the resulting deviation changes.','',
        f'Maximum archived power reproduction error: {legacy.difference_W.abs().max():.3g} W. '
        f'True-25-C powers differ from archived powers by '
        f'{explicit.Pmp_difference_pct.min():.6f} to {explicit.Pmp_difference_pct.max():.6f}%, '
        'including operating/nominal temperature settings and numerical grid refinement. '
        'This discrepancy is a reported finding, not a silently passed true-STC equivalence. '
        'Refining the grid from 0.02 V to 0.01 V for the global worst case, '
        'the annual-geometry worst case and the bypass-transition case changes ΔL by at most '
        f'{convergence.refinement_difference_pp.abs().max():.6f} pp; see grid_convergence.csv.','',
        f'Result rows: {len(good)} successful / {len(frame)} requested; '
        f'{int((frame.status!="ok").sum())} failed loss conditions.','',
        '## Numerical findings','',
        '| Module | Minimum ΔL (pp) | Maximum ΔL (pp) | Maximum absolute ΔL (pp) |',
        '|---|---:|---:|---:|']
    for name,r in ranges.iterrows():
        lines.append(f'| {module_label(name)} | {r["min"]:.6f} | {r["max"]:.6f} | {r.max_abs:.6f} |')
    lines += ['',f'Largest absolute deviation: **{abs(worst.Delta_L_pp):.6f} pp**, '
        f'{module_label(worst.module)}, {worst["shape"]}, {worst.orientation}, '
        f'h={worst.band_height_mm:g} mm, tau_I={worst.tau_I:g}, G={worst.irradiance_W_m2:g} W/m². '
        f'L_STC={worst.L_STC_percent:.6f}%, L(G)={worst.relative_power_loss_percent:.6f}%.','',
        f'Among {len(nonstc)} non-STC conditions, '
        f'{int((nonstc.Delta_L_pp.abs()<=1).sum())} are within 1 pp and '
        f'{int((nonstc.Delta_L_pp.abs()<=2).sum())} are within 2 pp. '
        f'{int((nonstc.Delta_L_pp.abs()>2).sum())} exceed 2 pp.','',
        f'{int(good.bypass_MPP_changed_from_STC.sum())} conditions have a different bypass-diode '
        'MPP mask from their own STC state; see bypass_transitions.csv.','']
    for keys,label in [(['shape'],'Shape'),(['orientation'],'Orientation'),(['tau_I'],'Transmittance'),(['band_height_mm'],'Band height')]:
        grouped=nonstc.groupby(keys[0]).Delta_L_pp.apply(lambda s:s.abs().max())
        lines.append(label+' maximum absolute deviations: '+', '.join(f'{key}: {value:.6f} pp' for key,value in grouped.items())+'.')
    annual=summary[summary.band_height_mm==60]
    lines += ['',f'The 60 mm annual-reference geometries have a worst deviation of '
              f'{annual.max_abs_Delta_L_pp.max():.6f} pp across the tested irradiances.','',
        '## Scientific interpretation and limitations','',
        'Changes with irradiance reflect the existing nonlinear circuit: photocurrent decreases '
        'with G, whereas fixed shunt conductance, series resistance, diode voltage and reverse '
        'current do not scale proportionally. MPP branches and bypass states can therefore change. '
        'Differences between modules cannot be attributed solely to architecture: fitted cell '
        'parameters and breakdown models also differ.',
        'Geometry and transmittance can change which cells are limiting and the preferred MPP. '
        'Use configuration_summary.csv rather than extrapolating one geometry to all cases.',
        'The one MPP bypass transition is SunPower, triangular opaque band, 80 mm, '
        'landscape: bypass diode D1 conducts at the STC MPP but not at the 200 W/m² MPP. '
        'The largest global exception occurs instead in portrait without an MPP '
        'bypass-mask change; reverse/shunt response and competing MPPs also matter.',
        'This is a numerical sensitivity test, **not independent experimental validation**. '
        'The discrepancy between the archived nominal reference and true 25 °C, the insufficient '
        'temperature interface, constant Rsh and unvalidated low-irradiance/reverse-bias laws '
        'limit inference. Small errors for some configurations do not validate the annual '
        'STC approximation; large exceptions are explicitly retained.',
        'No annual energy-error estimate is made. Raw weather data exist, but a defensible '
        'hourly estimate would require an independently controlled temperature model and the '
        'time-varying band transmittance, including irradiance below 200 W/m² and above '
        '1000 W/m². The fixed-tau grid here must not be substituted for that integration.','',
        'The annual model also interpolates an archived transmittance table, whereas '
        'this test directly simulates tau_I=0.44. The 60 mm results assess that '
        'representative geometry, not the full error of the hourly annual pipeline.','',
        '## Reproduction','',
        '    python chapters/stc_sensitivity.py --workers 4',
        'Netlists, RAW outputs, logs, and per-run metadata are retained in runs/. '
        'Successful runs with identical netlist hashes are reused. Windows source-placeholder '
        'and .lib path adaptations are confined to copies. The source/input/manuscript '
        'hash check is in verification.json.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    si=('Using the unchanged fitted cell-level circuits, we evaluated the relative module '
        'power loss at irradiances of 200–1000 W m−2 for 56 bottom-edge band configurations. '
        'Clean and shaded maximum powers were recomputed at the same irradiance and explicit '
        'operating temperature of 25 °C, with a global 0.02 V voltage sweep. '
        f'The maximum absolute departure from the corresponding STC relative loss was '
        f'{abs(worst.Delta_L_pp):.2f} percentage points; '
        f'{int((nonstc.Delta_L_pp.abs()>2).sum())} of {len(nonstc)} non-STC cases exceeded '
        f'2 percentage points. The largest exception was {module_label(worst.module)}, '
        f'{worst["shape"]} band, {worst.orientation}, h={worst.band_height_mm:g} mm, '
        f'tau_I={worst.tau_I:g}, at {worst.irradiance_W_m2:g} W m−2. '
        f'For the 60 mm rectangular geometry at tau_I=0.44, all eight module/orientation '
        f'configurations remained within {annual.max_abs_Delta_L_pp.max():.2f} percentage '
        f'points over the tested irradiances. A change in the MPP bypass-diode state '
        f'occurred in {int(good.bypass_MPP_changed_from_STC.sum())} condition. '
        'The archived calculations were reproduced separately using '
        'their original nominal-temperature settings; their native operating temperature '
        'was 27 °C. Because the existing parameter t also sets diode TNOM, a physically '
        'consistent independent temperature sweep was not supported and was not imposed. '
        'These calculations assess numerical sensitivity within the existing model and '
        'do not independently validate either its low-irradiance response or the annual '
        'STC-loss approximation. See report.md for per-module, geometry and bypass results.')
    (out/'SI_text.md').write_text(si+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=config.OUTPUT_DIR/'stc_sensitivity')
    parser.add_argument('--exe',type=Path,default=Path(os.environ.get('LOCALAPPDATA',''))/'Programs/ADI/LTspice/LTspice.exe')
    parser.add_argument('--library',type=Path,default=Path(os.environ.get('LOCALAPPDATA',''))/'LTspice/lib/cmp/standard.dio')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    protected={p:digest(p) for directory in (config.INPUT_DIR,config.OUTPUT_DIR,config.REPO_DIR/'soilspice')
               for p in directory.rglob('*') if p.is_file() and out not in p.parents
               and '__pycache__' not in p.parts}
    config.LTSPICE_EXE=args.exe.resolve()  # runtime only, never write config.py
    runner=Runner(out,args.library)
    temperature_audit(out)
    validate(runner,args.workers)
    jobs=[State(name,G=G) for name in ORDER for G in GS]
    jobs += [State(name,shape,ori,h,tau,G) for name,shape,ori,h,tau in cases() for G in GS]
    runner.batch(jobs,args.workers,'Irradiance sensitivity at explicit 25 C')
    frame,summary=compile_results(runner)
    grid_convergence(runner,frame,args.workers)
    figures(frame,out)
    report(frame,summary,out,protected)
    modified=[str(p) for p,h in protected.items() if digest(p)!=h]
    verification=dict(protected_file_count=len(protected),modified_protected_files=modified,
        archived_reproduction_passed=True,tau1_checks_passed=True,
        requested_loss_conditions=len(frame),successful_loss_conditions=int((frame.status=='ok').sum()),
        successful_ltspice_runs=sum(r['status']=='ok' for r in runner.results.values()),
        failed_ltspice_runs=sum(r['status']!='ok' for r in runner.results.values()),
        electrical_templates_sha256={str(p.relative_to(config.REPO_DIR)):h for p,h in protected.items() if p.parent==config.NETLIST_DIR},
        LTspice_executable=str(config.LTSPICE_EXE),diode_library_sha256=digest(args.library),
        temperature_sweep_supported=False,
        explicit_STC_differs_from_archived_nominal_reference=True,
        protected_files_sha256={str(p.relative_to(config.REPO_DIR)):h for p,h in protected.items()})
    (out/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    assert not modified,modified
    print(summary.sort_values('max_abs_Delta_L_pp',ascending=False).head(10).to_string(index=False),flush=True)
    print(f'Outputs: {out}; protected files unchanged: {len(protected)}',flush=True)


if __name__=='__main__':
    main()
