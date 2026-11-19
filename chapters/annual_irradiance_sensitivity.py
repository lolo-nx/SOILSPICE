"""Energy-weighted irradiance-only correction, PH15mm rain-only reference.

Run with pvlib and pandas 2.x installed:
    python chapters/annual_irradiance_sensitivity.py --workers 4

Uses saved band trajectories and the existing hourly clean-production model.
Does not change manuscript results, cleaning policies or physical-model files.
The delta from a common 25-C STC baseline is applied to the ORIGINAL archived
loss LUT, so nominal-temperature and interpolation offsets are not mistaken
for irradiance sensitivity. No electrical temperature sensitivity is claimed.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from soilspice import config, energy, modules
from soilspice.style import ORDER, ORIENTATIONS, SITE_COLOR, STYLE, save
from chapters import stc_sensitivity as electrical

TAUS=(0.,.1,.2,.3,.4,.44,.5,.6,.7,.71,.8,.9,.95,.98,.99,1.)
BASE_G=(1,5,20,50,100,150,200,300,400,500,600,700,800,900,1000)


def run_ltspice_with_fallback(net_path):
    """Same circuit; retain a timed-out normal-solver log and retry Alternate.

    This is a numerical-solver selection only, not a parameter or tolerance
    change. subprocess.run kills/waits for its own child on timeout.
    """
    command=[str(config.LTSPICE_EXE),'-b',str(net_path.resolve())]
    metadata=dict(normal_timeout_seconds=20,alternate_used=False)
    flags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0
    try:
        result=subprocess.run(command,capture_output=True,timeout=20,creationflags=flags)
    except subprocess.TimeoutExpired:
        log=net_path.with_suffix('.log')
        if log.exists(): shutil.copyfile(log,net_path.with_suffix('.normal_timeout.log'))
        metadata['alternate_used']=True
        result=subprocess.run([str(config.LTSPICE_EXE),'-alt','-b',str(net_path.resolve())],
                              capture_output=True,timeout=120,creationflags=flags)
    net_path.with_suffix('.solver.json').write_text(json.dumps(metadata,indent=2)+'\n')
    raw=net_path.with_suffix('.raw')
    if result.returncode or not raw.exists():
        raise RuntimeError(f'LTspice failed: {net_path.name}, return code {result.returncode}')
    return raw


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def effective_weather(hours):
    """Recover actual model G and Tc from its two linear-temperature yields.

    Existing model: Y_m = (G/1000)[1+gamma_m(Tc-25)]. Two different
    gamma values uniquely recover both inputs without duplicating/changing
    the transposition, IAM, or Faiman calculation in energy.pv_yield.
    Night-time Tc is undefined by zero power and is not used in weighting.
    """
    ga,gb=energy.GAMMA_T['JASOLAR'],energy.GAMMA_T['LONGI']
    ya,yb=hours.Y_JASOLAR.to_numpy(),hours.Y_LONGI.to_numpy()
    thermal=(ya-yb)/(ga-gb)
    g_kw=ya-ga*thermal
    if np.min(g_kw)<-1e-9: raise ValueError('Negative reconstructed irradiance')
    g_kw=np.maximum(0,g_kw)
    tc=25+np.divide(thermal,g_kw,out=np.full_like(g_kw,np.nan),where=g_kw>1e-10)
    for name,gamma in energy.GAMMA_T.items():
        np.testing.assert_allclose(g_kw+gamma*thermal,hours[f'Y_{name}'],rtol=1e-10,atol=1e-12)
    return 1000*g_kw,tc


def align_saved_series(hours,series):
    """Map each actual source hour to its saved event/midnight interval."""
    series=series.sort_values('date').copy()
    starts=pd.DatetimeIndex(pd.to_datetime(series.date,utc=True)).as_unit('ns')
    ends=pd.DatetimeIndex(pd.to_datetime(series.end_time,utc=True)).as_unit('ns')
    index=hours.index.as_unit('ns')
    left=index.get_indexer(starts)
    right=index.searchsorted(ends)
    assert (left>=0).all() and left[0]==0 and right[-1]==len(hours)
    np.testing.assert_array_equal(left[1:],right[:-1])
    counts=right-left
    assert (counts>0).all()
    np.testing.assert_allclose(counts/24,series.duration_days,rtol=1e-12,atol=1e-12)
    interval=np.repeat(np.arange(len(series)),counts)
    age_hours=np.arange(len(hours))-left[interval]+1
    return dict(interval=interval,age_hours=age_hours,left=left,right=right)


def prepare_sites(scenario,out):
    trajectory=pd.read_csv(scenario/'data/growth_series_reference.csv')
    settings_path=scenario/'data/model_settings.json'
    try:
        settings_text=settings_path.read_text(encoding='utf-8-sig')
    except UnicodeDecodeError:
        settings_text=settings_path.read_text(encoding='cp1252')
    settings=json.loads(settings_text)
    assert settings['washing_event_total_mm']==15
    source=pd.read_csv(scenario/'data/energy_maps.csv')
    raw_runs=pd.read_csv(scenario/'data/model_runs.csv')
    assert energy.GROWTH_REFERENCE==('JASOLAR','landscape',60.0,.1,0.0)
    sites={}
    checks=[]
    bins=[]
    support=[]
    for site in energy.SITES:
        hours=energy.site_hours(site)
        hours=hours[(hours.index.year>=energy.YEARS[0])&(hours.index.year<=energy.YEARS[1])].copy()
        G,Tc=effective_weather(hours)
        series=trajectory[trajectory.site==site].sort_values('date')
        mapping=align_saved_series(hours,series)
        tau_ref=series.tau.to_numpy()[mapping['interval']]
        clean_ref=hours.Y_JASOLAR.to_numpy()
        # All saved trajectories use JA Solar landscape, h=60mm, f=0.1.
        predicted=np.add.reduceat(clean_ref,mapping['left'])
        np.testing.assert_allclose(predicted,series.Y,rtol=2e-6,atol=1e-9)
        annual_years=series.duration_days.sum()/365.25
        D=energy.deposition(site)
        mu=0.
        mu_start=[]
        for row in series.itertuples(index=False):
            mu_start.append(mu)
            mu+=D*row.duration_days
            if row.rain>=energy.LIGHT: mu=0.
        uniform=np.exp(-energy.ALPHA*(np.asarray(mu_start)[mapping['interval']]+D*mapping['age_hours']/24))
        uniform_interval=np.add.reduceat(clean_ref*uniform,mapping['left'])
        observed=series.Y.to_numpy()*(1-series.L_u.to_numpy())
        np.testing.assert_allclose(uniform_interval,observed,rtol=2e-6,atol=1e-9)
        sites[site]=dict(hours=hours,G=G,Tc=Tc,uniform=uniform,tau_ref=tau_ref,
            annual_years=annual_years,series=series,mapping=mapping)
        for name in ORDER:
            y=hours[f'Y_{name}'].to_numpy()
            support.append(dict(site=site,module=name,hour_count=len(hours),
                irradiance_min_positive_W_m2=float(G[G>0].min()),irradiance_max_W_m2=float(G.max()),
                energy_weighted_cell_temperature_C=float(np.nansum(y*Tc)/y.sum()),
                clean_energy_below_20_pct=float(100*y[G<20].sum()/y.sum()),
                clean_energy_below_200_pct=float(100*y[G<200].sum()/y.sum()),
                clean_energy_below_400_pct=float(100*y[G<400].sum()/y.sum()),
                clean_energy_above_1000_pct=float(100*y[G>1000].sum()/y.sum())))
            for lo,hi in zip((0,1,20,50,100,200,400,600,800,1000,1200),
                             (1,20,50,100,200,400,600,800,1000,1200,np.inf)):
                mask=(G>=lo)&(G<hi)
                bins.append(dict(site=site,module=name,G_lower_W_m2=lo,G_upper_W_m2=hi,
                    hour_count=int(mask.sum()),clean_energy_kwh_kwp=float(y[mask].sum()),
                    fraction_clean_energy=float(y[mask].sum()/y.sum())))
            for orientation in ORIENTATIONS:
                ratio=energy.slope_length(name,orientation)/energy.slope_length('JASOLAR','landscape')
                tau=tau_ref**ratio
                tau_lookup=np.rint(tau*200)/200
                curve=energy.curve(name,orientation,60,'rectangle')
                stc=np.interp(tau_lookup,np.r_[energy.spice_table('rectangle')[0],1.],np.r_[curve,0.])
                predicted_band=100*np.sum(y*uniform*stc)/np.sum(y)
                expected=source[(source.site==site)&(source.module==name)&(source.orientation==orientation)
                                &np.isclose(source.h,60)&np.isclose(source.f,.1)&(source['shape']=='rectangle')].iloc[0]
                old_band=100*(expected.loss-expected.loss_u)
                saved_y=raw_runs[(raw_runs.site==site)&(raw_runs.module==name)
                    &(raw_runs.orientation==orientation)&(raw_runs.market=='ground')&(raw_runs['mode']=='none')
                    &np.isclose(raw_runs.h,60)&np.isclose(raw_runs.f,.1)&np.isclose(raw_runs.margin,0)].Y.iloc[0]
                np.testing.assert_allclose(y.sum()/annual_years,saved_y,rtol=2e-6,atol=1e-7)
                checks.append(dict(site=site,module=name,orientation=orientation,
                    recomputed_band_loss_pct=predicted_band,saved_band_loss_pct=old_band,
                    discrepancy_pp=predicted_band-old_band,
                    recomputed_clean_energy_kwh_kwp_year=float(y.sum()/annual_years),
                    saved_clean_energy_kwh_kwp_year=float(saved_y)))
                assert abs(predicted_band-old_band)<2e-5,(site,name,orientation,predicted_band,old_band)
        print(f'Weather/trajectory preflight {site}: {len(hours)} hours; Gmax={G.max():.1f} W/m2',flush=True)
    pd.DataFrame(checks).to_csv(out/'baseline_reproduction.csv',index=False)
    pd.DataFrame(bins).to_csv(out/'irradiance_energy_distribution.csv',index=False)
    pd.DataFrame(support).to_csv(out/'weather_energy_support.csv',index=False)
    return sites


def grid_states(Gs):
    return ([electrical.State(name,G=G) for name in ORDER for G in Gs]+
            [electrical.State(name,'rectangle',orientation,60,tau,G)
             for name in ORDER for orientation in ORIENTATIONS for tau in TAUS for G in Gs])


def verify_alternate_solver(runner):
    """Compare Normal and Alternate on a paired SunPower clean/shaded state."""
    records=[]
    directory=runner.out/'solver_checks'
    directory.mkdir(exist_ok=True)
    for state in (electrical.State('SUNPOWER'),
                  electrical.State('SUNPOWER','rectangle','portrait',60,.44,1000)):
        text,diodes,nodes=electrical.netlist(state,runner.library)
        path=directory/(state.key+'_alternate.net')
        path.write_text(text,encoding='utf-8')
        result=subprocess.run([str(config.LTSPICE_EXE),'-alt','-b',str(path.resolve())],
            capture_output=True,timeout=120,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
        assert result.returncode==0
        alternate=electrical.measure(path.with_suffix('.raw'),diodes,nodes)
        normal=runner.results[state]
        difference=alternate['Pmax_W']-normal['Pmax_W']
        records.append(dict(state=state.key,normal_Pmp_W=normal['Pmax_W'],alternate_Pmp_W=alternate['Pmax_W'],
                            difference_W=difference,pass_solver_comparison=bool(abs(difference)<.02)))
    check=pd.DataFrame(records)
    check.to_csv(runner.out/'alternate_solver_comparison.csv',index=False)
    assert check.pass_solver_comparison.all()


def native_grid(runner,Gs):
    records=[]
    lookup={}
    for name in ORDER:
        for orientation in ORIENTATIONS:
            values=np.zeros((len(Gs),len(TAUS)))
            for i,G in enumerate(Gs):
                clean=runner.results[electrical.State(name,G=G)]
                for j,tau in enumerate(TAUS):
                    shade=runner.results[electrical.State(name,'rectangle',orientation,60,tau,G)]
                    if clean['status']!='ok' or shade['status']!='ok':
                        raise RuntimeError('An annual LUT simulation failed; inspect failed_simulations.csv')
                    loss=1-shade['Pmax_W']/clean['Pmax_W']
                    if not -.0001<=loss<=1.0001: raise ValueError((name,G,tau,loss))
                    values[i,j]=loss
                    records.append(dict(module=name,orientation=orientation,band_height_mm=60,
                        tau_I=tau,irradiance_W_m2=G,cell_temperature_C=25,
                        clean_Pmp_W=clean['Pmax_W'],soiled_Pmp_W=shade['Pmax_W'],relative_band_loss_pct=100*loss))
            stc=values[Gs.index(1000)]
            lookup[name,orientation]=RegularGridInterpolator((Gs,TAUS),values-stc[None,:],
                                                            bounds_error=True)
    pd.DataFrame(records).to_csv(runner.out/'irradiance_tau_grid.csv',index=False)
    return lookup


def aggregate_sites(sites,lookup,Gs,out):
    records=[]
    yearly=[]
    for site,weather in sites.items():
        h,G,U=weather['hours'],weather['G'],weather['uniform']
        for name in ORDER:
            y=h[f'Y_{name}'].to_numpy()
            for orientation in ORIENTATIONS:
                exponent=energy.slope_length(name,orientation)/energy.slope_length('JASOLAR','landscape')
                tau=np.rint(weather['tau_ref']**exponent*200)/200
                curve=energy.curve(name,orientation,60,'rectangle')
                baseline=np.interp(tau,np.r_[energy.spice_table('rectangle')[0],1.],np.r_[curve,0.])
                for G_basis,label in ((G*U,'after_uniform_layer'),(G,'weather_effective')):
                    covered=(G_basis>=Gs[0])&(G_basis<=Gs[-1])
                    delta=np.zeros_like(y)
                    delta[covered]=lookup[name,orientation](np.column_stack([G_basis[covered],tau[covered]]))
                    raw=baseline+delta
                    corrected=np.clip(raw,0,1)
                    weights=y*U
                    original_energy=weights*baseline
                    corrected_energy=weights*corrected
                    # Keep STC contribution outside the solved domain. Also
                    # supply a conservative physical [0,100% band loss] bound
                    # for those hours, rather than silently extrapolating.
                    lower=np.where(covered,corrected_energy,0)
                    upper=np.where(covered,corrected_energy,weights)
                    row=dict(site=site,module=name,orientation=orientation,shape='rectangle',h_mm=60,
                        f=.1,manual_cleaning='none',irradiance_basis=label,
                        clean_energy_kwh_kwp_year=float(y.sum()/weather['annual_years']),
                        original_band_loss_pct=100*original_energy.sum()/y.sum(),
                        corrected_band_loss_pct=100*corrected_energy.sum()/y.sum(),
                        annual_difference_pp=100*(corrected_energy-original_energy).sum()/y.sum(),
                        annual_difference_kwh_kwp_year=float((corrected_energy-original_energy).sum()/weather['annual_years']),
                        uniform_loss_pct=100*np.sum(y*(1-U))/y.sum(),
                        supported_hour_fraction=float(covered.mean()),
                        supported_clean_energy_fraction=float(y[covered].sum()/y.sum()),
                        supported_after_uniform_energy_fraction=float(weights[covered].sum()/weights.sum()),
                        corrected_lower_bound_pct=100*lower.sum()/y.sum(),
                        corrected_upper_bound_pct=100*upper.sum()/y.sum(),
                        clipped_hour_count=int(((raw<0)|(raw>1)).sum()),
                        clipping_effect_pp=100*np.sum(weights*(corrected-raw))/y.sum())
                    row['relative_change_band_loss_pct']=100*row['annual_difference_pp']/row['original_band_loss_pct']
                    row['original_total_loss_pct']=row['uniform_loss_pct']+row['original_band_loss_pct']
                    row['corrected_total_loss_pct']=row['uniform_loss_pct']+row['corrected_band_loss_pct']
                    records.append(row)
                    for year in sorted(h.index.year.unique()):
                        mask=h.index.year==year
                        yearly.append(dict(site=site,module=name,orientation=orientation,year=int(year),irradiance_basis=label,
                            clean_energy_kwh_kwp=float(y[mask].sum()),
                            original_band_loss_pct=100*original_energy[mask].sum()/y[mask].sum(),
                            corrected_band_loss_pct=100*corrected_energy[mask].sum()/y[mask].sum(),
                            annual_difference_pp=100*(corrected_energy[mask]-original_energy[mask]).sum()/y[mask].sum()))
    result=pd.DataFrame(records)
    result.to_csv(out/'annual_results_per_configuration.csv',index=False)
    pd.DataFrame(yearly).to_csv(out/'annual_results_by_year.csv',index=False)
    summary=result.groupby(['site','irradiance_basis']).agg(
        original_median_loss_pct=('original_band_loss_pct','median'),
        corrected_median_loss_pct=('corrected_band_loss_pct','median'),
        median_paired_difference_pp=('annual_difference_pp','median'),
        min_difference_pp=('annual_difference_pp','min'),max_difference_pp=('annual_difference_pp','max'),
        median_relative_change_pct=('relative_change_band_loss_pct','median'),
        min_supported_clean_energy_fraction=('supported_clean_energy_fraction','min')).reset_index()
    summary['difference_of_site_medians_pp']=summary.corrected_median_loss_pct-summary.original_median_loss_pct
    summary.to_csv(out/'annual_site_summary.csv',index=False)
    return result,summary


def interpolation_resolution_check(sites,result,Gs,out):
    """Coarsen existing samples as a diagnostic, without extra LTspice runs.

    This is a resolution-sensitivity check, not a guaranteed error bound on
    the full grid. The main result always uses all available samples.
    """
    coarse_g=[g for g in (1,20,100,400,800,1000,1200) if g in Gs]
    coarse_tau=[0.,.2,.4,.6,.8,.95,1.]
    grid=pd.read_csv(out/'irradiance_tau_grid.csv')
    lookup={}
    for name in ORDER:
        for orientation in ORIENTATIONS:
            q=grid[(grid.module==name)&(grid.orientation==orientation)]
            values=q.pivot(index='irradiance_W_m2',columns='tau_I',values='relative_band_loss_pct')
            values=values.loc[coarse_g,coarse_tau]/100
            delta=values-values.loc[1000]
            lookup[name,orientation]=RegularGridInterpolator((coarse_g,coarse_tau),delta.to_numpy())
    directory=out/'coarsened_grid_check'
    directory.mkdir(exist_ok=True)
    coarse,_=aggregate_sites(sites,lookup,coarse_g,directory)
    keys=['site','module','orientation','irradiance_basis']
    comparison=result[keys+['corrected_band_loss_pct']].merge(
        coarse[keys+['corrected_band_loss_pct']],on=keys,suffixes=('_full','_coarse'),validate='one_to_one')
    comparison['coarsening_difference_pp']=comparison.corrected_band_loss_pct_coarse-comparison.corrected_band_loss_pct_full
    comparison.to_csv(out/'interpolation_resolution_check.csv',index=False)
    return comparison


def figures(result,distribution,out):
    primary=result[result.irradiance_basis=='after_uniform_layer']
    with plt.rc_context(STYLE):
        fig,axes=plt.subplots(1,3,figsize=(8.8,3.0))
        for ax,site in zip(axes,energy.SITES):
            q=primary[primary.site==site]
            for _,row in q.iterrows():
                ax.plot([0,1],[row.original_band_loss_pct,row.corrected_band_loss_pct],color=SITE_COLOR[site],
                        marker='o',ms=3,alpha=.65,lw=.8)
            ax.plot([0,1],[q.original_band_loss_pct.median(),q.corrected_band_loss_pct.median()],
                    color='black',marker='D',ms=4,lw=1.5,label='Median')
            ax.set_xticks([0,1],['Original STC LUT','Irradiance correction'])
            ax.tick_params(axis='x',labelsize=6)
            ax.set_title(site)
            ax.set_ylim(bottom=0)
        axes[0].set_ylabel('Additional annual band energy loss (%)')
        axes[0].legend(frameon=False)
        fig.suptitle('PH15 mm · rectangular band, h = 60 mm, f = 0.1 · rain only',fontsize=9)
        fig.tight_layout()
        save(fig,out/'fig_annual_irradiance_effect')
        fig,axes=plt.subplots(1,2,figsize=(8,3))
        for site in energy.SITES:
            d=distribution[distribution.site==site].sort_values(['module','G_upper_W_m2']).copy()
            d['cumulative_fraction']=d.groupby('module').fraction_clean_energy.cumsum()
            q=d.groupby('G_upper_W_m2').cumulative_fraction.median()
            q=q[np.isfinite(q.index)]
            axes[0].plot(np.r_[0,q.index],np.r_[0,100*q],color=SITE_COLOR[site],label=site)
            q=primary[primary.site==site]
            axes[1].plot(np.repeat(list(energy.SITES).index(site),len(q)),q.annual_difference_pp,
                         'o',color=SITE_COLOR[site],ms=4,mfc='none')
        axes[0].set_xlim(0,1200)
        axes[0].set_xlabel('Effective irradiance (W m$^{-2}$)')
        axes[0].set_ylabel('Cumulative clean energy (%)')
        axes[0].legend(frameon=False)
        axes[1].set_xticks(range(3),list(energy.SITES))
        axes[1].set_ylabel('Change in annual band loss (pp)')
        axes[1].axhline(0,color='.5',lw=.5)
        fig.tight_layout()
        save(fig,out/'fig_energy_weighting_and_changes')


def report(result,summary,Gs,out):
    primary=summary[summary.irradiance_basis=='after_uniform_layer']
    lines=['# Annual irradiance-only sensitivity, PH15 mm','','## Scope','',
        'Reference: rectangular band, h=60 mm, f=0.1, no added glass margin, '
        'rain only, all four modules in portrait and landscape, 2016–2025 weather. '
        'Original manuscript/economic outputs are not overwritten.','',
        '**No irradiance cutoff is used in the economic energy calculation.** '
        'The field-test power-ratio analysis uses 400 W/m²; field-test energy '
        'integration uses 20 W/m². Neither filter is applied to the economic model.','',
        'Hourly clean-production weights come from the unchanged energy.site_hours/pv_yield '
        'pipeline, retaining transposition, IAM and Faiman temperature weighting. '
        'Effective irradiance and cell temperature are recovered algebraically from the '
        'two existing clean-yield series with distinct temperature coefficients, and '
        'all four module yields are checked against the recovered inputs. No weather is invented.','',
        'The saved JA Solar landscape, 60 mm trajectory is reused exactly. For another '
        'orientation/module, tau is its power L_s/L_s,ref, as follows from the existing '
        'rain-only mass balance. Event closures, source hour boundaries and midnight '
        'intervals are unchanged. Uniform dust attenuation is reconstructed from the '
        'same saved event totals and deposition rate, and checked against the saved '
        'energy-weighted interval attenuation. All 24 original band-loss values must '
        'reproduce within 0.00002 pp before LTspice annual calculations proceed.','',
        '## Calculation','',
        'For each real hour i, Y_i is clean energy and U_i is uniform-layer transmission. '
        'The original additional band loss is 100 sum(Y_i U_i L_archive(tau_i))/sum(Y_i). '
        'The primary grid irradiance is G_i U_i, because uniform dust attenuates the '
        'irradiance reaching both covered and uncovered cells. An alternative using '
        'weather G_i before uniform attenuation is also exported.','',
        'At fixed electrical temperature 25 °C, paired clean/soiled LTspice calculations '
        'give ΔL(G,tau)=L_25(G,tau)−L_25(1000,tau). We apply '
        'L_corrected=L_archive(tau)+ΔL(G,tau), rather than replacing the archived baseline. '
        'This explicitly separates the irradiance effect from the previously identified '
        '27 °C versus 25 °C nominal-reference discrepancy and STC LUT interpolation differences. '
        'The original tau quantisation (0.005) is retained. Linear interpolation is used '
        'in irradiance and tau; no statistical confidence interval is implied.','',
        f'Electrical grid: G={list(Gs)} W/m²; tau={list(TAUS)}. '
        'All main-grid clean and soiled powers are paired at the same G and 25 °C. '
        'Previously computed identical netlists are reused where possible.','',
        'No hours or energy are silently dropped. Outside the solved G range, the '
        'original STC contribution is retained and a conservative 0–100% band-loss '
        'bound is exported. Values outside physical loss bounds are clipped only '
        'after baseline anchoring; counts and energy effect of any clipping are reported.','',
        '## Annual results','',
        '| Site | Original median (%) | Irradiance-corrected median (%) | Difference of medians (pp) | Paired differences, min–max (pp) |',
        '|---|---:|---:|---:|---:|']
    for r in primary.itertuples(index=False):
        lines.append(f'| {r.site} | {r.original_median_loss_pct:.6f} | {r.corrected_median_loss_pct:.6f} | '
                     f'{r.difference_of_site_medians_pp:+.6f} | {r.min_difference_pp:+.6f} to {r.max_difference_pp:+.6f} |')
    q=result[result.irradiance_basis=='after_uniform_layer']
    lines+=['',f'Minimum supported clean-energy fraction: {q.supported_clean_energy_fraction.min():.9%}. '
        f'Maximum unsupported-domain bound width: '
        f'{(q.corrected_upper_bound_pct-q.corrected_lower_bound_pct).max():.6f} pp. '
        f'Maximum clipping contribution: {q.clipping_effect_pp.abs().max():.9f} pp.','',
        'The difference of site medians is not necessarily the median of paired '
        'configuration differences; both are exported. Annual differences by year '
        'and the clean-energy fractions below 20, 200 and 400 W/m² can be derived '
        'directly from the CSVs.','',
        '## Contribution of low irradiance','',
        '| Site | Clean energy below 20 W/m² (%) | Below 200 W/m² (%) | Below 400 W/m² (%) |',
        '|---|---:|---:|---:|']
    distribution=pd.read_csv(out/'irradiance_energy_distribution.csv')
    for site in energy.SITES:
        fractions=[]
        for threshold in (20,200,400):
            bins=distribution[(distribution.site==site)&(distribution.G_upper_W_m2<=threshold)]
            fractions.append(100*bins.groupby('module').fraction_clean_energy.sum().median())
        lines.append(f'| {site} | {fractions[0]:.6f} | {fractions[1]:.6f} | {fractions[2]:.6f} |')
    lines += ['',
        '## Limits','',
        'This is an irradiance-only sensitivity of the annual energy estimator, '
        'not a validation or a complete replacement economic calculation. Cell '
        'temperature affects the existing clean-energy weights, but the relative '
        'band-loss correction is evaluated at fixed 25 °C because the original '
        'electrical temperature interface remains insufficient. Fixed Rsh and '
        'unvalidated low-irradiance parameters remain assumptions.','',
        'Cleaning is not reoptimised, and manual/robot costs, rain events and dust '
        'trajectories are not modified. These results must not be applied as a '
        'single percentage correction to every economic strategy or height.','',
        '## Reproduction','',
        '    python chapters/annual_irradiance_sensitivity.py --workers 4','',
        'See annual_results_per_configuration.csv, annual_site_summary.csv, '
        'annual_results_by_year.csv, irradiance_energy_distribution.csv, '
        'irradiance_tau_grid.csv, baseline_reproduction.csv and verification.json.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    values='; '.join(f'{r.site}: {r.original_median_loss_pct:.3f}% to {r.corrected_median_loss_pct:.3f}% '
        f'({r.difference_of_site_medians_pp:+.3f} percentage points)' for r in primary.itertuples(index=False))
    si=('We estimated the irradiance-only effect on annual additional band losses using the '
        'existing 2016–2025 hourly weather and simulated rain-driven band-transmittance '
        'trajectories for the PH15 mm scenario (rectangular band, h=60 mm, f=0.1, no manual '
        'cleaning). All positive clean-energy production was retained; no irradiance '
        'threshold was imposed. Each hour was weighted by its clean energy after the '
        'existing uniform-soiling attenuation. The correction L25(G,tau)−L25(1000,tau) '
        'from paired clean/shaded LTspice calculations was added to the original archived '
        'STC loss lookup, preserving the archived baseline and avoiding a nominal-temperature '
        'offset. Effective irradiance after uniform attenuation was used in the primary '
        'calculation. Median annual band losses over four modules and two orientations '
        'changed as follows: '+values+'. '
        'This is a conditional numerical estimate of the irradiance effect; the electrical '
        'loss correction is evaluated at 25 °C and does not assess the unvalidated electrical '
        'temperature dependence. Cleaning policies and economic results were not reoptimised.')
    (out/'SI_text.md').write_text(si+'\n',encoding='utf-8')
    resolution=out/'interpolation_resolution_check.csv'
    if resolution.exists():
        check=pd.read_csv(resolution)
        maximum=check[check.irradiance_basis=='after_uniform_layer'].coarsening_difference_pp.abs().max()
        with (out/'report.md').open('a',encoding='utf-8') as f:
            f.write(f'\n## Interpolation resolution diagnostic\n\nCoarsening the existing grid to '
                f'G=[1,20,100,400,800,1000,1200] and tau=[0,0.2,0.4,0.6,0.8,0.95,1] '
                f'changes the primary annual result by at most {maximum:.6f} pp across the '
                '24 configurations. This is a resolution-sensitivity diagnostic, not an '
                'independent upper bound on interpolation error. Main results use the full grid.\n')
    refinement=out/'refinement_comparison.csv'
    if refinement.exists():
        check=pd.read_csv(refinement)
        maximum=check[check.irradiance_basis=='after_uniform_layer'].refinement_difference_pp.abs().max()
        with (out/'report.md').open('a',encoding='utf-8') as f:
            f.write(f'\nAdding G=150,300,500,700,900 W/m² to the initial 11-level '
                f'grid changes annual band loss by at most {maximum:.6f} pp across '
                'the 24 configurations. Final main results use this refined 16-level grid. '
                'The comparison is retained in refinement_comparison.csv.\n')
    solver_files=list((out/'runs').glob('*.solver.json'))
    alternate_count=sum(json.loads(path.read_text()).get('alternate_used',False) for path in solver_files)
    with (out/'report.md').open('a',encoding='utf-8') as f:
        f.write(f'\n## Numerical solver handling\n\nThe Normal solver is limited to 20 seconds '
            f'per new run. On timeout, its log is retained and the identical netlist is '
            f'retried with the native LTspice Alternate solver (no parameter or tolerance change). '
            f'Alternate was used in {alternate_count} completed run attempts. The paired '
            'SunPower clean/shaded solver comparison is saved in alternate_solver_comparison.csv. '
            'Completed original-grid runs are reused; interrupted calculations are resumed from '
            'per-netlist hash-checked metadata.\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--exe',type=Path,default=Path(os.environ.get('LOCALAPPDATA',''))/'Programs/ADI/LTspice/LTspice.exe')
    parser.add_argument('--library',type=Path,default=Path(os.environ.get('LOCALAPPDATA',''))/'LTspice/lib/cmp/standard.dio')
    args=parser.parse_args()
    scenario=config.CH4_DIR/'total_scenarios/PH15mm'
    out=config.OUTPUT_DIR/'stc_sensitivity/annual_irradiance'
    out.mkdir(parents=True,exist_ok=True)
    initial=out/'initial_resolution_annual_results.csv'
    if (out/'annual_results_per_configuration.csv').exists() and not initial.exists():
        shutil.copyfile(out/'annual_results_per_configuration.csv',initial)
    protected={p:sha(p) for directory in (config.INPUT_DIR,scenario,config.REPO_DIR/'soilspice')
               for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    sites=prepare_sites(scenario,out)
    Gmax=max(w['G'].max() for w in sites.values())
    Gs=tuple(sorted(set(BASE_G+tuple(range(1200,int(math.ceil(Gmax/200)*200)+1,200)))))
    (out/'grid_settings.json').write_text(json.dumps(dict(irradiances=Gs,transmittances=TAUS,
        electrical_temperature_C=25,geometry='rectangle, h=60mm, all modules/orientations',
        rainfall_washing_threshold_mm=15,trapped_fraction=.1,cleaning='none'),indent=2)+'\n')
    if args.prepare_only:
        print('Weather preparation complete; electrical grid:',Gs,flush=True)
        return
    config.LTSPICE_EXE=args.exe.resolve()
    electrical.spice.run_ltspice=run_ltspice_with_fallback  # process-local numerical runner
    runner=electrical.Runner(out,args.library)
    states=grid_states(Gs)
    for state in states:
        previous=out.parent/'runs'/(state.key+'.json')
        if previous.exists():
            cached=json.loads(previous.read_text())
            text,_,_=electrical.netlist(state,args.library)
            if cached.get('status')=='ok' and cached.get('netlist_sha256')==hashlib.sha256(text.encode()).hexdigest():
                cached['reused_from']=str(previous.relative_to(config.REPO_DIR))
                runner.results[state]=cached
    print(f'Reusing {len(runner.results)} previous identical electrical simulations.',flush=True)
    runner.batch(states,args.workers,'Annual G/tau grid')
    verify_alternate_solver(runner)
    lookup=native_grid(runner,Gs)
    result,summary=aggregate_sites(sites,lookup,Gs,out)
    if initial.exists():
        keys=['site','module','orientation','irradiance_basis']
        old=pd.read_csv(initial)[keys+['corrected_band_loss_pct']]
        comparison=result[keys+['corrected_band_loss_pct']].merge(old,on=keys,suffixes=('_refined','_initial'),validate='one_to_one')
        comparison['refinement_difference_pp']=comparison.corrected_band_loss_pct_refined-comparison.corrected_band_loss_pct_initial
        comparison.to_csv(out/'refinement_comparison.csv',index=False)
    interpolation_resolution_check(sites,result,Gs,out)
    figures(result,pd.read_csv(out/'irradiance_energy_distribution.csv'),out)
    report(result,summary,Gs,out)
    modified=[str(p) for p,h in protected.items() if sha(p)!=h]
    (out/'verification.json').write_text(json.dumps(dict(modified_protected_files=modified,
        protected_file_count=len(protected),all_baseline_checks_passed=True,
        irradiance_cutoff_in_economics=None,field_power_cutoff_W_m2=400,field_energy_cutoff_W_m2=20,
        electrical_simulations=len(runner.results),failed_simulations=sum(r['status']!='ok' for r in runner.results.values()),
        protected_files_sha256={str(p.relative_to(config.REPO_DIR)):h for p,h in protected.items()}),indent=2)+'\n')
    assert not modified,modified
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':
    main()
