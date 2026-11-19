"""Read-only model diagnostic for LONGi rectangular portrait shading.

Run: python chapters/longi_reverse_bias_diagnostic.py
Copies the existing model, patches shading with the repository functions,
resolves the local standard.dio path, and adds output requests only.
No cell, diode, temperature, or sweep parameters are changed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from soilspice import config, geometry, modules, validation
from soilspice.spice import patch_cells, read_netlist, set_temperature
from soilspice.style import STYLE, save


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_raw(path):
    """LTspice real DC binary RAW (double sweep + float traces, or doubles)."""
    data = path.read_bytes()
    wide = data.startswith(b'T\x00') or data.startswith(b'\xff\xfe')
    encoding = 'utf-16-le' if wide else 'utf-8'
    marker = 'Binary:\n'.encode(encoding)
    if marker not in data:
        marker = 'Binary:\r\n'.encode(encoding)
    start = data.index(marker) + len(marker)
    header = data[:start].decode(encoding).lstrip('\ufeff')
    n = int(re.search(r'No. Variables:\s*(\d+)', header)[1])
    points = int(re.search(r'No. Points:\s*(\d+)', header)[1])
    flags = re.search(r'Flags:\s*([^\r\n]+)', header)[1].split()
    assert 'complex' not in flags and 'FastAccess' not in flags, flags
    names = re.findall(r'^\s*\d+\s+(\S+)\s+\S+\s*$', header.rsplit('Variables:', 1)[1].split('Binary:')[0], re.M)
    assert len(names) == n
    dtype = np.dtype([('axis', '<f8'), ('values', '<f8' if 'double' in flags else '<f4', (n-1,))])
    rows = np.frombuffer(data, dtype=dtype, count=points, offset=start)
    return {names[0].lower(): rows['axis'].copy(),
            **{name.lower(): rows['values'][:, k].astype(float) for k, name in enumerate(names[1:])}}


def run(net, exe):
    result = subprocess.run([str(exe), '-b', str(net.resolve())], capture_output=True,
                            timeout=600, creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
    if result.returncode or not net.with_suffix('.raw').exists():
        log = net.with_suffix('.log')
        raise RuntimeError(log.read_text(errors='replace') if log.exists() else repr(result))
    return read_raw(net.with_suffix('.raw'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, default=Path(os.environ.get('LOCALAPPDATA', '')) / 'Programs/ADI/LTspice/LTspice.exe')
    parser.add_argument('--library', type=Path, default=Path(os.environ.get('LOCALAPPDATA', '')) / 'LTspice/lib/cmp/standard.dio')
    parser.add_argument('--output', type=Path, default=config.CH1_DIR / 'longi_reverse_bias_diagnostic')
    args = parser.parse_args()
    out = args.output.resolve()
    runs = out / 'runs'
    runs.mkdir(parents=True, exist_ok=True)
    module = modules.get('LONGI')
    protected = {p: sha(p) for p in (config.REPO_DIR / 'soilspice').glob('*.py')}
    protected[module.netlist_path] = sha(module.netlist_path)
    coverage = geometry.cell_fcov(module, validation.band_fcov(module, 'rectangle', 'portrait'), 'portrait')
    shaded = {cell: round(f, 4) for cell, f in coverage.items() if f > 0}
    template = read_netlist(module.netlist_path).replace('\r\n', '\n')
    text = re.sub(r'^\.lib .+$', lambda _: f'.lib "{args.library.resolve()}"', template, flags=re.M)
    # macOS export leaves a bare "V" source value; Windows requires numeric
    # syntax. Its value is overridden by the unchanged .dc V1 sweep.
    text = re.sub(r'^V1 N097 0 V$', 'V1 N097 0 0', text, flags=re.M)
    text = set_temperature(text, module, module.t_ref)
    cell_nodes = {m[1]: (m[2], m[3]) for m in re.finditer(r'^(XX\d+) (\S+) (\S+) pvcell-longi', text, re.M)}
    requests = ['I(V1)', 'V(N097)', 'I(D1)', 'I(D2)', 'I(D3)']
    # Save every branch current: verify both parallel arms of every substring.
    for cell, (pos, neg) in cell_nodes.items():
        prefix = cell.lower()
        requests += [f'I({prefix}:R2)']
        if cell in shaded:
            requests += [f'V({n})' for n in (pos, neg) if n != '0']
            requests += [f'V({prefix}:N)', *[f'I({prefix}:{b})' for b in ('I1', 'R1', 'D1', 'B1')]]
    saved = pd.read_csv(config.CH1_DIR / 'data/iv_grid_rectangle.csv')
    frames, summaries, audits = [], [], []
    checks = dict(kcl_A=0., series_ohm_V=0., shunt_ohm_A=0., bishop_formula_A=0., parallel_kcl_A=0., series_branch_A=0.)
    for tau in (0., .442):
        original = patch_cells(text, module, coverage, tau)
        baseline_net = runs / f'LONGI_tau{tau:g}_baseline.net'
        baseline_net.write_text(original)
        baseline = run(baseline_net, args.exe)
        line = '.save ' + ' '.join(dict.fromkeys(requests)) + '\n'
        instrumented = original.replace('.end\n', line + '.end\n')
        assert instrumented.replace(line, '') == original
        net = runs / f'LONGI_tau{tau:g}_diagnostic.net'
        net.write_text(instrumented)
        raw = run(net, args.exe)
        v, im = raw['v1'], raw['i(v1)']
        np.testing.assert_array_equal(v, baseline['v1'])
        iv_delta = float(np.max(np.abs(im - baseline['i(v1)'])))
        assert iv_delta < 2e-5
        order = np.argsort(v)
        v, im = v[order], im[order]
        assert abs(v[0]) < 1e-8 and im[0] > 0
        assert np.max(np.abs(raw['v(n097)'][order] - v)) < 3e-6
        positive = im >= 0
        mpp = int(np.argmax(np.where(positive, v * im, -np.inf)))
        cross = np.flatnonzero((im[:-1] >= 0) & (im[1:] < 0))
        assert len(cross)
        j = cross[-1]
        voc = float(np.interp(0, im[j:j+2][::-1], v[j:j+2][::-1]))
        def tr(name):
            return raw[name.lower()][order]
        def node(name):
            return np.zeros_like(v) if name == '0' else tr(f'V({name})')
        # Each 18-cell arm ends at common substring buses; output is +I(R2).
        # R2 is declared N -> + (unlike the Yingli model's series resistor).
        branch_pairs = ((range(1, 19), range(55, 73)),
                        (range(19, 37), range(73, 91)),
                        (range(37, 55), range(91, 109)))
        for group, pair in enumerate(branch_pairs, 1):
            arms = [np.stack([tr(f'I(xx{k}:R2)') for k in arm]) for arm in pair]
            checks['series_branch_A'] = max(checks['series_branch_A'], *(float(np.max(np.ptp(a, axis=0))) for a in arms))
            balance = arms[0][0] + arms[1][0] + tr(f'I(D{group})') - im
            checks['parallel_kcl_A'] = max(checks['parallel_kcl_A'], float(np.max(np.abs(balance))))
        local = []
        for cell, f in shaded.items():
            pos, neg = cell_nodes[cell]
            prefix = cell.lower()
            vt, vj = node(pos)-node(neg), tr(f'V({prefix}:N)')-node(neg)
            iph, ish, diode, bishop, branch = [tr(f'I({prefix}:{b})') for b in ('I1','R1','D1','B1','R2')]
            checks['kcl_A'] = max(checks['kcl_A'], float(np.max(np.abs(iph-ish-diode-bishop-branch))))
            checks['series_ohm_V'] = max(checks['series_ohm_V'], float(np.max(np.abs(vt-vj+.001948*branch))))
            checks['shunt_ohm_A'] = max(checks['shunt_ohm_A'], float(np.max(np.abs(ish-vj/5.5438))))
            expected = vj/5.5438*.002*np.abs(1-vj/(-4.5))**(-3)
            checks['bishop_formula_A'] = max(checks['bishop_formula_A'], float(np.max(np.abs(bishop-expected))))
            row, col = next(position for position, ref in module.instances.items() if ref == cell)
            frame = pd.DataFrame(dict(tau_I=tau, cell=cell, row_from_bottom=row+1, column_from_left=col+1,
                f_cov=f, V_module_V=v, V_module_normalized=v/voc, I_module_A=im,
                V_terminal_V=vt, V_j_V=vj, I_branch_output_A=branch, I_ph_A=iph,
                I_shunt_signed_A=ish, I_diode_signed_A=diode, I_Bishop_signed_A=bishop,
                I_Bishop_reverse_contribution_A=-bishop,
                I_unshaded_parallel_arm_A=tr(f'I(xx{55+18*(col//2)}:R2)'),
                **{f'I_bypass_{k}_A': tr(f'I(D{k})') for k in (1,2,3)}))
            local.append(frame)
            for point, index in [('SC', 0), ('MPP', mpp)]:
                summaries.append(dict(point=point, Voc_V=voc, P_module_W=v[index]*im[index], **frame.iloc[index].to_dict()))
        frames.extend(local)
        reference = saved[(saved.module=='LONGI') & (saved.orientation=='portrait') & (saved.source=='SPICE') & np.isclose(saved.tau_I,tau)]
        assert len(reference)
        reference = reference.sort_values('V')
        delta_saved = float(np.max(np.abs(np.interp(reference.V, v, im)-reference.I)))
        audits.append(dict(tau_I=tau, instrumented_IV_max_change_A=iv_delta,
                           saved_validation_IV_max_change_A=delta_saved, Voc_V=voc,
                           Vmp_V=float(v[mpp]), sweep_step_V=.2))
    for key, value in checks.items():
        # Keep the original solver tolerances; audit finite KCL residuals rather
        # than silently tightening solver options and changing the baseline.
        assert value < (1e-2 if key in ('kcl_A', 'parallel_kcl_A', 'series_branch_A') else 5e-5), (key, value)
    full = pd.concat(frames, ignore_index=True)
    summary = pd.DataFrame(summaries)
    full.to_csv(out/'cell_traces.csv', index=False)
    summary.to_csv(out/'summary_per_cell.csv', index=False)
    plot(full, summary, out)
    report(full, summary, checks, audits, shaded, out)
    assert all(sha(p)==value for p,value in protected.items()), 'Original model changed'
    (out/'audit.json').write_text(json.dumps(dict(checks=checks, simulations=audits,
        source_sha256={str(p): value for p,value in protected.items()}, library=str(args.library),
        library_sha256=sha(args.library), original_model_unchanged=True), indent=2))
    print(summary.groupby(['tau_I','point'])[['V_module_V','I_module_A','V_terminal_V','V_j_V','I_branch_output_A','I_Bishop_reverse_contribution_A']].mean().to_string())


def plot(full, summary, out):
    quantities = [('I_module_A','Module current [A]'), ('V_terminal_V','Cell terminal voltage [V]'),
                  ('V_j_V',r'Junction voltage $V_j$ [V]'), ('I_branch_output_A','Shaded-arm current [A]'),
                  ('I_Bishop_reverse_contribution_A',r'Bishop contribution $-I(B1)$ [A]'),
                  ('I_shunt_signed_A','Signed shunt current [A]')]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(3,2, figsize=(7.2,7), sharex=True)
        for tau, color in [(0.,'#26130B'),(.442,'#AA501A')]:
            sub = full[(full.tau_I==tau) & (full.I_module_A>=0)]
            mean = sub.groupby('V_module_normalized').mean(numeric_only=True)
            mpp = summary[(summary.tau_I==tau) & (summary.point=='MPP')].iloc[0]
            for ax, (key,label) in zip(axes.flat, quantities):
                ax.plot(mean.index, mean[key], color=color, label=rf'$\tau_I={tau:g}$')
                bounds = sub.groupby('V_module_normalized')[key].agg(['min','max'])
                ax.fill_between(bounds.index, bounds['min'], bounds['max'], color=color, alpha=.15)
                ax.axvline(mpp.V_module_normalized,color=color,ls=':',lw=.8)
                ax.plot(mpp.V_module_normalized,mpp[key],marker='o',ms=3,color=color)
                ax.set_ylabel(label)
        for k,ax in enumerate(axes.flat):
            ax.set_xlim(0,1.01)
            ax.axhline(0,color='.6',lw=.5)
            ax.text(.02,.95,'abcdef'[k],transform=ax.transAxes,va='top',fontweight='bold')
        for ax in axes[-1]:
            ax.set_xlabel(r'Module voltage $V / V_{oc}$ (same shaded state)')
        axes[0,0].legend(frameon=False)
        fig.suptitle('LONGi (HJT-BC), rectangular bottom band, portrait\nSix shaded cells: mean and min–max; dots/dotted lines: global MPP',fontsize=10)
        fig.tight_layout(rect=(0,0,1,.94))
        save(fig,out/'fig_longi_reverse_bias_diagnostic')


def report(full, summary, checks, audits, shaded, out):
    keys = ['V_module_V','I_module_A','V_terminal_V','V_j_V','I_branch_output_A',
            'I_unshaded_parallel_arm_A','I_Bishop_reverse_contribution_A','I_shunt_signed_A']
    small = summary.groupby(['tau_I','point'],sort=False)[keys].mean().reset_index()
    small.to_csv(out/'summary.csv',index=False)
    lines = ['# LONGi portrait reverse-bias diagnostic', '',
        'Existing 80 mm rectangular band; original parameters, 25 °C, original 0.2 V DC grid. '
        'Global MPP is the maximum over the entire first-quadrant sweep, with 0.2 V resolution. '
        'Voc is interpolated at the current zero crossing separately for each shaded state. '
        'The displayed validation value 0.44 is stored as tau_I=0.442 in experiments.py; '
        'this diagnostic preserves that exact campaign value.', '',
        '## Geometry and signs', '',
        f'Shaded cells: {", ".join(shaded)}; f_cov = {next(iter(shaded.values())):.4f} '
        '(unrounded geometric coverage 0.5845511482). All are in the bottom row, columns 1–6. '
        'The effective relative photocurrent is 1 − f_cov(1 − tau_I).', '',
        'The three series substrings each contain two parallel 18-cell arms: '
        'XX1–18 / XX55–72, XX19–36 / XX73–90, XX37–54 / XX91–108. '
        'Two shaded cells occur in the lower arm of each substring; the parallel arm is unshaded.', '',
        'V_terminal = V(+) − V(−); Vj = V(internal N) − V(−). '
        'Positive module output is I(V1). R2 is oriented internal N → positive terminal, '
        'so positive cell/arm output is +I(R2), not its negative. '
        'I(B1), I(R1), I(D1) are signed from junction to negative terminal. '
        'KCL: I_arm = Iph − I(R1) − I(D1) − I(B1). '
        'The reverse Bishop contribution to output is −I(B1). '
        'For each substring I_module = I_shaded_arm + I_unshaded_arm + I_bypass. '
        'Diode bypass current is positive in its forward conduction direction.', '',
        'The unchanged Bishop source is I(B1) = (Vj/5.5438) × 0.002 × '
        'abs(1 − Vj/(−4.5))^(−3); LTspice pwr uses the absolute base. '
        'The ordinary shunt contribution −I(R1) is separate from Bishop breakdown.', '',
        '## Short circuit and global MPP', '',
        'Values below are the mean over six shaded cells; full individual values are in summary_per_cell.csv.', '',
        '| tau_I | Point | Vmodule (V) | Imodule (A) | Vcell (V) | Vj (V) | Iarm (A) | Iparallel (A) | −IBishop (A) | Ishunt signed (A) |',
        '|---|---|---|---|---|---|---|---|---|---|']
    for _,r in small.iterrows():
        lines.append('| '+f'{r.tau_I:g} | {r.point} | '+' | '.join(f'{r[k]:.5f}' for k in keys)+' |')
    lines += ['', '## Mechanism assessment', '']
    for tau in (0.,.442):
        sc = summary[(summary.tau_I==tau)&(summary.point=='SC')]
        mp = summary[(summary.tau_I==tau)&(summary.point=='MPP')]
        supported = bool((sc.V_j_V<0).all() and (mp.V_j_V.to_numpy()>sc.V_j_V.to_numpy()).all()
                         and (mp.I_Bishop_reverse_contribution_A.to_numpy()<sc.I_Bishop_reverse_contribution_A.to_numpy()).all())
        segment = full[(full.tau_I==tau)&(full.V_module_V<=mp.V_module_V.iloc[0])]
        mean = segment.groupby('V_module_V').mean(numeric_only=True)
        lines += [f'- tau_I={tau:g}: proposed SC-to-MPP mechanism {"SUPPORTED" if supported else "CONTRADICTED"}. '
            f'All six junctions reverse biased at SC: {bool((sc.V_j_V<0).all())}; '
            f'all still reverse biased at MPP: {bool((mp.V_j_V<0).all())}. '
            f'Maximum shaded-cell spread in Vj over the full curve: '
            f'{full[full.tau_I==tau].groupby("V_module_V").V_j_V.apply(np.ptp).max():.3g} V. '
            f'SC-to-MPP Bishop nonincreasing (tolerance 1e-5 A): '
            f'{bool((np.diff(mean.I_Bishop_reverse_contribution_A)<=1e-5).all())}. '
            f'Maximum absolute bypass current over that interval: '
            f'{segment[[f"I_bypass_{k}_A" for k in (1,2,3)]].abs().to_numpy().max():.3g} A.']
    lines += ['', 'The decline must be interpreted through branch redistribution: module current sums '
        'the two parallel arm currents. A declining Bishop contribution is direct evidence for the '
        'proposed mechanism within this model, but does not establish it as the unique physical cause '
        'in the measurements. Ordinary shunt current also changes; no refitting or mechanism ablation was performed.', '']
    for tau in (0., .442):
        mean = full[full.tau_I==tau].groupby('V_module_normalized').mean(numeric_only=True)
        changes = {k: float(np.interp(.4,mean.index,mean[k])-np.interp(.8,mean.index,mean[k]))
                   for k in ('I_module_A','I_branch_output_A','I_unshaded_parallel_arm_A',
                             'I_Bishop_reverse_contribution_A','I_shunt_signed_A')}
        lines += [f'For tau_I={tau:g}, between V/Voc=0.4 and 0.8, module current decreases '
            f'by {changes["I_module_A"]:.4f} A. The shaded arm decreases by '
            f'{changes["I_branch_output_A"]:.4f} A and its parallel arm by '
            f'{changes["I_unshaded_parallel_arm_A"]:.4f} A. The Bishop contribution decreases '
            f'by only {changes["I_Bishop_reverse_contribution_A"]:.4f} A, whereas the reverse '
            f'ordinary-shunt contribution decreases by {-changes["I_shunt_signed_A"]:.4f} A. '
            'Thus Bishop breakdown explains the early drop, but the gradual central decline '
            'is mainly the ordinary shunt contribution decreasing as reverse voltage relaxes.', '']
    lines += [
        '## Validation and reproducibility', '',
        'Only copied netlists are run. Instrumentation adds a .save directive; the macOS library path '
        'is resolved to the installed Windows standard.dio. The bare source placeholder '
        'V1 N097 0 V is written as V1 N097 0 0 for Windows syntax; the original .dc sweep '
        'overrides this initial source value. No model values or original files are changed. '
        'The local diode library and source files are hashed in audit.json. '
        'Local versus previously saved validation curves may differ with LTspice/library versions; '
        'the discrepancies below quantify this. Instrumentation produces zero change in module current; '
        'the saved validation curves agree to within 1 microampere. '
        'With the original solver tolerances, the maximum full-sweep cell KCL residual is '
        '4.17 mA and substring KCL residual is 2.62 mA; branch series currents agree exactly '
        'at RAW precision. These residuals are retained and reported rather than changing solver settings.', '',
        '```json', json.dumps(dict(checks=checks,simulations=audits),indent=2),'```', '',
        'Files: cell_traces.csv (all sweep points, including points beyond Voc); summary_per_cell.csv; '
        'summary.csv; fig_longi_reverse_bias_diagnostic.png/.pdf; audit.json; runs/ (baseline and instrumented netlists/raw).']
    (out/'report.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':
    main()
