"""Diagnostic reverse-only shunt sweep for Yingli, without changing the production model.

    MPLBACKEND=Agg python chapters/yingli_reverse_shunt_sensitivity.py

Original constant-Rsh netlists are copied. Only their shunt element is replaced
by a continuous piecewise-linear current source; forward resistance stays 5.99
ohm. All other circuit lines and electrical parameters are preserved verbatim.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PyLTSpice import RawRead

from soilspice import config, geometry, modules, validation
from soilspice.spice import iv_parameters, run_ltspice
from soilspice.style import STYLE, save

NEGLIGIBLE_OHM = 1.e6
LIMIT_CHECK_OHM = 1.e7
RESISTANCES = (5.99, 10., 20., 50., 100., 200., 500., NEGLIGIBLE_OHM)
ORIGINAL_SHUNT = "Rp N001 - {rp}"


def label(resistance):
    return "negligible" if resistance == NEGLIGIBLE_OHM else f"{resistance:g}"


def shunt_current(voltage, resistance):
    """Diagnostic law: continuous current at zero, with a conductance kink."""
    voltage = np.asarray(voltage, dtype=float)
    return np.maximum(voltage, 0.) / 5.99 + np.minimum(voltage, 0.) / resistance


def patch_shunt(text, resistance):
    if resistance not in (*RESISTANCES, LIMIT_CHECK_OHM):
        raise ValueError("Resistance outside the prescribed diagnostic sweep")
    assert text.count(ORIGINAL_SHUNT) == 1
    reverse = f"V(N001,-)/{resistance:g}"
    replacement = f"Bsh N001 - I=if(V(N001,-)>=0,V(N001,-)/{{rp}},{reverse})"
    result = text.replace(ORIGINAL_SHUNT, replacement)
    assert result.replace(replacement, ORIGINAL_SHUNT) == text
    return result


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trace(raw, name):
    names = {n.lower(): n for n in raw.get_trace_names()}
    return np.asarray(raw.get_trace(names[name.lower()]).get_wave(), dtype=float)


def instrument(text, cells, initial_guess=""):
    requests = ["I(V1)", "I(D1)", "I(D2)", "I(D3)"]
    for cell, negative in cells:
        prefix = cell[1:].lower()
        if negative != "0":
            requests.append(f"V({negative})")
        requests += [f"V({prefix}:N001)", *(f"I({prefix}:{branch})" for branch in ("Bsh", "D1", "I1", "Rs"))]
    # Tighten numerical convergence at the conductance kink, without changing any circuit parameter.
    directive = (".save " + " ".join(dict.fromkeys(requests)) +
                 "\n.options reltol=1e-5 abstol=1e-10 vntol=1e-8\n" + initial_guess)
    assert text.count("\n.end\n") == 1
    result = text.replace("\n.end\n", "\n" + directive + ".end\n")
    assert result.replace(directive, "") == text
    return result


def evaluate(measured, simulations, orientation):
    """Identical voltage samples across all eight resistance cases."""
    low = max(measured.v[0], *(curve[0][0] for curve in simulations.values()))
    high = min(measured.v[-1], *(curve[0][-1] for curve in simulations.values()))
    voltage = np.arange(low, high, .1)
    actual = np.interp(voltage, measured.v, measured.i)
    valid = actual >= .3
    plateau = valid & (voltage >= 5.) & (voltage <= 25.)
    limiting = np.interp(voltage, *simulations[NEGLIGIBLE_OHM])
    rows = []
    for resistance, (v, i) in simulations.items():
        predicted = np.interp(voltage, v, i)
        error = predicted - actual
        pmax = iv_parameters(v, i)["Pmax_W"]
        rows.append(dict(orientation=orientation, tau_I=measured.tau_i,
                         R_reverse_ohm=resistance, R_label=label(resistance),
                         n_voltage_points=int(valid.sum()),
                         nRMSE_IV_pct=100 * np.sqrt(np.mean(error[valid] ** 2)) / measured.params["Isc_A"],
                         nMBE_IV_pct=100 * np.mean(error[valid]) / measured.params["Isc_A"],
                         nRMSE_PV_pct=100 * np.sqrt(np.mean((voltage[valid] * error[valid]) ** 2)) /
                                       measured.params["Pmax_W"],
                         plateau_5_25V_RMSE_A=np.sqrt(np.mean(error[plateau] ** 2)),
                         plateau_5_25V_bias_A=np.mean(error[plateau]),
                         Pmax_measured_W=measured.params["Pmax_W"], Pmax_simulated_W=pmax,
                         Pmax_error_W=pmax - measured.params["Pmax_W"],
                         Pmax_error_pct=100 * (pmax / measured.params["Pmax_W"] - 1),
                         max_current_distance_to_limit_5_25V_A=np.max(np.abs(predicted[plateau] - limiting[plateau])),
                         mean_current_distance_to_limit_5_25V_A=np.mean(np.abs(predicted[plateau] - limiting[plateau]))))
    return rows


def make_figures(states, metrics, output):
    colors = [".35", *[plt.get_cmap("viridis")(v) for v in np.linspace(.08, .9, 6)], "#D55E00"]
    styles = ["--", "-", "-", "-", "-", "-", "-", "-."]
    handles = [Line2D([], [], color="black", marker="o", ls="none", ms=3, label="Measured")]
    handles += [Line2D([], [], color=c, ls=ls, label=("Negligible (1 MΩ)" if r == NEGLIGIBLE_OHM
                else f"{r:g} Ω")) for r, c, ls in zip(RESISTANCES, colors, styles)]
    with plt.rc_context(STYLE):
        for orientation in ("portrait", "landscape"):
            taus = sorted((tau for o, tau in states if o == orientation), reverse=True)
            fig, axes = plt.subplots(2, len(taus), figsize=(10.5, 4.4), squeeze=False)
            for col, tau in enumerate(taus):
                measured, simulations = states[orientation, tau]
                for row in (0, 1):
                    ax = axes[row, col]
                    for resistance, color, ls in zip(RESISTANCES, colors, styles):
                        ax.plot(*simulations[resistance], color=color, ls=ls, lw=1.)
                    ax.plot(measured.v, measured.i, color="black", marker="o", ls="none", ms=2.3,
                            markevery=max(1, len(measured.v) // 40), mew=0, zorder=5)
                    if row == 0:
                        ax.set_xlim(0, 39)
                        ax.set_ylim(0, 10)
                        ax.set_title(rf"$\tau_I={tau:g}$")
                        ax.axvspan(5, 25, color=".6", alpha=.06, zorder=0)
                    else:
                        ax.set_xlim(5, 25)
                        samples = np.linspace(5, 25, 101)
                        values = [np.interp(samples, measured.v, measured.i)]
                        values += [np.interp(samples, *curve) for curve in simulations.values()]
                        low, high = min(np.min(a) for a in values), max(np.max(a) for a in values)
                        padding = max(.04, .08 * (high - low))
                        ax.set_ylim(low - padding, high + padding)
                        ax.set_xlabel("Module voltage [V]")
                    if col == 0:
                        ax.set_ylabel("Current [A]" if row == 0 else "Current [A], 5–25 V detail")
            fig.suptitle(f"Yingli, rectangular 80 mm band, {orientation}: reverse-only shunt sensitivity", fontsize=9)
            fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=6.5)
            fig.tight_layout(rect=(0, .13, 1, .94), h_pad=1.1)
            save(fig, output / f"fig_reverse_shunt_{orientation}")

        fig, axes = plt.subplots(2, 3, figsize=(9.0, 5.0))
        for row, orientation in enumerate(("portrait", "landscape")):
            for tau in sorted(metrics.tau_I.unique(), reverse=True):
                q = metrics[(metrics.orientation == orientation) & (metrics.tau_I == tau)]
                if q.empty:
                    continue
                q = q.set_index("R_reverse_ohm").loc[list(RESISTANCES)]
                color = plt.get_cmap("plasma")({1.: .1, .706: .3, .536: .5, .442: .7, 0.: .9}[tau])
                for col, metric in enumerate(("nRMSE_IV_pct", "plateau_5_25V_RMSE_A", "Pmax_error_pct")):
                    axes[row, col].plot(range(8), q[metric], marker="o", ms=2.5, lw=1,
                                        color=color, label=rf"$\tau_I={tau:g}$")
            for col, title in enumerate(("I–V nRMSE [% of measured Isc]", "5–25 V current RMSE [A]",
                                         "Signed maximum-power error [%]")):
                ax = axes[row, col]
                ax.set_xticks(range(8))
                ax.set_xticklabels(["5.99", "10", "20", "50", "100", "200", "500", "negligible"], rotation=45)
                ax.set_xlabel(r"Reverse $R_{\rm sh}$ [Ω], discrete cases")
                ax.set_title(f"{orientation.capitalize()}: {title}", fontsize=7)
                if col == 2:
                    ax.axhline(0, color=".5", ls=":", lw=.5)
        fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", ncol=5)
        fig.tight_layout(rect=(0, .055, 1, 1))
        save(fig, output / "fig_reverse_shunt_error_sensitivity")


def write_report(metrics, summary, verification, output):
    lines = ["# Yingli reverse-only shunt sensitivity", "",
             "Diagnostic study, not a calibrated physical model. Reproduce with "
             "`MPLBACKEND=Agg python chapters/yingli_reverse_shunt_sensitivity.py`; "
             "add `--reuse-simulations` to reuse matching saved netlists/raw files.", "",
             "## Model and scope", "",
             "Only the shunt element in copies of the original constant-5.99-ohm Yingli netlists changes. "
             "Its current is Vj/5.99 for Vj >= 0 and Vj/Rreverse for Vj < 0. Current is continuous at zero; "
             "its derivative may have a kink. The law applies to every cell, so forward-biased cells retain "
             "the original conductance. Photocurrent, geometry, series resistance, diode parameters, "
             "bypass diodes and temperature handling remain unchanged. No Bishop term is added. "
             "The production module definition and all pre-existing results are preserved.", "",
             "Rectangular 80 mm bands are compared at every available portrait/landscape campaign "
             "transmittance, including the clean control. The 0.442 landscape case was not measured.", "",
             "Exact zero reverse conductance makes reverse-biased identical-cell voltage sharing "
             "numerically ill-conditioned. Negligible conduction is represented by 1 MΩ, explicitly "
             "checked against 10 MΩ. Their maximum module-current difference in the 5–25 V diagnostic "
             f"window is {verification['max_limit_check_current_change_A'] * 1e6:.3f} µA. "
             "The 1 MΩ case retains about 5 µA or less of reverse shunt current in sampled cell states; "
             "it is an approximation to the open-shunt limit, not a literal infinite-resistance simulation.", "",
             "Numerical convergence uses reltol=1e-5, abstol=1e-10 and vntol=1e-8, with initial node "
             "guesses from the original short-circuit solution. These affect numerical solution only, "
             "not the circuit equations. The 10 MΩ convergence check uses the original solver defaults.", "",
             "## Metric definitions", "",
             "Both orientations use the existing clean-flash temperature corrections. All eight models "
             "of a configuration are evaluated on one common 0.1 V grid over their shared measured "
             "voltage support. Measured currents below 0.3 A are excluded. I–V nRMSE is normalized "
             "by measured Isc. The 5–25 V RMSE is in amperes and includes bypass transitions in "
             "landscape, so it is a fixed current window rather than an isolated fitted plateau. "
             "Each configuration has equal weight in the shaded-curve mean; clean curves are excluded "
             "from that mean. Minor differences from earlier pairwise metrics can result from using "
             "the shared support of all eight cases.", "",
             "Maximum-power error means 100*(Pmax_sim/Pmax_measured−1), signed; positive means "
             "overprediction. Pmax is sampled on the original 0.05 V LTspice sweep without refitting.", ""]
    for orientation in ("portrait", "landscape"):
        q = metrics[metrics.orientation == orientation]
        taus = sorted(q.tau_I.unique(), reverse=True)
        lines += [f"## {orientation.capitalize()}: I–V nRMSE [%]", "",
                  "| Reverse Rsh [Ω] | " + " | ".join(f"tau={t:g}" for t in taus) + " |",
                  "|---|" + "---:|" * len(taus)]
        for r in RESISTANCES:
            a = q[q.R_reverse_ohm == r].set_index("tau_I")
            lines.append("| " + label(r) + " | " + " | ".join(f"{a.loc[t, 'nRMSE_IV_pct']:.6f}" for t in taus) + " |")
        lines += ["", f"## {orientation.capitalize()}: current RMSE in 5–25 V [A]", "",
                  "| Reverse Rsh [Ω] | " + " | ".join(f"tau={t:g}" for t in taus) + " |",
                  "|---|" + "---:|" * len(taus)]
        for r in RESISTANCES:
            a = q[q.R_reverse_ohm == r].set_index("tau_I")
            lines.append("| " + label(r) + " | " + " | ".join(f"{a.loc[t, 'plateau_5_25V_RMSE_A']:.6f}" for t in taus) + " |")
        lines.append("")
    lines += ["## Maximum power: invariant across all eight resistances", "",
              "The maximum-power point lies on the unchanged forward-bias branch in these cases. "
              "All eight sampled Pmax values are identical for each configuration; a lower plateau "
              "error therefore does not imply improved maximum-power agreement.", "",
              "| Orientation | tau | Measured Pmax [W] | Simulated Pmax [W] | Error [W] | Error [%] |",
              "|---|---:|---:|---:|---:|---:|"]
    for (orientation, tau), q in metrics.groupby(["orientation", "tau_I"], sort=False):
        assert q.Pmax_simulated_W.max() - q.Pmax_simulated_W.min() < 1e-6
        a = q.iloc[0]
        lines.append(f"| {orientation} | {tau:g} | {a.Pmax_measured_W:.6f} | {a.Pmax_simulated_W:.6f} | "
                     f"{a.Pmax_error_W:+.6f} | {a.Pmax_error_pct:+.6f} |")
    lines += ["", "## Findings", "",
              "Increasing reverse Rsh removes excess portrait plateau current, but agreement is "
              "not monotonic across all resistances. The translucent cases eventually undershoot "
              "the measured current when reverse leakage is suppressed too far. The opaque case "
              "continues to improve toward negligible conduction, while retaining a residual mismatch.", "",
              "Lowest errors among the tested values are descriptive sensitivity results, not "
              "calibration or recommended physical parameters:", "",
              "| Orientation | tau | Lowest I–V nRMSE at [Ω] | Lowest 5–25 V RMSE at [Ω] |",
              "|---|---:|---:|---:|"]
    for (orientation, tau), q in metrics[metrics.tau_I != 1].groupby(["orientation", "tau_I"], sort=False):
        a, b = q.loc[q.nRMSE_IV_pct.idxmin()], q.loc[q.plateau_5_25V_RMSE_A.idxmin()]
        lines.append(f"| {orientation} | {tau:g} | {a.R_label} | {b.R_label} |")
    lines += ["", "Convergence of simulated current toward negligible leakage is quantified separately "
              "from improvement against measurements:", "",
              "| Reverse Rsh [Ω] | Maximum distance to negligible case, portrait [mA] | Landscape [mA] |",
              "|---:|---:|---:|"]
    for r in (100., 200., 500.):
        maxima = [metrics[(metrics.orientation == o) & (metrics.tau_I != 1) &
                         (metrics.R_reverse_ohm == r)].max_current_distance_to_limit_5_25V_A.max() * 1000
                  for o in ("portrait", "landscape")]
        lines.append(f"| {r:g} | {maxima[0]:.6f} | {maxima[1]:.6f} |")
    lines += ["", "Portrait mean nRMSE is lowest at 50 Ω among the tested cases. Landscape mean "
              "nRMSE changes only slightly; none of its shaded cases has a higher whole-curve nRMSE "
              "than the original baseline. The clean curve is unchanged. These observations support "
              "reverse leakage as an important source of the portrait discrepancy, without establishing "
              "a unique or physically calibrated reverse shunt resistance.", "",
              "## Artifacts and verification", "",
              "- `fig_reverse_shunt_portrait.pdf/png`: full curves and 5–25 V detail for all portrait transmittances.",
              "- `fig_reverse_shunt_landscape.pdf/png`: corresponding landscape checks.",
              "- `fig_reverse_shunt_error_sensitivity.pdf/png`: nRMSE, current-window RMSE and signed Pmax error.",
              "- `curve_metrics.csv`: full-precision errors for all 72 resistance/configuration combinations.",
              "- `iv_curves.csv`: measured and simulated curve samples.",
              "- `mean_shaded_curve_metrics.csv`: arithmetic means across shaded transmittances.",
              "- `sampled_cell_states.csv`: direct internal branch-current checks at module 0, 5, 15, 25 and 35 V.",
              "- `convergence_to_limit.csv` and `negligible_conduction_verification.csv`: endpoint checks.",
              "- `verification.json`: unchanged-input hashes and numerical checks.", "",
              f"Maximum difference from the original 5.99 Ω module-current arrays: {verification['max_baseline_current_difference_A']:.9g} A. "
              f"Maximum KCL residual: {verification['max_KCL_residual_A'] * 1000:.6f} mA. "
              f"Maximum piecewise-shunt current-law residual: {verification['max_piecewise_shunt_law_residual_A'] * 1e6:.6f} µA. "
              f"Preserved pre-existing scientific/library files: {verification['preserved_existing_files']}.", "",
              "The inherited `V1 N001 0 V` line generates the same missing-unswept-value warning as "
              "the original model; the explicit DC sweep is retained and the original baseline is reproduced. "
              "`runs/exploratory/` contains preliminary convergence probes excluded from all reported results.", ""]
    (output / "README.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=config.CH1_DIR / "yingli_shunt_comparison/runs/constant")
    parser.add_argument("--output", type=Path, default=config.CH1_DIR / "yingli_reverse_shunt_sensitivity")
    parser.add_argument("--reuse-simulations", action="store_true", help="reuse raw outputs only if the generated netlist is identical")
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    (output / "runs").mkdir(parents=True, exist_ok=True)
    protected = {p: digest(p) for p in config.OUTPUT_DIR.rglob("*") if p.is_file() and output not in p.parents}
    for directory in (config.REPO_DIR / "soilspice", config.NETLIST_DIR):
        protected.update({p: digest(p) for p in directory.rglob("*") if p.is_file()})
    module = replace(modules.get("YINGLI"), rp_law=None)
    states, curves, metric_rows, branch_rows, source_manifest, cache = {}, [], [], [], [], {}
    baseline_max_change, kcl_max, shunt_max, forward_max = 0., 0., 0., 0.
    limit_checks = []
    for orientation in ("portrait", "landscape"):
        coverage = geometry.cell_fcov(module, validation.band_fcov(module, "rectangle", orientation), orientation)
        selected = [cell for cell, f in coverage.items() if abs(f - max(coverage.values())) < 1e-9]
        for measured in validation.load_measurements(module, "rectangle", orientation):
            tau = measured.tau_i
            stem = "clean_YINGLI" if tau == 1 else f"YINGLI_rectangle_{orientation}_tau{tau:g}_constant"
            net_path, raw_path = source / f"{stem}.net", source / f"{stem}.raw"
            text = net_path.read_text()
            baseline = RawRead(str(raw_path))
            initial_nodes = {name: trace(baseline, name)[0] for name in baseline.get_trace_names()
                             if name.startswith("V(")}
            initial_guess = ".nodeset " + " ".join(f"{name}={value:.12g}" for name, value in initial_nodes.items()) + "\n"
            assert not re.search(r"^XX.*\brp=", text, re.MULTILINE)
            cells = [(cell, re.search(rf"^{cell}\s+\S+\s+(\S+)\s", text, re.MULTILINE)[1]) for cell in selected]
            simulations = {}
            source_manifest.append(dict(orientation=orientation, tau_I=tau, net_path=str(net_path),
                                        net_sha256=digest(net_path), raw_sha256=digest(raw_path),
                                        measurement_path=str(measured.path), measurement_sha256=digest(measured.path)))
            for resistance in RESISTANCES:
                key = (stem, resistance, tuple(cells))
                generated = instrument(patch_shunt(text, resistance), cells, initial_guess)
                path = output / "runs" / f"{stem}_Rrev{label(resistance)}_{orientation}.net"
                if key not in cache:
                    if not (args.reuse_simulations and path.exists() and path.read_text() == generated
                            and path.with_suffix(".raw").exists() and path.with_suffix(".raw").stat().st_size > 8192):
                        path.write_text(generated)
                        run_ltspice(path)
                    cache[key] = RawRead(str(path.with_suffix(".raw")))
                raw = cache[key]
                voltage, current = trace(raw, "v1"), trace(raw, "I(V1)")
                if resistance == 5.99:
                    np.testing.assert_array_equal(voltage, trace(baseline, "v1"))
                    change = float(np.max(np.abs(current - trace(baseline, "I(V1)"))))
                    assert change < 2e-5, f"Behavioural 5.99-ohm baseline differs by {change} A"
                    baseline_max_change = max(baseline_max_change, change)
                keep = current >= 0
                simulations[resistance] = (voltage[keep], current[keep])
                curves.append(pd.DataFrame(dict(orientation=orientation, tau_I=tau, source="SPICE",
                                                R_reverse_ohm=resistance, V=voltage[keep], I=current[keep])))
                for cell, negative in cells:
                    prefix = cell[1:].lower()
                    vminus = np.zeros_like(voltage) if negative == "0" else trace(raw, f"V({negative})")
                    junction = trace(raw, f"V({prefix}:N001)") - vminus
                    ish, iph, diode, series = [trace(raw, f"I({prefix}:{name})") for name in ("Bsh", "I1", "D1", "Rs")]
                    kcl_max = max(kcl_max, float(np.max(np.abs(iph - diode - ish + series))))
                    shunt_max = max(shunt_max, float(np.max(np.abs(ish - shunt_current(junction, resistance)))))
                    forward = junction > 1e-5
                    if forward.any():
                        forward_max = max(forward_max, float(np.max(np.abs(ish[forward] - junction[forward] / 5.99))))
                    for point in (0., 5., 15., 25., 35.):
                        branch_rows.append(dict(orientation=orientation, tau_I=tau, cell=cell,
                                                R_reverse_ohm=resistance, V_module_V=point,
                                                V_junction_V=np.interp(point, voltage, junction),
                                                I_ph_A=np.interp(point, voltage, iph),
                                                I_diode_A=np.interp(point, voltage, diode),
                                                I_shunt_A=np.interp(point, voltage, ish)))
                print(f"{orientation} tau={tau:g} reverse Rsh={label(resistance)} Ω done", flush=True)
            check_path = output / "runs" / f"{stem}_limit_check_10Mohm_{orientation}.net"
            # Verify negligible leakage in the diagnostic 5-25 V window. The more
            # ill-conditioned 10-Mohm check uses the original solver defaults.
            check_text = instrument(patch_shunt(text, LIMIT_CHECK_OHM), cells).replace(
                ".options reltol=1e-5 abstol=1e-10 vntol=1e-8\n", "")
            if not (args.reuse_simulations and check_path.exists() and check_path.read_text() == check_text
                    and check_path.with_suffix(".raw").exists() and check_path.with_suffix(".raw").stat().st_size > 8192):
                check_path.write_text(check_text)
                run_ltspice(check_path)
            check_raw = RawRead(str(check_path.with_suffix(".raw")))
            cv, ci = trace(check_raw, "v1"), trace(check_raw, "I(V1)")
            rv, ri = simulations[NEGLIGIBLE_OHM]
            common = (cv >= 5. - 1e-8) & (cv <= 25. + 1e-8)
            limit_delta = float(np.max(np.abs(ci[common] - np.interp(cv[common], rv, ri))))
            assert limit_delta < 2e-5, f"Negligible-shunt endpoint not converged: {limit_delta} A"
            limit_checks.append(dict(orientation=orientation, tau_I=tau, R_primary_ohm=NEGLIGIBLE_OHM,
                                     R_check_ohm=LIMIT_CHECK_OHM, voltage_window="5-25 V",
                                     max_current_change_A=limit_delta))
            if tau == 1.:
                for curve in simulations.values():
                    np.testing.assert_array_equal(curve[0], simulations[5.99][0])
                    np.testing.assert_allclose(curve[1], simulations[5.99][1], atol=2e-5, rtol=0)
            states[orientation, tau] = (measured, simulations)
            curves.append(pd.DataFrame(dict(orientation=orientation, tau_I=tau, source="measured",
                                            R_reverse_ohm=np.nan, V=measured.v, I=measured.i)))
            metric_rows.extend(evaluate(measured, simulations, orientation))
    # KCL tolerance matches reltol=1e-5 for the approximately 10-A circuit.
    assert kcl_max < 1e-4 and max(shunt_max, forward_max) < 3e-5, (kcl_max, shunt_max, forward_max)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output / "curve_metrics.csv", index=False)
    pd.concat(curves, ignore_index=True).to_csv(output / "iv_curves.csv", index=False)
    pd.DataFrame(branch_rows).to_csv(output / "sampled_cell_states.csv", index=False)
    pd.DataFrame(limit_checks).to_csv(output / "negligible_conduction_verification.csv", index=False)
    summary = metrics[metrics.tau_I != 1.].groupby(["orientation", "R_reverse_ohm", "R_label"], sort=False).mean(numeric_only=True).reset_index()
    summary.to_csv(output / "mean_shaded_curve_metrics.csv", index=False)
    metrics[["orientation", "tau_I", "R_reverse_ohm", "max_current_distance_to_limit_5_25V_A",
             "mean_current_distance_to_limit_5_25V_A"]].to_csv(output / "convergence_to_limit.csv", index=False)
    make_figures(states, metrics, output)
    for path, sha in protected.items():
        assert digest(path) == sha, f"Pre-existing file modified: {path}"
    verification = dict(diagnostic_only=True, original_models_unchanged=True, other_modules_unchanged=True,
                        R_forward_ohm=5.99, R_reverse_ohm=[label(r) for r in RESISTANCES],
                        shunt_law="Ish(Vj)=Vj/5.99 if Vj>=0; Vj/Rrev if Vj<0",
                        negligible_conduction_ohm=NEGLIGIBLE_OHM, limit_convergence_check_ohm=LIMIT_CHECK_OHM,
                        max_limit_check_current_change_A=max(r['max_current_change_A'] for r in limit_checks),
                        exact_zero_reverse_conductance_not_used="Ill-conditioned voltage sharing between identical reverse-biased cells; use verified negligible-current limit",
                        continuous_at_zero=True, derivative_continuity_required=False, bishop_terms_added=False,
                        all_other_electrical_lines_preserved_verbatim=True,
                        numerical_solver_options="reltol=1e-5 abstol=1e-10 vntol=1e-8; tighten convergence at the conductance kink",
                        solver_initial_guess="nodeset from original constant-shunt short-circuit node voltages; no final equation or electrical parameter changed",
                        max_baseline_current_difference_A=baseline_max_change,
                        max_KCL_residual_A=kcl_max, max_piecewise_shunt_law_residual_A=shunt_max,
                        max_forward_conductance_residual_A=forward_max, clean_curves_unchanged=True,
                        measurement_count=len(states), shaded_measurement_count=sum(t != 1 for _, t in states),
                        common_grid_step_V=.1, measured_current_cutoff_A=.3, plateau_window_V=[5, 25],
                        Pmax_error_definition="100*(Pmax_simulated/Pmax_measured-1), signed",
                        preserved_existing_files=len(protected), sources=source_manifest,
                        generator_sha256=digest(Path(__file__)))
    (output / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    write_report(metrics, summary, verification, output)
    print(summary[["orientation", "R_label", "nRMSE_IV_pct", "plateau_5_25V_RMSE_A", "Pmax_error_pct"]].to_string(index=False))


if __name__ == "__main__":
    main()
