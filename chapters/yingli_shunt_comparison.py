"""Compare constant and PVsyst Rsh against every Yingli campaign I-V curve.

    MPLBACKEND=Agg python chapters/yingli_shunt_comparison.py

Uses the existing campaign geometry, temperature corrections and plotting style.
All new simulations and plots are isolated from the published validation exports.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from soilspice import config, experiments, geometry, modules, plots, validation
from soilspice.spice import Simulator, shaded_rp
from soilspice.style import STYLE, save


def paired_errors(measured, simulations):
    """Same voltage samples and measured normalization for both models."""
    lo = max(measured.v[0], *(s.v[0] for s in simulations.values()))
    hi = min(measured.v[-1], *(s.v[-1] for s in simulations.values()))
    voltage = np.arange(lo, hi, .1)
    current = np.interp(voltage, measured.v, measured.i)
    valid = current >= .3  # Existing whole-curve convention: exclude noisy near-Voc tail.
    plateau = valid & (voltage >= 5.) & (voltage <= 25.)
    rows = []
    for model, simulation in simulations.items():
        error = np.interp(voltage, simulation.v, simulation.i) - current
        row = dict(model=model, n_grid_pts=int(valid.sum()), plateau_n_grid_pts=int(plateau.sum()))
        for label, delta, reference in (("IV", error, measured.params["Isc_A"]),
                                         ("PV", voltage * error, measured.params["Pmax_W"])):
            row[f"nRMSE_{label}_pct"] = 100 * np.sqrt(np.mean(delta[valid] ** 2)) / reference
            row[f"nMBE_{label}_pct"] = 100 * np.mean(delta[valid]) / reference
        row["plateau_5_25V_RMSE_A"] = np.sqrt(np.mean(error[plateau] ** 2))
        row["plateau_5_25V_bias_A"] = np.mean(error[plateau])
        row["Pmax_measured_W"] = measured.params["Pmax_W"]
        row["Pmax_simulated_W"] = simulation.params["Pmax_W"]
        rows.append(row)
    return rows


def draw_curves(ax, measured, simulations, *, separate_models=False):
    color, marker = plots.tau_style(measured.tau_i)
    ax.plot(measured.v, measured.i, ls="none", marker=marker, ms=2.4,
            markevery=max(1, len(measured.v) // 28), color=color, mew=.35, zorder=4)
    # Draw the dashed constant model last so a nearly identical new curve cannot hide it.
    for model, style in (("pvsyst", "-"), ("constant", "--")):
        curve = simulations[model]
        line_color = (".35" if model == "constant" else "#0072B2") if separate_models else color
        ax.plot(curve.v, curve.i, ls=style, color=line_color,
                lw=.85 if model == "constant" else 1.25)
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 10)
    ax.set_xlabel("Voltage [V]")
    ax.set_ylabel("Current [A]")


def make_figures(states, output):
    model_handles = [Line2D([], [], color=".25", marker="o", ls="none", ms=3, label="Measured"),
                     Line2D([], [], color=".25", ls="--", label=r"Constant $R_{\rm sh}$"),
                     Line2D([], [], color=".25", ls="-", label=r"PVsyst $R_{\rm sh}(g)$")]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True, sharey=True)
        for row, shape in enumerate(("rectangle", "triangle")):
            for col, orientation in enumerate(("portrait", "landscape")):
                ax = axes[row, col]
                for key, (measured, simulations) in states.items():
                    if key[:2] == (shape, orientation):
                        draw_curves(ax, measured, simulations)
                ax.set_title(f"{shape.capitalize()}, {orientation}")
                ax.text(.03, .04, "abcd"[2 * row + col], transform=ax.transAxes, fontweight="bold")
        taus = sorted({key[2] for key in states}, reverse=True)
        tau_handles = [Line2D([], [], color=plots.tau_style(t)[0], lw=1.2,
                              label=rf"$\tau_I={t:g}$") for t in taus]
        fig.legend(handles=tau_handles, loc="upper center", ncol=len(taus), bbox_to_anchor=(.5, 1.02))
        fig.legend(handles=model_handles, loc="lower center", ncol=3)
        fig.tight_layout(rect=(0, .055, 1, .96))
        save(fig, output / "fig_yingli_shunt_comparison")

        for shape in ("rectangle", "triangle"):
            fig, axes = plt.subplots(2, len(taus), figsize=(10.5, 4.2), sharex=True, sharey=True)
            for row, orientation in enumerate(("portrait", "landscape")):
                for col, tau in enumerate(taus):
                    ax = axes[row, col]
                    state = states.get((shape, orientation, tau))
                    if state is None:
                        ax.text(.5, .5, "Not measured", transform=ax.transAxes, ha="center")
                    else:
                        draw_curves(ax, *state, separate_models=True)
                    if row == 0:
                        ax.set_title(rf"$\tau_I={tau:g}$")
                    if col == 0:
                        ax.set_ylabel(f"{orientation.capitalize()}\nCurrent [A]")
                    else:
                        ax.set_ylabel("")
                    if row == 0:
                        ax.set_xlabel("")
            fig.suptitle(f"Yingli: {shape}, {experiments.BAND_HEIGHT_MM:g} mm band", fontsize=9)
            detailed_handles = [model_handles[0], model_handles[1],
                                Line2D([], [], color="#0072B2", label=r"PVsyst $R_{\rm sh}(g)$")]
            fig.legend(handles=detailed_handles, loc="lower center", ncol=3)
            fig.tight_layout(rect=(0, .065, 1, .93))
            save(fig, output / f"fig_yingli_shunt_{shape}_by_transmittance")

        # Direct model difference exposes changes too small to resolve in a full-scale I-V panel.
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.5), sharex=True)
        for row, shape in enumerate(("rectangle", "triangle")):
            for col, orientation in enumerate(("portrait", "landscape")):
                ax = axes[row, col]
                for key, (measured, simulations) in states.items():
                    if key[:2] != (shape, orientation) or key[2] == 1.:
                        continue
                    constant, pvsyst = simulations["constant"], simulations["pvsyst"]
                    voltage = np.arange(0, min(constant.v[-1], pvsyst.v[-1]), .1)
                    delta = np.interp(voltage, pvsyst.v, pvsyst.i) - np.interp(voltage, constant.v, constant.i)
                    ax.plot(voltage, 1000 * delta, color=plots.tau_style(measured.tau_i)[0], lw=1.1)
                ax.axhline(0, color=".5", ls=":", lw=.6)
                ax.set_title(f"{shape.capitalize()}, {orientation}")
                ax.set_xlabel("Voltage [V]")
                ax.set_ylabel(r"$I_{\rm PVsyst}-I_{\rm constant}$ [mA]")
        fig.legend(handles=tau_handles[1:], loc="lower center", ncol=len(taus) - 1)
        fig.tight_layout(rect=(0, .055, 1, 1))
        save(fig, output / "fig_yingli_shunt_current_change")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=config.CH1_DIR / "yingli_shunt_comparison")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Retain a manifest of existing scientific exports and verify they stay untouched.
    protected = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in config.OUTPUT_DIR.rglob("*")
                 if p.is_file() and output not in p.parents}
    # This archived comparison opts in explicitly; the production model stays constant-Rsh.
    module = replace(modules.get("YINGLI"), rp_law=modules.PVsystShuntLaw(r_stc=5.99))
    models = {"constant": replace(module, rp_law=None), "pvsyst": module}
    simulators = {name: Simulator(m, output / "runs" / name) for name, m in models.items()}
    states, curve_rows, metric_rows, cell_rows, inputs = {}, [], [], [], {}
    for shape in ("rectangle", "triangle"):
        for orientation in ("portrait", "landscape"):
            measured = validation.load_measurements(module, shape, orientation)
            entries = experiments.flashes(module.name, shape, orientation)
            if len(measured) != len(entries):
                raise RuntimeError(f"Missing campaign measurements: {shape}, {orientation}")
            fcov = validation.band_fcov(module, shape, orientation)
            cells = geometry.cell_fcov(module, fcov, orientation)
            for measurement in measured:
                inputs[str(measurement.path)] = hashlib.sha256(measurement.path.read_bytes()).hexdigest()
                simulations = {}
                for model, simulator in simulators.items():
                    tag = f"YINGLI_{shape}_{orientation}_tau{measurement.tau_i:g}_{model}"
                    simulations[model] = (simulator.clean() if measurement.tau_i == 1.
                                          else simulator.run(cells, measurement.tau_i, tag))
                    print(f"{shape} {orientation} tau={measurement.tau_i:g}: {model} done", flush=True)
                key = (shape, orientation, measurement.tau_i)
                if measurement.tau_i == 1.:
                    np.testing.assert_array_equal(simulations["constant"].v, simulations["pvsyst"].v)
                    np.testing.assert_array_equal(simulations["constant"].i, simulations["pvsyst"].i)
                states[key] = (measurement, simulations)
                identifiers = dict(shape=shape, orientation=orientation, tau_I=measurement.tau_i)
                metric_rows.extend({**identifiers, **r} for r in paired_errors(measurement, simulations))
                for source, curve in {"measured": measurement, **simulations}.items():
                    curve_rows.append(pd.DataFrame({**identifiers, "source": source, "V": curve.v, "I": curve.i}))
                for cell, coverage in cells.items():
                    f = 0. if measurement.tau_i == 1. else round(coverage, 4)
                    cell_rows.append(dict(**identifiers, cell=cell, f_cov=f,
                                          g=1 - f * (1 - measurement.tau_i),
                                          Rsh_constant_ohm=module.rp_law.r_stc,
                                          Rsh_pvsyst_ohm=shaded_rp(module, f, measurement.tau_i)))
    metrics = pd.DataFrame(metric_rows)
    pd.concat(curve_rows, ignore_index=True).to_csv(output / "iv_curves.csv", index=False)
    metrics.to_csv(output / "curve_errors.csv", index=False)
    pd.DataFrame(cell_rows).to_csv(output / "cell_irradiance_shunt.csv", index=False)
    shaded = metrics[metrics.tau_I != 1.]
    summary = shaded.groupby(["shape", "orientation", "model"]).mean(numeric_only=True).reset_index()
    summary.to_csv(output / "error_summary.csv", index=False)
    make_figures(states, output)
    summary_text = summary[["shape", "orientation", "model", "nRMSE_IV_pct", "nRMSE_PV_pct",
                            "plateau_5_25V_RMSE_A"]].to_string(index=False)
    (output / "README.md").write_text(
        "# Yingli shunt-resistance comparison\n\n"
        "Reproduce: `MPLBACKEND=Agg python chapters/yingli_shunt_comparison.py`.\n\n"
        "Only the cell shunt resistance changes. The template, diode, series resistance, bypass diodes, "
        "photocurrent, temperatures and band geometry are identical between the two models. "
        "No Bishop breakdown term is introduced. Constant Rsh is 5.99 ohm/cell. The PVsyst law uses "
        "R0=4*5.99, exponent 5.5 and g=1-f_cov*(1-tau_I), including partial cell coverage.\n\n"
        "All 18 campaign states are plotted: 14 shaded measurements and the same clean flash repeated "
        "for four geometry/orientation comparisons. Transmittances are 1, .706, .536, .442 and 0; "
        ".442 was not measured in landscape. Both rectangular and triangular 80 mm bands are included. "
        "Measured shaded curves use the existing clean-flash temperature correction.\n\n"
        "`fig_yingli_shunt_comparison` is the overview; the two `by_transmittance` figures use blue "
        "solid curves for PVsyst and grey dashed curves for constant Rsh. `fig_yingli_shunt_current_change` "
        "shows the small direct current change in mA. PDF and PNG versions are supplied.\n\n"
        "`curve_errors.csv` compares both simulations on identical 0.1 V samples of their shared "
        "measured voltage support, excluding measured current below 0.3 A, as in the existing validation. "
        "IV and PV errors are normalized by measured Isc and Pmax. The 5-25 V current-window metric "
        "is predefined and includes transitions for triangular bands; it is not a fitted plateau region. "
        "`error_summary.csv` gives arithmetic means over shaded transmittances, excluding clean curves.\n\n"
        "The improvement is modest and does not resolve the portrait plateau discrepancy. IV nRMSE "
        "decreases for all 14 shaded curves, but PV error and maximum-power error do not improve "
        "uniformly. The clean simulations are identical. No parameters were fitted.\n\n"
        "Existing publication exports were preserved. They still describe the previous constant-Rsh "
        "simulations; this comparison does not regenerate chapter 2, energy or economic results.\n\n"
        "Mean shaded-curve errors:\n\n```text\n" + summary_text + "\n```\n")
    for p, digest in protected.items():
        assert hashlib.sha256(p.read_bytes()).hexdigest() == digest, f"Existing output modified: {p}"
    provenance = dict(module="YINGLI", law=asdict(module.rp_law),
                      formula="Rbase + (R0-Rbase)*exp(-5.5*g); R0=4*Rstc; g=1-f_cov*(1-tau_I)",
                      band_height_mm=experiments.BAND_HEIGHT_MM,
                      measurement_count=len(states), shaded_measurement_count=int(len(shaded) / 2),
                      model_temperature_C=module.t_ref, temperature_normalization=config.NORMALIZE_TO_REF_TEMPERATURE,
                      plateau_window_V=[5, 25], paired_voltage_step_V=.1, measured_current_cutoff_A=.3,
                      other_electrical_parameters_unchanged=True, bishop_breakdown_added=False,
                      clean_curves_identical=True,
                      preserved_existing_exports=len(protected), measurement_sha256=inputs,
                      netlist_sha256=hashlib.sha256(module.netlist_path.read_bytes()).hexdigest(),
                      source_sha256={p: hashlib.sha256((config.REPO_DIR / p).read_bytes()).hexdigest()
                                     for p in ("soilspice/spice.py", "soilspice/modules.py", "chapters/yingli_shunt_comparison.py")})
    (output / "verification.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(summary_text)


if __name__ == "__main__":
    main()
