"""Instrument the original constant-Rsh Yingli portrait simulations; change no circuit parameters.

    MPLBACKEND=Agg python chapters/yingli_reverse_bias_diagnostic.py

Adds only .save requests to copies of the original netlists, then checks that
the module I-V curves reproduce the saved baseline. Internal branch currents
are read directly from LTspice rather than inferred from a fitted diode law.
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

from soilspice import config, experiments, geometry, modules, plots, validation
from soilspice.spice import run_ltspice
from soilspice.style import STYLE, save

POINTS = (0., .1, 1., 5., 10., 15., 20., 25., 30., 35.)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trace(raw, name):
    lookup = {n.lower(): n for n in raw.get_trace_names()}
    return np.asarray(raw.get_trace(lookup[name.lower()]).get_wave(), dtype=float)


def terminal(raw, name, voltage):
    return np.zeros_like(voltage) if name == "0" else trace(raw, f"V({name})")


def instrument(text, cells):
    """Only add output requests, leaving every original circuit line intact."""
    requests = ["I(V1)", "I(D1)", "I(D2)", "I(D3)"]
    for cell, positive, negative in cells:
        name = cell[1:].lower()  # LTspice XX10 becomes internal instance x10.
        requests += [f"V({node})" for node in (positive, negative) if node != "0"]
        requests += [f"V({name}:N001)", *(f"I({name}:{branch})" for branch in ("Rp", "D1", "I1", "Rs"))]
    line = ".save " + " ".join(dict.fromkeys(requests)) + "\n"
    assert text.count("\n.end\n") == 1
    output = text.replace("\n.end\n", "\n" + line + ".end\n")
    assert output.replace(line, "") == text
    return output


def make_figures(states, output):
    legend = [Line2D([], [], color=".25", marker="o", ls="none", ms=3, label="Measured"),
              Line2D([], [], color=".25", ls="--", label=r"SPICE, constant $R_{\rm sh}$")]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
        for tau, state in states.items():
            frame, measured = state
            mean = frame.groupby("V_module_V", sort=True).mean(numeric_only=True)
            v = mean.index.to_numpy()
            color, marker = plots.tau_style(tau)
            axes[0, 0].plot(measured.v, measured.i, ls="none", marker=marker, color=color,
                            ms=2.4, markevery=max(1, len(measured.v) // 25), mew=0)
            axes[0, 0].plot(v, mean.I_module_A, color=color, ls="--", lw=1)
            axes[0, 1].plot(v, mean.V_junction_V, color=color, lw=1.1)
            axes[1, 0].plot(v, -mean.I_shunt_A, color=color, lw=1.1)
            valid = (v >= measured.v[0]) & (v <= measured.v[-1])
            delta = mean.I_module_A.to_numpy()[valid] - np.interp(v[valid], measured.v, measured.i)
            axes[1, 1].plot(v[valid], delta, color=color, lw=1.1)
        axes[0, 0].set_ylabel("Module current [A]")
        axes[0, 0].legend(handles=legend, loc="lower left", fontsize=6)
        axes[0, 0].set_ylim(0, 10)
        axes[0, 1].set_ylabel("Shaded-cell junction voltage [V]")
        axes[0, 1].axhline(0, color=".4", ls=":", lw=.7)
        axes[0, 1].set_ylim(-6, 1)
        axes[1, 0].set_ylabel(r"Reverse shunt contribution $-I_{\rm sh}$ [A]")
        axes[1, 1].set_ylabel(r"$I_{\rm simulated}-I_{\rm measured}$ [A]")
        for k, ax in enumerate(axes.flat):
            ax.set_xlabel("Module voltage [V]")
            ax.set_xlim(0, 39)
            ax.axvspan(5, 25, color=".6", alpha=.07, zorder=0)
            ax.text(.02, .96, "abcd"[k], transform=ax.transAxes, va="top", fontweight="bold")
            if k >= 2:
                ax.axhline(0, color=".4", ls=":", lw=.7)
        taus = list(states)
        tau_legend = [Line2D([], [], color=plots.tau_style(t)[0], label=rf"$\tau_I={t:g}$") for t in taus]
        fig.legend(handles=tau_legend, loc="upper center", ncol=len(taus))
        fig.tight_layout(rect=(0, 0, 1, .94))
        save(fig, output / "fig_yingli_reverse_bias_diagnostic")

        fig, axes = plt.subplots(2, len(states), figsize=(10.5, 4.2), sharex=True, sharey="row")
        for col, (tau, (frame, measured)) in enumerate(states.items()):
            mean = frame.groupby("V_module_V", sort=True).mean(numeric_only=True)
            v = mean.index.to_numpy()
            color, marker = plots.tau_style(tau)
            ax = axes[0, col]
            ax.plot(measured.v, measured.i, ls="none", marker=marker, color=color, ms=2.5,
                    markevery=max(1, len(measured.v) // 28), mew=0)
            ax.plot(v, mean.I_module_A, color="#0072B2", ls="--", lw=1.1)
            ax.set_title(rf"$\tau_I={tau:g}$")
            ax.set_ylim(0, 10)
            ax = axes[1, col]
            ax.plot(v, mean.V_junction_V, color="#0072B2", lw=1.2)
            ax.plot(v, mean.V_terminal_V, color=".4", ls=":", lw=1.)
            ax.axhline(0, color=".5", lw=.6)
            samples = [point for point in POINTS if point <= v[-1]]
            ax.plot(samples, np.interp(samples, v, mean.V_junction_V), ls="none", marker="o", ms=2,
                    color="#0072B2")
            ax.set_ylim(-6, 1)
            for ax in axes[:, col]:
                ax.set_xlim(0, 39)
                ax.axvspan(5, 25, color=".6", alpha=.07, zorder=0)
            axes[1, col].set_xlabel("Module voltage [V]")
        axes[0, 0].set_ylabel("Module current [A]")
        axes[1, 0].set_ylabel("Shaded-cell voltage [V]")
        fig.legend(handles=[*legend,
                            Line2D([], [], color="#0072B2", label="Junction voltage"),
                            Line2D([], [], color=".4", ls=":", label="Cell terminal voltage")],
                   loc="lower center", ncol=4)
        fig.suptitle("Yingli, 80 mm rectangular band, portrait: original constant-shunt model", fontsize=9)
        fig.tight_layout(rect=(0, .08, 1, .94))
        save(fig, output / "fig_yingli_reverse_bias_by_transmittance")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=config.CH1_DIR / "yingli_shunt_comparison/runs/constant")
    parser.add_argument("--output", type=Path, default=config.CH1_DIR / "yingli_reverse_bias_diagnostic")
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    runs = output / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    module = replace(modules.get("YINGLI"), rp_law=None)
    protected = {path: digest(path) for path in config.OUTPUT_DIR.rglob("*")
                 if path.is_file() and output not in path.parents}
    protected.update({path: digest(path) for path in (config.REPO_DIR / "soilspice").glob("*.py")})
    protected[module.netlist_path] = digest(module.netlist_path)
    coverage = geometry.cell_fcov(module, validation.band_fcov(module, "rectangle", "portrait"), "portrait")
    most_shaded = [cell for cell, f in coverage.items() if abs(f - max(coverage.values())) < 1e-9]
    states, all_cells, points, summaries, manifests = {}, [], [], [], []
    max_kcl, max_ohm, max_series, max_iv_change = 0., 0., 0., 0.
    for measured in validation.load_measurements(module, "rectangle", "portrait"):
        tau = measured.tau_i
        stem = "clean_YINGLI" if tau == 1 else f"YINGLI_rectangle_portrait_tau{tau:g}_constant"
        net_path, raw_path = source / f"{stem}.net", source / f"{stem}.raw"
        text = net_path.read_text()
        assert not re.search(r"^XX.*\brp=", text, re.MULTILINE), "Expected constant-shunt source"
        cells = [(cell, *re.search(rf"^{cell}\s+(\S+)\s+(\S+)\s", text, re.MULTILINE).groups())
                 for cell in most_shaded]
        diagnostic_net = runs / f"{stem}_diagnostic.net"
        diagnostic_net.write_text(instrument(text, cells))
        raw = RawRead(str(run_ltspice(diagnostic_net)))
        baseline = RawRead(str(raw_path))
        v, current = trace(raw, "v1"), trace(raw, "I(V1)")
        np.testing.assert_array_equal(v, trace(baseline, "v1"))
        delta = float(np.max(np.abs(current - trace(baseline, "I(V1)"))))
        assert delta < 2e-5, f"Instrumented I-V changed: {delta} A"
        max_iv_change = max(max_iv_change, delta)
        keep = current >= 0
        bypass = np.stack([trace(raw, f"I(D{k})") for k in (1, 2, 3)])
        frames = []
        for cell, positive, negative in cells:
            prefix = cell[1:].lower()
            v_negative = terminal(raw, negative, v)
            v_terminal = terminal(raw, positive, v) - v_negative
            v_junction = trace(raw, f"V({prefix}:N001)") - v_negative
            iph, ish, diode, series = [trace(raw, f"I({prefix}:{branch})")
                                         for branch in ("I1", "Rp", "D1", "Rs")]
            max_kcl = max(max_kcl, float(np.max(np.abs(iph - ish - diode + series))))
            max_ohm = max(max_ohm, float(np.max(np.abs(v_junction - 5.99 * ish))))
            max_series = max(max_series, float(np.max(np.abs(v_terminal - v_junction - .00491 * series))))
            f = 0. if tau == 1 else round(coverage[cell], 4)
            frame = pd.DataFrame(dict(tau_I=tau, cell=cell, f_cov=f, g=1 - f * (1 - tau),
                                      V_module_V=v[keep], I_module_A=current[keep],
                                      V_terminal_V=v_terminal[keep], V_junction_V=v_junction[keep],
                                      I_ph_A=iph[keep], I_shunt_A=ish[keep], I_diode_A=diode[keep],
                                      I_cell_output_A=-series[keep],
                                      P_cell_absorbed_W=v_terminal[keep] * series[keep],
                                      bypass_current_max_abs_A=np.max(np.abs(bypass[:, keep]), axis=0)))
            frames.append(frame)
            for point in POINTS:
                if point > frame.V_module_V.max():
                    continue
                row = dict(tau_I=tau, cell=cell, V_module_V=point)
                for column in ("I_module_A", "V_terminal_V", "V_junction_V", "I_ph_A", "I_shunt_A", "I_diode_A",
                               "I_cell_output_A", "P_cell_absorbed_W", "bypass_current_max_abs_A"):
                    row[column] = float(np.interp(point, v, frame[column])) if keep.all() else float(
                        np.interp(point, v[keep], frame[column]))
                row["I_measured_A"] = (float(np.interp(point, measured.v, measured.i))
                                       if measured.v[0] <= point <= measured.v[-1] else np.nan)
                points.append(row)
        frame = pd.concat(frames, ignore_index=True)
        mean = frame.groupby("V_module_V", sort=True).mean(numeric_only=True)
        crossing = np.flatnonzero(np.diff(np.sign(mean.V_junction_V)) != 0)
        crossing_v = (float(np.interp(0, mean.V_junction_V.iloc[crossing[0]:crossing[0]+2],
                                     mean.index[crossing[0]:crossing[0]+2])) if len(crossing) else np.nan)
        row = dict(tau_I=tau, shaded_cell_count=len(cells), f_cov=0. if tau == 1 else round(max(coverage.values()), 4),
                   V_junction_min_V=frame.V_junction_V.min(), reverse_to_forward_module_voltage_V=crossing_v,
                   max_spread_between_shaded_cells_V=frame.groupby("V_module_V").V_junction_V.apply(np.ptp).max())
        for point in (0., 5., 10., 15., 20., 25.):
            for column in ("V_junction_V", "I_module_A", "I_ph_A", "I_shunt_A", "I_diode_A"):
                row[f"{column}_at_{point:g}V"] = float(np.interp(point, mean.index, mean[column]))
        grid = np.arange(5., 25.001, .1)
        residual = np.interp(grid, mean.index, mean.I_module_A) - np.interp(grid, measured.v, measured.i)
        row["mean_current_excess_5_25V_A"] = float(np.mean(residual))
        summaries.append(row)
        states[tau] = (frame, measured)
        all_cells.append(frame)
        manifests.append(dict(tau_I=tau, source_net=str(net_path), source_net_sha256=digest(net_path),
                              source_raw_sha256=digest(raw_path), measurement=str(measured.path),
                              measurement_sha256=digest(measured.path), only_save_directive_added=True,
                              max_module_current_change_A=delta))
        print(f"tau={tau:g}: {len(cells)} cells, minimum junction voltage {row['V_junction_min_V']:.6f} V", flush=True)
    assert max_kcl < 2e-5 and max_ohm < 2e-5 and max_series < 2e-5, (max_kcl, max_ohm, max_series)
    pd.concat(all_cells, ignore_index=True).to_csv(output / "cell_states.csv", index=False)
    point_table = pd.DataFrame(points)
    point_table.to_csv(output / "diagnostic_points_by_cell.csv", index=False)
    point_table.groupby(["tau_I", "V_module_V"]).mean(numeric_only=True).reset_index().to_csv(
        output / "diagnostic_points_mean.csv", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(output / "diagnostic_summary.csv", index=False)
    make_figures(states, output)
    for path, sha in protected.items():
        assert digest(path) == sha, f"Existing file changed: {path}"
    verification = dict(model="original constant Rsh, 5.99 ohm/cell", electrical_parameters_changed=False,
                        bishop_terms_added=False, most_shaded_cells=most_shaded, source_cases=manifests,
                        max_module_current_change_A=max_iv_change, max_KCL_residual_A=max_kcl,
                        max_shunt_voltage_residual_V=max_ohm, max_series_voltage_residual_V=max_series,
                        preserved_preexisting_files=len(protected),
                        source_sha256=digest(Path(__file__)),
                        signs="Vj=V(internal diode anode)-V(cell minus); Ish and Id positive anode-to-minus; "
                              "Iph positive minus-to-anode; Icell_out=Iph-Id-Ish")
    (output / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(summary[["tau_I", "V_junction_V_at_0V", "V_junction_V_at_5V", "V_junction_V_at_15V",
                   "V_junction_V_at_25V", "I_shunt_A_at_15V", "mean_current_excess_5_25V_A"]].to_string(index=False))


if __name__ == "__main__":
    main()
