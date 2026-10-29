"""Chapter 1 - Experimental validation of the model.

Paper: Experimental validation of the model (Fig. iv_grid_rectangle, errors_grid, Table uncertainty),
Uniform soiling against cell-level modelling (fig_model_hierarchy), Methods: modules and measurements
(pv_module_schematics). SI: iv_grid_triangle.

    python chapters/1_validation.py            # from the saved simulations (no LTspice)
    python chapters/1_validation.py --full     # simulate every flash test with LTspice first

Steps
  1. SPICE against the flash tests, per module (--full): I-V parameters, whole-curve errors, loss residuals
     -> data/validation/<MODULE>/, curve_errors_summary.csv, residuals.csv
  2. I-V curves of the four modules, measured and simulated (--full: LTspice; else the saved curves)
  3. whole-curve errors (nRMSE on P-V, nMBE on I-V)
  4. model uncertainty from the loss residuals: bias and 95 % interval (Table uncertainty)
  5. three model levels on the 60 shaded flashes: uniform soiling (L = S), ideal bypass, SPICE
  6. electrical layout of the modules (Methods)

Outputs: results/1_validation/ (paper figures and tables), results/1_validation/data/ (everything else)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from soilspice import (analysis, config, experiments, fig_validation, modules, plots, schematics,  # noqa: E402
                       tables, uncertainty, validation)
from soilspice.bypass import bypass_loss  # noqa: E402
from soilspice.simulation import CONFIGS  # noqa: E402
from soilspice.spice import Simulator  # noqa: E402
from soilspice.style import ORDER  # noqa: E402

OUT = config.CH1_DIR
DATA = config.CH1_DIR / "data"
TOLERANCE_PCT = 5.0


# ---- 1. SPICE against the flash tests ----------------------------------------------------------------------------
def validate_module(module):
    """Measured and simulated curves of every configuration of one module, with their errors."""
    out = module.output_dir
    out.mkdir(parents=True, exist_ok=True)
    sim = Simulator.for_module(module)
    results = [validation.run_configuration(module, shape, orientation, sim=sim, with_spice=True)
               for shape, orientation in CONFIGS]
    results = [r for r in results if r.measured]
    for r in results:
        tag = f"{r.shape}_{r.orientation}"
        plots.iv_comparison(r, out / f"iv_{tag}")
        losses = validation.loss_comparison([r])
        if not losses.empty:
            plots.loss_bars(losses, f"{module.label} - {r.shape}, {r.orientation}", out / f"losses_{tag}")
    plots.power_ratio(validation.power_ratio(module, CONFIGS), f"{module.label} - measured power ratio",
                      out / "power_ratio")
    errors = validation.parameter_errors(results)
    if errors.empty:
        return
    errors.to_csv(out / f"validation_detail_{module.name}.csv", index=False)
    summary = validation.summarize(errors, TOLERANCE_PCT)
    summary.to_csv(out / f"validation_summary_{module.name}.csv", index=False)
    print(f"\n{module.name}: relative error of the I-V parameters\n"
          + summary[summary["parameter"].isin(["Pmp", "all"])].round(2).to_string(index=False))
    pmp = summary[summary["parameter"].isin(["Pmp", "Vmp", "Imp"])]
    wide = pmp.pivot_table(index=["shape", "orientation"], columns="parameter", values=["MAE_%", "RMSE_%"])
    wide.columns = [f"{p} {stat.rstrip('_%')} [%]" for stat, p in wide.columns]
    tables.to_latex(wide.reset_index().round(2), out / f"mae_rmse_{module.name}.tex",
                    caption=f"MAE and RMSE between SPICE and measurement - {module.label}",
                    label=f"tab:mae_rmse_{module.name.lower()}")
    validation.whole_curve_errors(results).to_csv(out / f"curve_errors_{module.name}.csv", index=False)
    plots.error_bars(errors, f"{module.label} - SPICE vs measurement", out / f"error_bars_{module.name}")
    plots.parity(validation.loss_comparison(results), f"{module.label} - SPICE vs measured losses",
                 out / f"parity_{module.name}", tol_pct=TOLERANCE_PCT)
    res = uncertainty.residuals(results)
    if res is not None and not res.empty:
        res.to_csv(out / f"residuals_{module.name}.csv", index=False)


def _collect(pattern):
    """Concatenates the per-module files `pattern` of every module on disk."""
    frames = [pd.read_csv(p) for name in modules.MODULE_ORDER
              if (p := modules.get(name).output_dir / pattern.format(name=name)).exists()]
    frames = [f for f in frames if not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def update_summaries():
    """residuals.csv and curve_errors_summary.csv from the latest results of every module."""
    residuals = _collect("residuals_{name}.csv")
    if not residuals.empty:
        residuals.to_csv(config.VALIDATION_DIR / "residuals.csv", index=False)
    curves = _collect("curve_errors_{name}.csv")
    if not curves.empty:
        summary = curves.groupby(["module", "shape", "orientation"]).mean(numeric_only=True).reset_index()
        summary.to_csv(config.VALIDATION_DIR / "curve_errors_summary.csv", index=False)
        plots.curve_error_summary(curves, config.VALIDATION_DIR / "curve_errors_summary")


def compare_modules(selected):
    """Measured losses of all modules together, against tau_I and tau_eff, and the resilience slope gamma."""
    out = config.CROSS_MODULE_DIR
    out.mkdir(parents=True, exist_ok=True)
    for shape, orientation in CONFIGS:
        df = validation.measured_losses(selected, shape, orientation, with_spice=True)
        if df.empty:
            continue
        base = out / f"cross_module_loss_{shape}_{orientation}_EXPvsSPICE"
        df.to_csv(base.with_suffix(".csv"), index=False)
        plots.losses_vs_x(df, "tau_I", "$\\tau_I$", f"Measured losses - {shape}, {orientation}", f"{base}_vs_tauI")
        x_label = ("$\\tau_{eff}$ (cell, $=1-f_{cov}(1-\\tau_I)$)" if shape == "rectangle"
                   else "$\\tau_{eff}$ (most shaded sub-string)")
        plots.losses_vs_x(df, "tau_eff", x_label, f"Measured losses - {shape}, {orientation}", f"{base}_vs_taueff")
    for shape in ("rectangle", "triangle"):
        frames = [df.assign(shape=shape, orientation=o) for o in ("portrait", "landscape")
                  if not (df := validation.measured_losses(selected, shape, o)).empty]
        if frames:
            analysis.gamma_table(pd.concat(frames), "tau_eff").to_csv(out / f"gamma_{shape}.csv", index=False)


def run_validation():
    selected = [modules.get(n) for n in experiments.measured_modules(modules.MODULE_ORDER)]
    for module in selected:
        print(f"\n=== {module.label} ===")
        try:
            validate_module(module)
        except FileNotFoundError as err:
            print(f"  [!] {err}")
    update_summaries()
    compare_modules(selected)


# ---- 4. model uncertainty -----------------------------------------------------------------------------------------
def run_uncertainty(metric="Pmp"):
    out = config.UNCERTAINTY_DIR
    out.mkdir(parents=True, exist_ok=True)
    res = pd.read_csv(config.VALIDATION_DIR / "residuals.csv")
    res = res[res[f"e_{metric}"].notna()]
    glob = uncertainty.summary(res, metric)
    per_config = uncertainty.summary(res, metric, by=["module", "shape", "orientation"])
    summary = pd.concat([per_config, glob.assign(module="all", shape="", orientation="").drop(columns="group")],
                        ignore_index=True)
    summary.to_csv(OUT / "table_uncertainty.csv", index=False)
    tables.to_latex(summary[["module", "shape", "orientation", "n", "MBE", "RMSE", "sigma", "half_width"]].round(2),
                    OUT / "table_uncertainty.tex",
                    caption=f"Residuals of the {metric} loss, SPICE minus measurement [\\% points]. "
                            r"half\_width = 1.96 $\sigma$.",
                    label=f"tab:residuals_{metric.lower()}")
    print(f"\n{metric} loss residuals, SPICE - EXP [% points], n = {len(res)}\n"
          + glob.drop(columns="group").round(3).to_string(index=False))
    homog = uncertainty.homogeneity(res, metric)
    homog.to_csv(out / f"homogeneity_{metric}.csv", index=False)
    uncertainty.trends(res, metric).to_csv(out / f"trends_{metric}.csv", index=False)
    plots.residual_diagnostics(res, metric, out / f"residual_diagnostics_{metric}")
    plots.residual_overview(res, metric, out / f"residual_overview_{metric}")
    unc = uncertainty.model(res, metric)
    unc.to_csv(out / f"model_{metric}.csv", index=False)
    sweep_path = config.SWEEP_DIR / "same_acov_tauI_sweep.csv"            # loss maps of chapter 2, if there
    if sweep_path.exists():
        sweep = uncertainty.apply(pd.read_csv(sweep_path), unc, metric)
        sweep.to_csv(out / "same_acov_tauI_sweep_uncertainty.csv", index=False)
        for module in plots.sort_modules(sweep["module"]):
            for shape in ("rectangle", "triangle"):
                plots.loss_with_uncertainty(sweep, metric, module, shape, out / f"loss_uncertainty_{module}_{shape}")


# ---- 5. three model levels ------------------------------------------------------------------------------------------
def three_models():
    """Loss of uniform soiling (L = S), of an ideal-bypass model and of SPICE, on the 60 shaded flashes (h = 80 mm)
    and on the equal-height sweep of chapter 2."""
    r = pd.read_csv(config.VALIDATION_DIR / "residuals.csv")
    r["S"] = 100 * r.A_cov * (1 - r.tau_I)
    r["L_meas"] = r["Pmp_loss_exp_%"].abs(); r["L_spice"] = r["Pmp_loss_sim_%"].abs()
    r["L_bypass"] = [100 * bypass_loss(m, s, o, 80.0, t) for m, s, o, t in zip(r.module, r["shape"], r.orientation, r.tau_I)]
    r.to_csv(DATA / "three_models_measured.csv", index=False)
    for c in ("S", "L_bypass", "L_spice"):
        e = r[c] - r.L_meas
        print(f"{c:9s} bias {e.mean():+.1f}  RMSE {np.sqrt((e ** 2).mean()):.1f}  max {e.abs().max():.1f} pp")
    if config.EQUAL_HEIGHT_CSV.exists():
        sw = pd.read_csv(config.EQUAL_HEIGHT_CSV)
        sw = sw[sw.module.isin(ORDER)].copy()
        sw["S"] = 100 * sw.A_cov * (1 - sw.tau_I); sw["L_spice"] = sw["Pmp_loss_%"].abs()
        sw["L_bypass"] = [100 * bypass_loss(m, s, o, h, t) for m, s, o, h, t in
                          zip(sw.module, sw["shape"], sw.orientation, sw.height_mm, sw.tau_I)]
        sw.to_csv(DATA / "three_models_sweep.csv", index=False)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--full", action="store_true", help="simulate the flash tests with LTspice first")
    full = ap.parse_args().full
    DATA.mkdir(parents=True, exist_ok=True)
    if not full:
        # Check before writing any figures: the saved curves and residuals may
        # still use the old per-module transmittances from an earlier campaign.
        try:
            for shape in ("rectangle", "triangle"):
                csv = DATA / f"iv_grid_{shape}.csv"
                if csv.exists():
                    validation.check_saved_campaign(
                        pd.read_csv(csv, usecols=["module", "orientation", "tau_I"]), shape=shape)
            residuals = config.VALIDATION_DIR / "residuals.csv"
            if residuals.exists():
                validation.check_saved_campaign(
                    pd.read_csv(residuals, usecols=["module", "shape", "orientation", "tau_I"]),
                    include_clean=False)
        except ValueError as err:
            ap.error(str(err))
    if full:
        run_validation()
    for shape, name in (("rectangle", "iv_grid_rectangle"), ("triangle", "iv_grid_triangle_SI")):
        csv = DATA / f"iv_grid_{shape}.csv"
        if full or not csv.exists():
            fig_validation.collect(shape).to_csv(csv, index=False)
        fig_validation.iv_grid(pd.read_csv(csv), OUT / name, shape)
    fig_validation.errors_grid(pd.read_csv(config.VALIDATION_DIR / "curve_errors_summary.csv"), OUT / "errors_grid")
    run_uncertainty()
    fig_validation.model_hierarchy(three_models(), OUT / "fig_model_hierarchy")
    schematics.draw(OUT / "pv_module_schematics")


if __name__ == "__main__":
    main()
