"""Chapter 2 - Electrical resilience across modules.

Paper: Electrical resilience across modules (optical_loss_modules_grid, regime_map, regime_atlas,
average_power_windows, initial_amplification), Methods: simulation framework and soiling scenarios. Also writes the
losses at equal band height used by chapters 3 and 4.

    python chapters/2_simulation.py            # from the saved loss maps (no LTspice)
    python chapters/2_simulation.py --full     # simulate every band with LTspice first (hours)

Steps
  1. loss maps at equal covered area: 4 shapes/orientations x A_cov 2-10 % x tau_I 0-0.9 (--full)
  2. electrical regimes and resilience metrics (M0, S*, P_an,S) of every module and configuration
  3. losses at equal band height h (margins included), every module (--full) -> equal_height_sweep.csv
  4. paper figures

Outputs: results/2_simulation/ (paper figures and table), results/2_simulation/data/ (loss maps, metrics, the
diagnostic figures of each module)
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from soilspice import analysis, config, fig_simulation, modules, plots, resilience, simulation, tables  # noqa: E402
from soilspice import uncertainty  # noqa: E402
from soilspice.style import ORDER  # noqa: E402

OUT = config.CH2_DIR
DATA = config.SWEEP_DIR
SWEEP_CSV = DATA / "same_acov_tauI_sweep.csv"


# ---- 1-2. loss maps at equal A_cov, regimes and resilience metrics ---------------------------------------------------
def loss_maps(full):
    if full:
        sweep = simulation.loss_maps([modules.get(n) for n in modules.MODULE_ORDER], configs=simulation.CONFIGS,
                                     acov_targets=simulation.ACOV_TARGETS, taus=simulation.TAU_VALUES)
        if sweep.empty:
            sys.exit("No simulation succeeded: check the netlists and LTSPICE_EXE in soilspice/config.py")
        sweep.to_csv(SWEEP_CSV, index=False)
    sweep = pd.read_csv(SWEEP_CSV)
    unc = uncertainty.global_model("Pmp")              # from chapter 1, None if not run
    sweep = analysis.add_regime(sweep)
    analysis.regime_census(sweep).to_csv(DATA / "regime_census.csv", index=False)
    analysis.resilience_summary(sweep).to_csv(DATA / "gamma_summary.csv", index=False)
    states = resilience.add_optical_loss(sweep)
    states.to_csv(DATA / "resilience_states.csv", index=False)
    metrics = resilience.summary(sweep)
    metrics.to_csv(DATA / "resilience_metrics.csv", index=False)
    cols = ["module", "shape", "orientation", "M0", "S_star_%", "max_loss_%", "P_an_S"]
    print("\nResilience metrics\n" + metrics[cols].round(2).to_string(index=False))
    tables.to_latex(metrics[cols].round(2), OUT / "table_resilience_metrics.tex",
                    caption="Resilience metrics of the SPICE loss maps: initial amplification $M_0$, "
                            "critical optical loss $S^*$, maximum power loss and $P_{an,S}$",
                    label="tab:resilience_metrics")
    # diagnostic figures of every module (data/)
    for (name, shape), _ in sweep.groupby(["module", "shape"]):
        plots.loss_map_grid(sweep, name, shape, DATA / f"heatmap_grid_{name}_{shape}", unc=unc)
    return sweep, states, metrics


# ---- 3. losses at equal band height ---------------------------------------------------------------------------------
def equal_height(full):
    out = config.EQUAL_HEIGHT_DIR
    out.mkdir(parents=True, exist_ok=True)
    if full:
        sweep = simulation.height_maps([modules.get(n) for n in modules.MODULE_ORDER], heights=simulation.HEIGHTS_MM)
        sweep = sweep.sort_values(["module", "shape", "orientation", "height_mm", "tau_I"]).reset_index(drop=True)
        sweep.to_csv(config.EQUAL_HEIGHT_CSV, index=False)
    sweep = pd.read_csv(config.EQUAL_HEIGHT_CSV)
    metrics = resilience.height_summary(sweep)
    metrics.to_csv(out / "equal_height_metrics.csv", index=False)
    unc = uncertainty.global_model("Pmp")
    plots.loss_vs_height(sweep, out / "loss_vs_height", taus=(0.5, 0.0), unc=unc)
    return sweep


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--full", action="store_true", help="simulate every band with LTspice first")
    full = ap.parse_args().full
    DATA.mkdir(parents=True, exist_ok=True)
    sweep, states, metrics = loss_maps(full)
    equal_height(full)
    states4, metrics4 = states[states["module"].isin(ORDER)], metrics[metrics["module"].isin(ORDER)]
    fig_simulation.optical_loss_grid(sweep, OUT / "optical_loss_modules_grid")
    fig_simulation.regime_map(states4, OUT / "regime_map")
    fig_simulation.regime_atlas(states4, OUT / "regime_atlas")
    fig_simulation.average_power_windows(metrics4, OUT / "average_power_windows")
    fig_simulation.initial_amplification(metrics4, OUT / "initial_amplification")
    combined = fig_simulation.resilience_overview(metrics4, OUT / "fig_resilience")
    combined.to_csv(DATA / "fig_resilience_data.csv", index=False)


if __name__ == "__main__":
    main()
