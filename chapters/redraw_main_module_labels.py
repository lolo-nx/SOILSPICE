"""Refresh main figure labels from saved results without running simulations."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import pandas as pd
from soilspice import config, fig_validation, fig_simulation
from soilspice.style import ORDER


def main():
    out = config.CH1_DIR
    data = out / "data"
    for shape, name in (("rectangle", "iv_grid_rectangle"), ("triangle", "iv_grid_triangle_SI")):
        options = dict(figsize=(7.2, 3.15), compact_layout=True, panel_letters=False) if shape == "triangle" else {}
        fig_validation.iv_grid(pd.read_csv(data / f"iv_grid_{shape}.csv"), out / name, shape, **options)
    normalized = data / "iv_grid_rectangle_normalized.csv"
    if normalized.exists():
        fig_validation.iv_grid(pd.read_csv(normalized), out / "iv_grid_rectangle_normalized", normalized=True)
    fig_validation.errors_grid(pd.read_csv(config.VALIDATION_DIR / "curve_errors_summary.csv"), out / "errors_grid")
    fig_validation.model_hierarchy(pd.read_csv(data / "three_models_measured.csv"), out / "fig_model_hierarchy")

    out, data = config.CH2_DIR, config.SWEEP_DIR
    sweep = pd.read_csv(data / "same_acov_tauI_sweep.csv")
    states = pd.read_csv(data / "resilience_states.csv")
    states = states[states.module.isin(ORDER)]
    metrics = pd.read_csv(data / "resilience_metrics.csv")
    metrics = metrics[metrics.module.isin(ORDER)]
    fig_simulation.optical_loss_grid(sweep, out / "optical_loss_modules_grid")
    fig_simulation.regime_atlas(states, out / "regime_atlas")
    fig_simulation.average_power_windows_horizontal(metrics, out / "average_power_windows_horizontal")
    fig_simulation.initial_amplification(metrics, out / "initial_amplification")


if __name__ == "__main__":
    main()
