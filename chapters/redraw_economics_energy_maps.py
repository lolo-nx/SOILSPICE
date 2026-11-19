"""Redraw PH15mm energy maps, with a Neuchatel inset, using saved results."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import pandas as pd
from soilspice import config, fig_energy

if __name__ == "__main__":
    scenario = config.CH4_DIR / "total_scenarios/PH15mm"
    maps = pd.read_csv(scenario / "data/energy_maps.csv")
    fig_energy.energy_maps(maps, scenario / "fig_energy_maps", neuchatel_zoom=True)
