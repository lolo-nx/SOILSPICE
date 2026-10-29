"""Chapter 4 - Energy yield and cleaning costs.

Paper: Three climates, three ways a band grows (fig_band_growth), Energy lost to the band (fig_energy_maps),
Lever 1 margin (fig_design_height, fig_design_levers), Lever 2 orientation, Lever 3 cleaning (fig_cleaning,
Table clean), Inverse diagnosis (fig_inverse_chain), rectangular against triangular band (fig_band_shape), Methods: sites (Table sites), band growth model, sensitivity
(Extended Data fig_sensitivity). SI: fig_pan_h_maps. No LTspice: uses the equal-height sweep of chapter 2.

    python chapters/4_economics.py                     # reference R_H=15 mm/event
    python chapters/4_economics.py --all-thresholds    # all five saved scenarios
    python chapters/4_economics.py --comparison-only   # compare all five thresholds
    python chapters/4_economics.py --washing-total 10   # one sensitivity scenario
    python chapters/4_economics.py --band-growth-only   # reference trajectories

Steps
  1. site table: yield, rain, rainy days, soiling rate and dust deposition (ERA5 2016-2025, pvlib)
  2. band growth model on every module, orientation, band height h, trapped fraction f, glass margin and cleaning
     strategy, both markets (soilspice.energy)
  3. maps of the energy lost on the full f grid; event-based series of the reference module (JA Solar, landscape)
  4. one-at-a-time sensitivity; band height at 1 % loss; smallest loss for which one cleaning pays
  5. paper figures, and the numbers quoted in the text

Outputs: results/4_economics/ (comparisons and shared data), total_scenarios/PH{value}mm/.
The user-selected illustrative reference is PH15mm; other thresholds are sensitivity scenarios.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pandas as pd  # noqa: E402

from soilspice import config, energy, fig_energy, tables, rain_scenarios  # noqa: E402


def run_economics(threshold, *, band_growth_only=False):
    """All main figures in the event-total scenario folder."""
    OUT = config.CH4_DIR / "total_scenarios" / rain_scenarios.scenario_tag(threshold)
    DATA = OUT / "data"
    DATA.mkdir(parents=True, exist_ok=True)
    rain_scenarios.write_scenario_settings(threshold, OUT)
    rain_kw = rain_scenarios.options(threshold)
    if band_growth_only:
        series = energy.growth_series(**rain_kw)
        fig_energy.band_growth(series, OUT / "fig_band_growth", threshold=threshold)
        return
    sites = energy.site_table(**rain_kw)
    sites.round(3).to_csv(OUT / "table_sites.csv", index=False)
    print(sites.round(2).to_string(index=False))

    runs = energy.model_runs(**rain_kw); runs.to_csv(DATA / "model_runs.csv", index=False)
    maps = energy.tau_grid(**rain_kw);   maps.to_csv(DATA / "energy_maps.csv", index=False)
    series = energy.growth_series(**rain_kw); series.to_csv(DATA / "growth_series_reference.csv", index=False)
    # Refresh an existing compatibility filename too; its rows are event/midnight intervals, not daily samples.
    if (DATA / "daily_series_reference.csv").exists():
        series.to_csv(DATA / "daily_series_reference.csv", index=False)
    sens = energy.sensitivity(**rain_kw); sens.to_csv(DATA / "sensitivity.csv", index=False)
    h1 = energy.design_h1();            h1.to_csv(DATA / "design_h1.csv", index=False)
    bound = energy.cleaning_bound(**rain_kw); bound.to_csv(OUT / "table_clean_bound.csv", index=False)
    tables.sites_table(sites, bound, OUT / "table_sites.tex", **rain_kw)
    tables.band_params_table(OUT / "table_band_parameters.tex", **rain_kw)

    numbers = fig_energy.paper_figures(runs, maps, series, sens, h1, OUT, DATA, **rain_kw)
    numbers.to_csv(OUT / "numbers_in_text.csv", index=False)

    # rectangle against triangle (f = 0.1)
    tri = energy.tau_grid("triangle", fs=[0.1], **rain_kw); tri.to_csv(DATA / "energy_triangle_f0.1.csv", index=False)
    fig_energy.band_shape(pd.concat([maps, tri]), OUT / "fig_band_shape")

    # Cycle-based diagnostics use the same rain-event criterion as the mass model.
    e_cycle = energy.washing_cycles(**rain_kw)
    cp = energy.cleaning_kwh()
    fig_energy.inverse_chain(OUT / "fig_inverse_chain", config.INPUT_DIR / "inverse_example" / "inverse_example_iv.csv",
                             e_cycle, cp)
    pan = energy.pan_h();               pan.to_csv(DATA / "pan_h.csv", index=False)
    fig_energy.pan_h_maps(pan, cp, e_cycle, OUT / "fig_pan_h_maps_SI")


def main():
    ap = argparse.ArgumentParser(description="Event-total washing: reference R_H=15 mm/event and sensitivity at 5,10,20,25 mm/event.")
    selected = ap.add_mutually_exclusive_group()
    selected.add_argument("--washing-total", type=float, choices=energy.TOTAL_THRESHOLDS,
                          help="run one scenario (default: 15 mm/event)")
    selected.add_argument("--all-thresholds", action="store_true",
                          help="run the five thresholds and their comparison")
    ap.add_argument("--comparison-only", action="store_true",
                    help="export the comparison of all five thresholds without full economics figures")
    ap.add_argument("--band-growth-only", action="store_true", help="draw only the requested growth trajectories")
    args = ap.parse_args()
    if args.comparison_only and (args.washing_total is not None or args.band_growth_only or args.all_thresholds):
        ap.error("--comparison-only compares all five thresholds and cannot be combined with another run mode")
    config.ECON_DIR.mkdir(parents=True, exist_ok=True)
    if args.comparison_only or (args.all_thresholds and not args.band_growth_only):
        rain_scenarios.run_comparison()
    if args.comparison_only:
        return
    thresholds = energy.TOTAL_THRESHOLDS if args.all_thresholds else (
        args.washing_total if args.washing_total is not None else energy.REFERENCE_TOTAL,)
    for threshold in thresholds:
        print(f"\nEconomics event-total scenario: {threshold:g} mm/event", flush=True)
        run_economics(threshold, band_growth_only=args.band_growth_only)
    if args.all_thresholds and not args.band_growth_only:
        rain_scenarios.finish_economics_exports()


if __name__ == "__main__":
    main()
