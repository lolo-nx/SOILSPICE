"""Build two Extended Data figures from corrected saved simulations; run no new model.

    MPLBACKEND=Agg python chapters/extended_data_sensitivity.py
"""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from soilspice import config, energy, energy_sensitivity as S, rain
from soilspice import fig_extended_sensitivity as F
from soilspice.style import ORDER, ORIENTATIONS

KEY = ["site", "market", "module", "orientation", "R_H_mm_event"]
STEMS = ("fig_ed1_parameter_sensitivity", "fig_ed2_washing_threshold_sensitivity")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = config.REPO_DIR
    out = config.ECON_DIR / "extended_data_sensitivity"
    manifest_path = config.ECON_DIR / "corrected_timestamp_regeneration.json"
    provenance = json.loads(manifest_path.read_text())
    # Refuse stale pre-correction results or inputs changed since the full regeneration.
    for name, sha in provenance["model_source_sha256"].items():
        assert digest(root / "soilspice" / name) == sha, f"Model changed since corrected regeneration: {name}"
    economic_source = "soilspice/energy_sensitivity.py"
    assert digest(root / economic_source) == provenance["pipeline_source_sha256"][economic_source]
    for site, sha in provenance["climate_source_sha256"].items():
        assert digest(config.CLIMATE_DIR / energy.FILES[site]) == sha, f"Climate changed: {site}"
    assert digest(config.EQUAL_HEIGHT_CSV) == provenance["spice_sweep_sha256"]
    assert energy.REFERENCE_TOTAL == 15 and energy.LIGHT == 1 and rain.DRY_HOURS == 3
    assert energy.TOTAL_THRESHOLDS == S.THRESHOLDS == (5., 10., 15., 20., 25.)
    assert energy.MARKETS["ground"] == (.06, .41)
    np.testing.assert_allclose(energy.AF, (1 - 1.05 ** -25) / .05, atol=1e-12)
    for relative, sha in provenance["export_sha256"].items():
        assert digest(root / relative) == sha, f"Saved corrected export changed: {relative}"
    print("Corrected model, climate, SPICE and saved-export hashes verified.", flush=True)

    source_cases = config.CH4_DIR / "total_scenarios/PH15mm/data/sensitivity.csv"
    source_summary = config.ECON_DIR / "threshold_sensitivity/summary_reference.csv"
    cases = pd.read_csv(source_cases)
    expected_cases = {name for name, _ in energy.SENS_CASES}
    assert expected_cases == {"reference", *(f"{p} {v}" for p, _ in F.PARAMETERS for v in ("low", "high"))}
    assert len(cases) == 33 and not cases.duplicated(["site", "case"]).any()
    assert set(cases.site) == set(energy.SITES)
    for _, group in cases.groupby("site"):
        assert set(group.case) == expected_cases
    cases["R_H_mm_event"] = 15.
    cases["h_mm"], cases["shape"], cases["margin_mm"], cases["market"] = 60., "rectangle", 0., "ground"
    cases["configuration_count"] = 8
    # Record actual inputs, particularly alpha sensitivity at unchanged deposition D.
    settings = dict(energy.SENS_CASES)
    for i, r in cases.iterrows():
        p = settings[r["case"]]
        cases.loc[i, "f"] = p.get("f", .1)
        cases.loc[i, "soiling_rate_pct_day"] = energy.RATE[r.site] * p.get("rate", 1.)
        cases.loc[i, "D_g_m2_day"] = energy.deposition(r.site, p.get("rate", 1.))
        cases.loc[i, "alpha_m2_g"] = energy.ALPHA * p.get("alpha", 1.)
        cases.loc[i, "R_L_mm_event"] = p.get("light", energy.LIGHT)
        cases.loc[i, "electricity_eur_kwh"] = energy.MARKETS["ground"][0]
        cases.loc[i, "cleaning_eur_kwp_per_intervention"] = energy.MARKETS["ground"][1] * p.get("cost", 1.)
    for site, group in cases.groupby("site"):
        alpha_cases = group[group.case.str.startswith("alpha")]
        np.testing.assert_allclose(alpha_cases.D_g_m2_day, energy.deposition(site), atol=1e-14)
    cases["energy_none_band_extra_loss_pct_median"] = cases.band_loss
    cases["trigger_band_extra_cost_pct_clean_revenue_median"] = cases.band_cost

    summary = pd.read_csv(source_summary)
    data = summary[(summary.market == "ground") & (summary.h_mm == 60) & np.isclose(summary.f, .1) &
                   (summary.base_margin_mm == 0) & (summary["shape"] == "rectangle")].copy()
    assert len(data) == 120 and not data.duplicated(KEY).any()
    expected_configs = {(m, o) for m in ORDER for o in ORIENTATIONS}
    for _, group in data.groupby(["site", "R_H_mm_event"]):
        assert set(zip(group.module, group.orientation)) == expected_configs
    sources, raw_reference, recomputed = [], None, []
    for scenario in provenance["scenarios"]:
        threshold = scenario["R_H_mm_event"]
        path = config.CH4_DIR / "total_scenarios" / f"PH{threshold:g}mm/data/model_runs.csv"
        assert digest(path) == scenario["model_runs_sha256"]
        runs = pd.read_csv(path)
        runs = runs[(runs.market == "ground") & (runs.h == 60) & np.isclose(runs.f, .1) & runs.margin.isin((0, 10))]
        recomputed.append(S.summarize_runs(runs, threshold))
        if threshold == 15:
            raw_reference = runs[runs.margin == 0].copy()
        sources.append(dict(path=str(path.relative_to(root)), sha256=digest(path)))
    expected = pd.concat(recomputed, ignore_index=True).set_index(KEY).sort_index()
    observed = data.set_index(KEY).sort_index()[expected.columns]
    pd.testing.assert_frame_equal(observed, expected, check_dtype=False, atol=1e-10, rtol=1e-11)
    np.testing.assert_allclose(data.energy_none_band_extra_loss_pct,
                              data.energy_none_total_loss_pct - data.energy_none_uniform_loss_pct, atol=1e-11)
    np.testing.assert_allclose(data.margin10_npv25_eur_kwp,
                              energy.AF * data.margin10_savings_eur_kwp_year - data.margin10_capex_eur_kwp, atol=1e-11)
    # Manual-cleaning results are exported in full, including negative gains, but never plotted.
    assert not data.manual_total_policy.isin(("none", "robot")).any()
    np.testing.assert_allclose(data.manual_gain_pct_clean_revenue,
                              100 * (data.cost_none_total_eur_kwp_year - data.cost_manual_total_eur_kwp_year) /
                              data.clean_revenue_eur_kwp_year, atol=1e-11)
    assert set(data.study_start_utc) == {"2016-01-01 00:00:00+00:00"}
    assert set(data.study_stop_utc) == {provenance["source_coverage_stop_utc"]}
    assert (data.missing_trailing_source_hours == 1).all()

    metrics = [*F.METRICS, "manual_gain_pct_clean_revenue"]
    aggregates = data.groupby(["site", "R_H_mm_event"], sort=False)[metrics].agg(["min", "median", "max"])
    aggregates.columns = [f"{metric}_{stat}" for metric, stat in aggregates.columns]
    aggregates = aggregates.reset_index()
    aggregates["configuration_count"] = 8
    checks = []
    for site in energy.SITES:
        base = cases[(cases.site == site) & (cases.case == "reference")].iloc[0]
        threshold = aggregates[(aggregates.site == site) & (aggregates.R_H_mm_event == 15)].iloc[0]
        raw = raw_reference[raw_reference.site == site]
        none, trigger = raw[raw["mode"] == "none"], raw[raw["mode"] == "trigger"]
        assert len(none) == len(trigger) == 8
        raw_energy = np.median(100 * (none.loss - none.loss_u))
        raw_cost = np.median(100 * (trigger.cost - trigger.cost_u) / trigger.revenue)
        ed2_energy = threshold.energy_none_band_extra_loss_pct_median
        np.testing.assert_allclose([base.band_loss, ed2_energy], raw_energy, atol=1e-11, rtol=1e-11)
        np.testing.assert_allclose(base.band_cost, raw_cost, atol=1e-11, rtol=1e-11)
        checks.append(dict(site=site, R_H_mm_event=15, ed1_band_energy_pct=base.band_loss,
                           ed2_band_energy_pct=ed2_energy, raw_band_energy_pct=raw_energy,
                           difference_percentage_points=base.band_loss - ed2_energy,
                           ed1_trigger_band_cost_pct=base.band_cost, raw_trigger_band_cost_pct=raw_cost))
    checks = pd.DataFrame(checks)
    print("All figure metrics reconcile with raw policy runs and synthesis tables; RH=15 references agree.", flush=True)

    out.mkdir(parents=True, exist_ok=True)
    cases.to_csv(out / "ed1_parameter_cases.csv", index=False)
    data.to_csv(out / "ed2_configurations.csv", index=False)
    aggregates.to_csv(out / "ed2_plot_data.csv", index=False)
    checks.to_csv(out / "reference_checks.csv", index=False)
    layouts = [F.parameter_sensitivity(cases, config.CH4_DIR / STEMS[0]),
               F.washing_threshold_sensitivity(aggregates, config.CH4_DIR / STEMS[1])]
    captions = (
        (r"\caption{\textbf{Sensitivity to band-growth and cleaning parameters.} "
         r"Additional band energy loss with rain-only cleaning (a--c) and additional soiling cost "
         r"under the trigger rule (d--f), relative to uniform-soiling simulations. Bars show "
         r"medians over four modules and two orientations for independently varied lower (blue) "
         r"and higher (orange) values; black lines mark the reference. Ground-market costs, a "
         r"60~mm rectangular band, no added margin and $R_H=15$~mm/event are fixed. Reference "
         r"$f=0.1$ and $R_L=1$~mm/event; $D$ remains fixed when varying $\alpha$. Horizontal "
         r"scales differ between sites. Events end after three consecutive dry hours; complete "
         r"washing at $P_e\geq R_H$ is assumed and uncalibrated for bands."
         r"\label{fig:ed-parameter-sensitivity}}"),
        (r"\caption{\textbf{Sensitivity to the natural-washing threshold.} Additional band energy "
         r"loss with rain-only cleaning (a--c) and net present value of a 10~mm additional glass "
         r"margin (d--f), using ground-market costs, a 60~mm rectangular band, $f=0.1$ and no "
         r"initial added margin. Lines show medians and shading the minimum--maximum range of "
         r"eight configurations, not confidence intervals. Margin values include investment "
         r"over 25~years at 5\% discounting, with independent cleaning optimisations. Panel a "
         r"uses a separate scale; b--c and d--f share scales within each row. Dotted lines "
         r"mark the illustrative reference $R_H=15$~mm/event. Events end after three consecutive "
         r"dry hours; complete washing at $P_e\geq R_H$ is assumed, with thresholds uncalibrated "
         r"for bands.\label{fig:ed-washing-sensitivity}}"),
    )
    for i, caption in enumerate(captions, 1):
        (out / f"caption_ed{i}.tex").write_text(caption + "\n")
    # A lightweight manifest makes the existing-corrected-results provenance reviewable.
    import pymupdf
    from PIL import Image
    figure_checks = []
    for stem, layout in zip(STEMS, layouts):
        pdf, png = [config.CH4_DIR / f"{stem}.{suffix}" for suffix in ("pdf", "png")]
        with pymupdf.open(pdf) as doc:
            assert len(doc) == 1 and not doc[0].get_images()
            np.testing.assert_allclose([doc[0].rect.width, doc[0].rect.height],
                                      [F.WIDTH_MM / 25.4 * 72, F.HEIGHT_MM / 25.4 * 72], atol=1e-3)
            assert len(doc[0].get_fonts()) > 0
        with Image.open(png) as image:
            assert min(image.info["dpi"]) > 599
            layout["png_pixels"] = list(image.size)
        figure_checks.append(dict(stem=stem, **layout, vector_pdf=True,
                                  pdf_sha256=digest(pdf), png_sha256=digest(png)))
    # Previous detailed figures and numeric results are all preserved.
    for relative, sha in provenance["export_sha256"].items():
        assert digest(root / relative) == sha, f"Pre-existing export was modified: {relative}"
    verification = dict(source_manifest=str(manifest_path.relative_to(root)),
                        source_manifest_sha256=digest(manifest_path),
                        source_sensitivity_sha256=digest(source_cases), source_summary_sha256=digest(source_summary),
                        model_source_sha256=provenance["model_source_sha256"], source_simulations=sources,
                        generator_source_sha256={str(p.relative_to(root)): digest(p) for p in
                                                 (Path(__file__).resolve(), root / "soilspice/fig_extended_sensitivity.py")},
                        corrected_timestamp_alignment_verified=True, event_separation_dry_hours=3,
                        washing_criterion="completed event total >= R_H, mm/event; no intensity criterion",
                        reference_R_H_mm_event=15, thresholds_mm_event=list(energy.TOTAL_THRESHOLDS),
                        existing_results_reused_without_new_simulations=True,
                        alpha_sensitivity_preserves_deposition=True, independent_economic_optimisations_preserved=True,
                        all_metrics_reconcile_with_summary_and_raw_runs=True,
                        max_reference_difference_pct_points=float(checks.difference_percentage_points.abs().max()),
                        manual_results_retained_without_clipping=True,
                        preserved_preexisting_exports=len(provenance["export_sha256"]), figures=figure_checks,
                        inconsistencies=[])
    (out / "verification.json").write_text(json.dumps(verification, indent=2, ensure_ascii=False) + "\n")
    readme = """# Extended Data sensitivity figures

Reproduce: `MPLBACKEND=Agg python chapters/extended_data_sensitivity.py`.
This script reads existing corrected simulations; it runs no new physical model.
Figures: `../../fig_ed1_parameter_sensitivity.pdf` and
`../../fig_ed2_washing_threshold_sensitivity.pdf`, with matching 600 dpi PNGs.
Both vector PDFs have exact dimensions 183 × 125 mm, editable embedded fonts and six panels.
The earlier figures, source results and manuscript are unchanged.

## Source data and conventions

- `ed1_parameter_cases.csv`: all 33 site/case medians used for Fig. 1, including
  reference values, actual inputs and units. The five parameters vary one at a time
  at RH=15 mm/event. Alpha changes with D held fixed; a change of soiling rate
  changes D according to the existing mass conversion. Blue/orange bars terminate
  at the lower/higher case median; both extend from the same black reference line.
- `ed2_configurations.csv`: all 120 ground-market configuration/threshold rows,
  including the full original manual-cleaning results and signed gains, whether
  or not plotted. No policy is reselected by this figure script.
- `ed2_plot_data.csv`: all 15 site/threshold groups, each with exactly eight
  configurations. Minimum, median and maximum are calculated from the individual
  configuration metrics, not from differences of separately aggregated quantities.
  The manual-cleaning aggregates remain in this file but are not plotted.
- `reference_checks.csv`: independent comparisons of Fig. 1, Fig. 2 and raw policy
  runs at RH=15. The energy metric is comparable; the two economic metrics have
  different definitions and units and must not be equated.
- `caption_ed1.tex` and `caption_ed2.tex`: standalone LaTeX captions.
- `verification.json`: source hashes, data reconciliation and output checks.

Reference geometry: 60 mm rectangular band, f=0.1, no initial added glass margin.
Ground-market prices: 0.06 EUR/kWh and 0.41 EUR/kWp per manual intervention.
Energy percentages are extra lost energy with band minus without band under
natural rainfall only, divided by the identical clean production. They are not
total-soiling losses or percentages of soiled production.
Fig. 1 cost percentages compare the two independently simulated trigger policies
with and without band, including cleaning charges, divided by clean revenue.
They do not minimise over fixed-interval strategies or include a robot.
Fig. 2 margin NPV uses the existing independently optimised costs with/without
band at each margin. It equals the 25-year annuity factor at 5% times the annual
reduction in additional band cost, less the added-glass investment
(100 EUR/kWp × 10 mm / module length along the slope). The existing choice among
no manual cleaning, trigger and fixed intervals is preserved; robot excluded.
The negative values of net benefit or manual-cleaning gain are never clipped.

## Timing and washing

The corrected regeneration manifest is checked against the current rain/energy
code, all weather and SPICE inputs, and the source exports before plotting.
Open-Meteo hour-end labels are shifted once to interval starts after PV yield
calculation. Three consecutive dry hours close an event, including across
midnight/year boundaries; one rain action occurs at closure. Decisions use only
the completed event total, never intensity: total >= RH assumes full removal,
otherwise total >= RL transfers dust; smaller totals have no rain action.
UTC daily manual decisions and all cleaning assumptions remain unchanged.
The physical source coverage ends at 2025-12-31 23:00 UTC: the last source hour
is absent and no synthetic weather is added. Existing annualisation by actual
duration / 365.25 and clean-yield ratios are retained.

RH=15 mm/event is an illustrative reference, not selected from these curves.
Complete washing is an assumption; none of the washing thresholds is calibrated
for actual bottom-edge bands. Shaded ranges describe variation across the eight
configurations, not confidence intervals. Lines connect only the five simulated
thresholds; no additional scenario or exact crossing threshold is inferred.
Fig. 1 horizontal scales vary by site. Fig. 2 panel a has a separate energy scale;
panels b–c share a scale, and d–f share a scale.

All comparable RH=15 references and all raw-run/synthesis-table identities agree
within floating-point tolerance (1e-11 percentage points). No inconsistency found.
"""
    (out / "README.md").write_text(readme)
    print(checks.to_string(index=False), flush=True)
    print(f"Saved both Extended Data figures and source data to {out}", flush=True)


if __name__ == "__main__":
    main()
