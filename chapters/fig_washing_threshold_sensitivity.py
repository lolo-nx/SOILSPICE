"""Additional 3 × 3 threshold-sensitivity figures, without regenerating main figures.

    python chapters/fig_washing_threshold_sensitivity.py
"""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from soilspice import config, energy, energy_sensitivity as S
from soilspice import fig_threshold_sensitivity as F

STEMS = {"ground": "fig_washing_threshold_sensitivity_ground", "rooftop": "fig_washing_threshold_sensitivity_rooftop"}
KEY = ["site", "market", "module", "orientation", "R_H_mm_event"]
COLUMNS = KEY + ["shape", "h_mm", "f", "base_margin_mm", *F.METRICS,
    "energy_none_total_loss_pct", "energy_none_uniform_loss_pct", "clean_revenue_eur_kwp_year",
    "cost_none_total_eur_kwp_year", "cost_manual_total_eur_kwp_year", "manual_gain_eur_kwp_year",
    "manual_total_policy", "manual_total_interval_days", "manual_total_cleanings_per_year",
    "cost_opt_band_extra_eur_kwp_year", "margin10_savings_eur_kwp_year", "margin10_capex_eur_kwp",
    "margin10_net_eur_kwp_year", "margin10_opt_policy", "margin10_opt_interval_days",
    "opt_total_policy", "opt_uniform_policy"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = config.CH4_DIR
    out = config.ECON_DIR / "washing_threshold_sensitivity_3x3"
    out.mkdir(parents=True, exist_ok=True)
    # Protect all pre-existing main, supplementary and earlier sensitivity images.
    protected = {p: digest(p) for p in root.rglob("*") if p.suffix in (".pdf", ".png") and p.stem not in STEMS.values()}
    source = config.ECON_DIR / "threshold_sensitivity" / "summary_reference.csv"
    table = pd.read_csv(source)
    table = table[(table.h_mm == 60) & np.isclose(table.f, .1) & (table.base_margin_mm == 0) & (table['shape'] == 'rectangle')].copy()
    assert len(table) == 240 and not table.duplicated(KEY).any()
    np.testing.assert_allclose(energy.AF, (1 - 1.05 ** -25) / .05, atol=1e-12)
    source_runs, recomputed = [], []
    for threshold in S.THRESHOLDS:
        path = root / "total_scenarios" / f"PH{threshold:g}mm" / "data" / "model_runs.csv"
        runs = pd.read_csv(path)
        runs = runs[(runs.h == 60) & np.isclose(runs.f, .1) & runs.margin.isin((0, 10))]
        for market in ("ground", "rooftop"):
            policies = set(runs[(runs.market == market) & ~runs['mode'].isin(('none', 'robot'))]['mode'])
            assert policies == ({'trigger', 'interval'} if market == 'ground' else {'trigger'})
            if market == 'ground':
                assert set(runs[(runs.market == market) & (runs['mode'] == 'interval')].k) == set(energy.INTERVALS)
        recomputed.append(S.summarize_runs(runs, threshold))
        source_runs.append(dict(R_H_mm_event=threshold, path=str(path.relative_to(config.REPO_DIR)), sha256=digest(path)))
    expected = pd.concat(recomputed, ignore_index=True).set_index(KEY).sort_index()
    data = table[COLUMNS].set_index(KEY).sort_index()
    pd.testing.assert_index_equal(data.index, expected.index)
    expected = expected[data.columns]
    pd.testing.assert_frame_equal(data, expected, check_dtype=False, rtol=1e-11, atol=1e-11)
    data = data.reset_index()
    assert not data.manual_total_policy.isin(('none', 'robot')).any()
    np.testing.assert_allclose(data.manual_gain_pct_clean_revenue,
                              100 * (data.cost_none_total_eur_kwp_year - data.cost_manual_total_eur_kwp_year) / data.clean_revenue_eur_kwp_year, atol=1e-11)
    np.testing.assert_allclose(data.margin10_npv25_eur_kwp,
                              energy.AF * data.margin10_savings_eur_kwp_year - data.margin10_capex_eur_kwp, atol=1e-11)
    np.testing.assert_allclose(data.energy_none_band_extra_loss_pct,
                              data.energy_none_total_loss_pct - data.energy_none_uniform_loss_pct, atol=1e-11)
    data["illustrative_reference_R_H_mm_event"] = F.REFERENCE_THRESHOLD
    data.to_csv(out / "plotted_data.csv", index=False)
    layouts = []
    for market, stem in STEMS.items():
        q = data[data.market == market]
        q.to_csv(out / f"plotted_data_{market}.csv", index=False)
        layouts.append(F.draw(data, market, root / stem))
    # Conclusions retain signs and adjacent intervals from the same source tables.
    orientations = S.orientation_comparison(table)
    robust, _, transitions = S.classify_across_thresholds(table, orientations)
    transitions = transitions[transitions.conclusion.isin(("margin10", "manual_cleaning"))].copy()
    transitions.to_csv(out / "sign_changes.csv", index=False)
    robust.to_csv(out / "robustness_by_configuration.csv", index=False)
    captions = []
    for market in STEMS:
        policies = ("economic-trigger and fixed-interval manual policies (14, 30, 60, 90, 120, 180, 270 and 365 days)"
                    if market == 'ground' else "the economic-trigger manual policy, the only manual policy evaluated for this market")
        scales = ("Panels a, g use a separate Neuchâtel scale; b–c and h–i share scales, and d–f use a common scale."
                  if market == 'ground' else "Neuchâtel uses a separate vertical scale in each row; Almería and Doha share a scale within each row.")
        caption = (f"Sensitivity to the assumed washing threshold, {market} market. Rectangular band, h = 60 mm, f = 0.1, initial added margin 0 mm; all other assumptions fixed. "
                   "Columns: Neuchâtel, Almería and Doha. a–c, additional annual energy loss due to the band without manual cleaning: loss with band minus loss without band, as a percentage of identical clean production. "
                   "d–f, net present value of an additional 10 mm margin over 25 years at a 5% discount rate, in EUR/kWp, including the extra-glass investment; costs with and without band are independently optimised at each margin, as in the existing economic analysis. "
                   "g–i, annual net gain from manual cleaning for total soiling (uniform deposit plus band), as a percentage of clean revenue: total cost with no manual cleaning minus the minimum total cost of the evaluated manual policies, including cleaning charges. "
                   f"This selection considers {policies}, excluding both no manual cleaning and the robot; negative gains are retained. "
                   "Colours identify modules; solid lines/circles denote portrait and dashed lines/triangles landscape. "
                   "The dotted vertical line marks R_H = 15 mm/event, an illustrative reference specified for the main analysis, not a threshold chosen by optimising these curves. Horizontal lines mark zero economic benefit. "
                   f"{scales} Curves show individual configurations; no uncertainty bars or confidence intervals are used. "
                   "Complete removal of uniform and band deposits is assumed when event precipitation P_e >= R_H. Events are separated by three consecutive dry hours. These washing thresholds are not calibrated for bottom-edge bands. "
                   "The annual normalisation and available policy sets follow the saved simulations.")
        captions.append((market, caption))
        (out / f"caption_{market}.tex").write_text(caption.replace("%", r"\%").replace("R_H", r"$R_H$").replace("P_e", r"$P_e$").replace(">=", r"$\geq$") + "\n")
    lines = ["# Sensibilité au seuil de lavage — nouvelles figures 3 × 3", "",
             "Reproduction : `python chapters/fig_washing_threshold_sensitivity.py`. Les figures principales et les anciennes figures de sensibilité sont conservées.", "",
             "Les PDF vectoriels font exactement 183 × 170 mm ; les PNG sont exportés à 600 dpi. Une couleur par module et un style par orientation ; panneaux a–i, sans barre d’incertitude. Les hypothèses détaillées restent dans les légendes ci-dessous, hors de l’image.", "",
             "Données : plotted_data_ground.csv et plotted_data_rooftop.csv (120 lignes chacun), ou plotted_data.csv (240 lignes). Les trois indicateurs tracés et leurs composantes sont vérifiés par rapport à summary_reference.csv et recalculés depuis les cinq model_runs.csv. Les gains manuels négatifs sont conservés ; la ligne nettoyage concerne la salissure totale, pas uniquement la bande.", "",
             "Le scénario R_H=15 mm/événement est une référence illustrative demandée pour les figures principales. Aucun seuil n’est déterminé ici à partir de l’apparence des courbes ou d’un optimum économique.", "",
             "## Conclusions à paramètres fixés", "",
             "Nombre de configurations, sur 8, dont le bénéfice reste strictement positif aux cinq seuils :", "",
             "| Site | Marché | Marge +10 mm | Nettoyage manuel |", "|---|---|---|---|"]
    for (site, market), q in robust.groupby(['site', 'market'], sort=False):
        lines.append(f"| {site} | {market} | {q.margin10_all_profitable.sum()}/8 | {q.manual_cleaning_all_profitable.sum()}/8 |")
    lines += ["", "0/8 signifie qu’aucune configuration n’est rentable aux cinq seuils, pas qu’elle ne l’est jamais. L’énergie supplémentaire de bande reste positive pour les 8 configurations à chaque site et seuil. Pour les signes de marge et nettoyage, zero signifie absence de gain net, non un gain positif.", "",
              "### Intervalles de changement de signe", "",
              "Les intervalles ci-dessous encadrent uniquement des changements entre deux simulations adjacentes. Aucun seuil exact de bascule n’est interpolé. Les valeurs détaillées et les unités sont dans sign_changes.csv ; les signes aux cinq seuils pour chaque configuration sont dans robustness_by_configuration.csv.", "",
              "| Site / marché | Module / orientation | Indicateur | Intervalle R_H [mm/événement] | Signe avant → après |", "|---|---|---|---|---|"]
    for r in transitions.itertuples():
        lines.append(f"| {r.site} / {r.market} | {r.module} / {r.orientation} | {r.conclusion} | {r.R_H_low_mm:g}–{r.R_H_high_mm:g} | {r.state_low} → {r.state_high} |")
    for market, caption in captions:
        lines += ["", f"## Figure caption — {market}", "", caption]
    lines += ["", "[Trajectoires temporelles complémentaires](../../fig_band_growth_total_comparison.pdf) · [Conventions complètes des tableaux de synthèse](../threshold_sensitivity/README.md)", ""]
    (out / "README.md").write_text("\n".join(lines))
    unchanged = all(p.is_file() and digest(p) == sha for p, sha in protected.items())
    assert unchanged, 'A pre-existing figure was modified'
    verification = dict(figure_layouts=layouts, summary_table=str(source.relative_to(config.REPO_DIR)), summary_sha256=digest(source),
                        source_simulations=source_runs, rows_per_market=120, all_plotted_metrics_reconcile_with_raw_runs_and_summary=True,
                        manual_policy_excludes_none_and_robot=True, negative_cleaning_gains_retained=int((data.manual_gain_pct_clean_revenue < -1e-10).sum()),
                        economic_normalisation='same as existing tables: duration / 365.25 years; fractions relative to clean revenue',
                        reference_R_H_mm_event=15, reference_provenance='user-specified illustrative reference, not automatically selected',
                        preexisting_figure_count=len(protected), all_preexisting_figures_unchanged=unchanged)
    (out / "verification.json").write_text(json.dumps(verification, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(verification, indent=2, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
