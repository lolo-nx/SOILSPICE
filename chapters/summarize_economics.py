"""Summarise saved energy/economic sensitivity to R_H=5,10,15,20,25 mm/event.

    python chapters/summarize_economics.py

No simulation is rerun. The illustrative reference is R_H=15 mm/event.
"""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from soilspice import config, energy
from soilspice import energy_sensitivity as S


def site_statistics(frame):
    cols = ["energy_none_total_loss_pct", "energy_none_band_extra_loss_pct", "cost_none_band_extra_pct_clean_revenue",
            "cost_opt_band_extra_pct_clean_revenue", "margin10_npv25_eur_kwp", "margin10_net_pct_clean_revenue",
            "manual_gain_pct_clean_revenue", "selected_gain_pct_clean_revenue", "opt_total_cleanings_per_year"]
    rows = []
    for (site, market, threshold), q in frame.groupby(["site", "market", "R_H_mm_event"], sort=False):
        r = dict(site=site, market=market, R_H_mm_event=threshold, configurations=len(q),
                 margin10_profitable_configurations=int((q.margin10_state == "positive").sum()),
                 manual_cleaning_profitable_configurations=int((q.manual_gain_state == "positive").sum()))
        for c in cols:
            r.update({c + "_median": q[c].median(), c + "_min": q[c].min(), c + "_max": q[c].max()})
        rows.append(r)
    return pd.DataFrame(rows)


def write_readme(out, ref, robust_ref, rankings_ref, transitions_ref):
    lines = ["# Sensibilité au seuil de lavage par événement", "",
        "Calcul reproductible : `python chapters/summarize_economics.py`. Synthèse des simulations sauvegardées ; aucune simulation ni trajectoire existante n’est modifiée.", "",
        "## Périmètre", "",
        "- Seuils : 5, 10, 15, 20, 25 mm/événement ; référence illustrative R_H=15 mm/événement (choix utilisateur). P_e >= R_H suppose un lavage complet des deux dépôts ; sinon transfert si P_e >= 1 mm. Une mise à jour à la clôture après 3 heures consécutives sèches.",
        "- Référence de la figure et des petits tableaux : bande rectangulaire de 60 mm, f=0,1, marge initiale 0 mm, 4 modules × 2 orientations × 3 sites. Les exports complets conservent toutes les hauteurs et les trois fractions déjà simulées, sans agréger ces paramètres entre seuils.",
        "- Deux marchés séparés : ground (0,06 €/kWh ; 0,41 €/kWp/nettoyage), rooftop (0,28 €/kWh ; 3,3 €/kWp/nettoyage). Les données existantes proposent au sol : aucun nettoyage, déclenchement économique, intervalles 14/30/60/90/120/180/270/365 jours ; en toiture : aucun nettoyage ou déclenchement économique seulement. Le robot est exclu de l’optimisation manuelle, conformément aux figures existantes (son investissement n’est pas inclus dans le coût brut sauvegardé).", "",
        "## Conventions et unités", "",
        "Soit E0 la production propre et Eperdue l’énergie perdue, pour le même site, module, orientation, seuil et paramètres. Les comparaisons entre seuils gardent les mêmes données météo, E0, prix, coûts de nettoyage, dépôt, alpha, h, f et grille de stratégies ; input_audit.json en vérifie les invariants et les empreintes des fichiers.", "",
        "- Énergie sans nettoyage manuel : L_total = Eperdue_avec_bande/E0 ; L_uniforme = Eperdue_sans_bande/E0 ; L_bande_supplémentaire = L_total − L_uniforme. Colonnes energy_none_* en % de la production propre. La bande n’est donc pas ajoutée deux fois à la perte uniforme. Les fractions du bilan 2016–2025 sont pondérées par la production propre, pas des moyennes arithmétiques des pourcentages annuels. annual_energy_reference.csv donne aussi les résultats année par année.",
        "- Coût annuel total C = p Eperdue + c Nnettoyages ; revenu propre V0 = p E0. Les colonnes *_pct_clean_revenue donnent 100 C/V0 ou 100 ΔC/V0. Un écart entre ces pourcentages s’exprime en points de pourcentage de revenu propre.",
        "- Sans nettoyage : coût supplémentaire de bande ΔC_none = C_avec,none − C_sans,none. Après optimisation : ΔC_opt = min_s C_avec,s − min_t C_sans,t ; s et t sont choisis indépendamment. Les politiques, intervalles et fréquences retenus sont exportés séparément. Le coût de référence sans bande reste celui du modèle uniforme aux mêmes seuils et météo. On ne soustrait pas le coût sans bande à la seule stratégie choisie avec bande.",
        "- Gain manuel signé = C_avec,none − min_{s manuel} C_avec,s, coûts de nettoyage déjà inclus. Les politiques manuelles excluent none et robot ; un gain négatif ou nul est conservé. Le gain de la politique sélectionnée parmi none et les stratégies manuelles est exporté séparément (selected_gain_*), ainsi que le gain manuel du cas sans bande et la réduction du surcoût de bande due à l’optimisation. Un gain total de nettoyage ne doit pas être interprété comme un gain dû à la seule bande.",
        "- Fréquence associée : opt_total_cleanings_per_year pour la politique finalement choisie avec bande ; manual_total_cleanings_per_year pour la meilleure politique manuelle, même lorsqu’elle n’est pas rentable. La politique sans bande possède sa propre fréquence. Une politique manuelle peut ne déclencher aucun nettoyage. Les fréquences minimales/maximales des politiques co-optimales sont également exportées ; le minimum exact est retenu, départagé par l’ordre des lignes sauvegardées en cas d’égalité (none, trigger, puis intervalles croissants).", "",
        "### Marge de 10 mm", "",
        "La convention existante est conservée : coût du module 100 €/kWp (0,10 €/Wp), investissement K10 = 100 × 10/L_pente en €/kWp, L_pente en mm. Le module et sa puissance active restent identiques ; seule la hauteur effectivement ombrée passe de h à max(h−10,0), sans changer la masse de bande ni son modèle de croissance. Chaque marge possède sa propre optimisation, avec et sans bande.",
        "Économie annuelle S10 = ΔC_opt,marge0 − ΔC_opt,marge10. Facteur d’actualisation AF = (1−1,05^(−25))/0,05. Bénéfice net sur 25 ans = AF S10 − K10, en €/kWp ; bénéfice annuel équivalent = S10 − K10/AF, en €/kWp/an et en % du revenu propre. La figure utilise le bénéfice net sur 25 ans, investissement déduit, et non la seule valeur brute de l’énergie économisée. La rentabilité correspond au signe strictement positif de ce bénéfice net.", "",
        "### Années et périodes censurées", "",
        "Le nombre moyen d’événements est le compte sur la période 2016–2025 divisé par 10. Les coûts annuels, rendements annuels et fréquences manuelles conservent la normalisation durée couverte/365,25. Après correction des timestamps Open-Meteo, la dernière heure source manque : la couverture finit le 31 décembre 2025 à 23:00 UTC, soit 3652 jours et 23 heures, sans météo synthétique. Les colonnes study_start_utc, study_stop_utc et missing_trailing_source_hours explicitent cette couverture. Les colonnes *_per_calendar_year donnent aussi les fréquences manuelles ramenées à 10 années nominales. Les pertes et coûts en % ne dépendent pas de cette différence de normalisation. La définition économique de la marge reste identique.",
        "Les intervalles sont mesurés entre clôtures d’événements lavants, sans remise à zéro au changement d’année. max_complete_washing_gap_days exclut les bornes censurées. Le début d’étude → premier lavage et le dernier lavage → fin d’étude sont exportés séparément avec dates et durées ; leur durée observée ne constitue pas un intervalle complet entre deux lavages. Tous les temps sont UTC et suivent le catalogue utilisé par les simulations principales.", "",
        "## Robustesse", "",
        "Chaque module/orientation est conservé. Les signes positif/négatif/zéro sont évalués aux cinq seuils (tolérance numérique 10⁻⁹ €/kWp/an ; orientation 10⁻⁸ point de pourcentage). Zéro est une indifférence numérique, pas une preuve de rentabilité. Ces tolérances ne sont pas des seuils de significativité physique ou statistique.",
        "Le classement portrait/paysage est évalué séparément pour la perte supplémentaire de bande sans nettoyage et pour le surcoût de bande après optimisations indépendantes. Un écart portrait − paysage positif favorise le paysage, négatif favorise le portrait. Les passages par une égalité sont conservés comme changements de classement.",
        "transitions_*.csv indique chaque intervalle de seuils adjacents où change un signe ou un classement, avec les deux valeurs et états. La position exacte de la bascule n’est pas interpolée : aucune simulation n’est disponible entre ces seuils. Les lignes de robustesse indiquent aussi tous les états dans l’ordre 5/10/15/20/25, les seuils rentables et si le signe/classement reste constant.", "",
        "## Résultats à h=60 mm, f=0,1", "",
        "| Site | Marché | Configurations marge rentable aux 5 seuils / 8 | Nettoyage avantageux aux 5 seuils / 8 | Classement énergétique stable / 4 modules | Classement économique stable / 4 modules |",
        "|---|---|---|---|---|---|",
    ]
    for (site, market), q in robust_ref.groupby(["site", "market"], sort=False):
        rank = rankings_ref[(rankings_ref.site == site) & (rankings_ref.market == market)]
        lines.append(f"| {site} | {market} | {q.margin10_all_profitable.sum()} | {q.manual_cleaning_all_profitable.sum()} | {rank.energy_ranking_stable.sum()} | {rank.cost_ranking_stable.sum()} |")
    lines += ["", "Les changements individuels de signe/classement sont dans transitions_reference.csv ; ces comptes ne remplacent pas les résultats par configuration.", "",
        "[Conclusions par module/orientation et liste des bascules](conclusions_reference.md) : seuils rentables, séquences de classement et valeurs de chaque changement. Une inversion de très faible amplitude ne démontre pas un avantage matériel ou une précision physique du modèle ; ses deux valeurs sont conservées.", "",
        "## Fichiers", "",
        "- summary_reference.csv : 240 lignes, 3 sites × 5 seuils × 8 configurations × 2 marchés, h=60 mm et f=0,1.",
        "- summary_all_configurations.csv : mêmes métriques pour toutes les hauteurs et fractions sauvegardées, deux marchés ; margin10 compare toujours 0 et +10 mm.",
        "- annual_energy_reference.csv : pertes totales, uniformes et supplémentaires de bande pour chaque année civile, sans nettoyage manuel, h=60 mm/f=0,1.",
        "- washing_summary.csv, washing_intervals.csv, censored_washing_intervals.csv : comptes, intervalles complets et bornes censurées, communs aux configurations d’un site/seuil.",
        "- robustness_reference.csv, robustness_all_configurations.csv : conclusions par module/orientation, avec classements énergétiques et économiques joints.",
        "- orientation_reference.csv, orientation_all_configurations.csv : écarts portrait − paysage et rang à chaque seuil ; orientation_robustness_*.csv : stabilité du classement.",
        "- transitions_reference.csv, transitions_all_configurations.csv : changements de signe ou rang entre seuils adjacents.",
        "- site_reference_statistics.csv : médiane, minimum, maximum et nombre de configurations rentables ; il s’agit d’une dispersion entre configurations, jamais d’intervalles de confiance.",
        "- input_audit.json, conventions.json, column_dictionary.csv : invariants vérifiés, paramètres/formules et unités des colonnes.", "",
        "[Figure de synthèse au sol](../../fig_energy_economics_threshold_sensitivity.pdf) · [Figure en toiture](../../fig_energy_economics_threshold_sensitivity_rooftop.pdf)", "",
        "Chaque courbe représente un module/orientation, sans intervalle de confiance. Les axes verticaux sont propres à chaque site.", "",
        "Les trajectoires temporelles des cinq scénarios sont conservées comme figures complémentaires : [comparaison des trajectoires](../../fig_band_growth_total_comparison.pdf) et [index des figures par seuil](../../total_scenarios/README.md).", "",
        "Le lavage complet reste une hypothèse et les seuils ne sont pas calibrés pour les bandes. Aucun seuil central n’est choisi selon la courbe ou le résultat économique.", ""]
    (out / "README.md").write_text("\n".join(lines))
    detail = ["# Conclusions par configuration : h=60 mm, f=0,1", "",
              "R_H en mm/événement. Rentable signifie bénéfice net strictement positif. Les zéros et pertes sont distingués dans robustness_reference.csv. Les marchés ont des prix et des ensembles de stratégies différents ; voir README.md.", ""]
    for (site, market), q in robust_ref.groupby(["site", "market"], sort=False):
        detail += [f"## {site} — {market}", "", "| Module / orientation | Marge +10 mm rentable pour R_H | Nettoyage manuel avantageux pour R_H |", "|---|---|---|"]
        for r in q.itertuples():
            detail.append(f"| {r.module} / {r.orientation} | {r.margin10_profitable_thresholds_mm or 'aucun'} | {r.manual_cleaning_profitable_thresholds_mm or 'aucun'} |")
        detail += ["", "Classement dans l’ordre R_H = 5 / 10 / 15 / 20 / 25 ; tie = égalité numérique.", "",
                   "| Module | Orientation préférable : pertes de bande | Orientation préférable : coût optimisé de bande |", "|---|---|---|"]
        for r in rankings_ref[(rankings_ref.site == site) & (rankings_ref.market == market)].itertuples():
            detail.append(f"| {r.module} | {r.energy_rankings_5_10_15_20_25.replace(';', ' / ')} | {r.cost_rankings_5_10_15_20_25.replace(';', ' / ')} |")
        detail.append("")
    detail += ["## Changements entre seuils adjacents", "", "La bascule est localisée dans le segment testé, sans interpolation d’un seuil précis. Les signes de marge et nettoyage portent sur des bénéfices annuels équivalents en €/kWp/an ; les classements portent sur l’écart portrait − paysage en points de pourcentage. Les passages par zéro/égalité sont signalés.", "",
               "| Site / marché | Module / orientation | Conclusion | R_H bas → haut | État bas → haut | Valeur basse → haute | Unité |", "|---|---|---|---|---|---|---|"]
    for r in transitions_ref.itertuples():
        detail.append(f"| {r.site} / {r.market} | {r.module} / {r.orientation} | {r.conclusion} | {r.R_H_low_mm:g} → {r.R_H_high_mm:g} | {r.state_low} → {r.state_high} | {r.value_low:.7g} → {r.value_high:.7g} | {r.unit} |")
    tiny = transitions_ref[transitions_ref.conclusion.str.startswith("orientation_") & (transitions_ref[["value_low", "value_high"]].abs().max(axis=1) < .001)]
    if len(tiny):
        detail += ["", "Attention : certaines inversions ci-dessus ont des écarts aux deux bornes inférieurs à 0,001 point de pourcentage. Elles sont conservées dans les classements, mais ne constituent pas une preuve d’un avantage matériel ou d’une précision numérique/physique à cette échelle."]
    (out / "conclusions_reference.md").write_text("\n".join(detail) + "\n")


def main():
    root, out = config.CH4_DIR, config.ECON_DIR / "threshold_sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    summary, audit = S.load_and_audit(root)
    wash, intervals = S.washing_summary(pd.read_csv(config.ECON_DIR / "rain_total_washing_intervals.csv"))
    summary = summary.merge(wash, on=["site", "R_H_mm_event"], validate="many_to_one")
    for col in ("opt_total", "opt_uniform", "manual_total", "manual_uniform", "margin10_opt"):
        summary[f"{col}_cleanings_per_calendar_year"] = summary[f"{col}_cleanings_per_year"] * summary.economic_years_365p25 / 10
    # Independent rain-only reference outputs must give the same band losses.
    ref = S.reference(summary)
    existing = pd.read_csv(config.ECON_DIR / "rain_total_reference_metrics.csv").rename(columns={"P_H_mm": "R_H_mm_event"})
    check = ref.merge(existing, on=["site", "R_H_mm_event", "module", "orientation"], suffixes=("", "_reference"), validate="many_to_one")
    assert len(check) == 240
    np.testing.assert_allclose(check.energy_none_total_loss_pct, check.total_energy_loss_pct, atol=1e-10)
    np.testing.assert_allclose(check.energy_none_band_extra_loss_pct, check.additional_band_energy_loss_pct, atol=1e-10)
    np.testing.assert_array_equal(check.washing_events, check.washing_events_reference)
    orientations = S.orientation_comparison(summary)
    robust, ranking_robust, transitions = S.classify_across_thresholds(summary, orientations)
    frames = {"summary_all_configurations": summary, "summary_reference": ref,
              "washing_summary": wash, "washing_intervals": intervals,
              "censored_washing_intervals": intervals[intervals.left_censored | intervals.right_censored],
              "robustness_all_configurations": robust, "robustness_reference": S.reference(robust),
              "orientation_all_configurations": orientations, "orientation_reference": S.reference(orientations),
              "orientation_robustness_all_configurations": ranking_robust, "orientation_robustness_reference": S.reference(ranking_robust),
              "transitions_all_configurations": transitions, "transitions_reference": S.reference(transitions),
              "site_reference_statistics": site_statistics(ref)}
    annual = pd.read_csv(config.ECON_DIR / "rain_total_annual_metrics.csv").rename(columns={"P_H_mm": "R_H_mm_event"})
    annual["uniform_energy_loss_pct"] = annual.total_energy_loss_pct - annual.additional_band_energy_loss_pct
    frames["annual_energy_reference"] = annual
    for name, frame in frames.items():
        frame.to_csv(out / (name + ".csv"), index=False)
    audit["rain_only_reference_reconciles"] = True
    audit["washing_gaps_span_calendar_years"] = True
    audit["rows"] = {name: len(frame) for name, frame in frames.items()}
    audit["supporting_inputs"] = {name: hashlib.sha256((config.ECON_DIR / name).read_bytes()).hexdigest()
                                 for name in ("rain_total_reference_metrics.csv", "rain_total_washing_intervals.csv", "rain_total_annual_metrics.csv")}
    (out / "input_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    conventions = dict(thresholds_mm_event=list(S.THRESHOLDS), selected_central_threshold=energy.REFERENCE_TOTAL,
                       reference=dict(shape="rectangle", h_mm=60, f=.1, base_margin_mm=0, compared_margin_mm=10),
                       markets=energy.MARKETS, annual_economic_normalization="duration_days / 365.25; identical to saved economic figures",
                       washing_count_normalization="event count / 10 nominal years, 2016–2025; final source hour unavailable", economic_years=float(wash.economic_years_365p25.iloc[0]),
                       hourly_time_convention=audit['assumptions']['hourly_time_convention'],
                       source_coverage=dict(start_utc=str(wash.study_start_utc.min()), stop_utc=str(wash.study_stop_utc.max()),
                                            missing_trailing_hours=float(wash.missing_trailing_source_hours.max())),
                       energy_normalization="clean-production-weighted aggregate fractions, annual values separately exported",
                       independent_optimization=True, include_none_in_optimum=True, include_none_in_manual_only=False,
                       robot_included=False, discount_rate=.05, lifetime_years=25, annuity_factor=energy.AF,
                       module_cost_eur_kwp=S.MODULE_PRICE_EUR_KWP, added_margin_mm=10,
                       margin_capex_formula="100 * 10 / slope_length_mm [EUR/kWp]",
                       margin_net_npv_formula="AF * (extra_band_cost_opt_margin0 - extra_band_cost_opt_margin10) - capex",
                       cost_opt_band_formula="min(cost_with_band) - min(cost_without_band), independently",
                       manual_gain_formula="cost_with_band_none - min(cost_with_band_manual_only); cleaning charges included",
                       orientation_metrics=["extra band energy loss without manual cleaning", "extra band cost after independent optimization"],
                       EUR_numeric_tolerance=S.EPS_EUR, percentage_point_numeric_tolerance=S.EPS_PCT,
                       transitions="adjacent simulated threshold brackets; no inferred exact crossing", washing="assumed complete; uncalibrated for bands")
    conventions["figure_representation"] = "one line per module/orientation; no confidence intervals; statistics file min/max is dispersion between configurations"
    (out / "conventions.json").write_text(json.dumps(conventions, indent=2, ensure_ascii=False) + "\n")
    dictionary = []
    for col in summary:
        if "pct" in col:
            unit = "% clean production" if col.startswith("energy") or "energy_loss" in col else "% clean revenue"
        elif "npv25" in col or col == "margin10_capex_eur_kwp":
            unit = "EUR/kWp (25-year net present value)" if "npv25" in col else "EUR/kWp (initial investment)"
        elif "eur_kwp_year" in col:
            unit = "EUR/kWp/year (365.25-day normalization)"
        elif "kwh_kwp_year" in col:
            unit = "kWh/kWp/year (365.25-day normalization)"
        elif "per_calendar_year" in col:
            unit = "events or cleanings / calendar year"
        elif "frequency" in col or "cleanings_per_year" in col:
            unit = "manual cleanings / year (365.25-day normalization)"
        elif "_days" in col:
            unit = "days"
        elif col == "missing_trailing_source_hours":
            unit = "hours"
        elif col == "R_H_mm_event":
            unit = "mm/event"
        elif col.endswith("_mm"):
            unit = "mm"
        else:
            unit = "dimensionless, count, category or UTC date; see README"
        dictionary.append(dict(column=col, unit=unit))
    pd.DataFrame(dictionary).to_csv(out / "column_dictionary.csv", index=False)
    S.plot_summary(summary, root / "fig_energy_economics_threshold_sensitivity", "ground")
    S.plot_summary(summary, root / "fig_energy_economics_threshold_sensitivity_rooftop", "rooftop")
    write_readme(out, ref, S.reference(robust), S.reference(ranking_robust), S.reference(transitions))
    print("Rows:", audit["rows"], flush=True)
    print(S.reference(robust).groupby(["site", "market"])[["margin10_all_profitable", "manual_cleaning_all_profitable"]].sum().to_string(), flush=True)
    print(S.reference(transitions).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
