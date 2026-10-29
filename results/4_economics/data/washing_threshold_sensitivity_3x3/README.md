# Sensibilité au seuil de lavage — nouvelles figures 3 × 3

Reproduction : `python chapters/fig_washing_threshold_sensitivity.py`. Les figures principales et les anciennes figures de sensibilité sont conservées.

Les PDF vectoriels font exactement 183 × 170 mm ; les PNG sont exportés à 600 dpi. Une couleur par module et un style par orientation ; panneaux a–i, sans barre d’incertitude. Les hypothèses détaillées restent dans les légendes ci-dessous, hors de l’image.

Données : plotted_data_ground.csv et plotted_data_rooftop.csv (120 lignes chacun), ou plotted_data.csv (240 lignes). Les trois indicateurs tracés et leurs composantes sont vérifiés par rapport à summary_reference.csv et recalculés depuis les cinq model_runs.csv. Les gains manuels négatifs sont conservés ; la ligne nettoyage concerne la salissure totale, pas uniquement la bande.

Le scénario R_H=15 mm/événement est une référence illustrative demandée pour les figures principales. Aucun seuil n’est déterminé ici à partir de l’apparence des courbes ou d’un optimum économique.

## Conclusions à paramètres fixés

Nombre de configurations, sur 8, dont le bénéfice reste strictement positif aux cinq seuils :

| Site | Marché | Marge +10 mm | Nettoyage manuel |
|---|---|---|---|
| Neuchâtel | ground | 0/8 | 0/8 |
| Neuchâtel | rooftop | 3/8 | 0/8 |
| Almería | ground | 7/8 | 8/8 |
| Almería | rooftop | 8/8 | 8/8 |
| Doha | ground | 6/8 | 8/8 |
| Doha | rooftop | 8/8 | 8/8 |

0/8 signifie qu’aucune configuration n’est rentable aux cinq seuils, pas qu’elle ne l’est jamais. L’énergie supplémentaire de bande reste positive pour les 8 configurations à chaque site et seuil. Pour les signes de marge et nettoyage, zero signifie absence de gain net, non un gain positif.

### Intervalles de changement de signe

Les intervalles ci-dessous encadrent uniquement des changements entre deux simulations adjacentes. Aucun seuil exact de bascule n’est interpolé. Les valeurs détaillées et les unités sont dans sign_changes.csv ; les signes aux cinq seuils pour chaque configuration sont dans robustness_by_configuration.csv.

| Site / marché | Module / orientation | Indicateur | Intervalle R_H [mm/événement] | Signe avant → après |
|---|---|---|---|---|
| Neuchâtel / ground | YINGLI / portrait | margin10 | 15–20 | negative → positive |
| Neuchâtel / ground | YINGLI / landscape | margin10 | 20–25 | negative → positive |
| Neuchâtel / ground | SUNPOWER / portrait | margin10 | 15–20 | negative → positive |
| Neuchâtel / ground | SUNPOWER / portrait | manual_cleaning | 20–25 | zero → positive |
| Neuchâtel / ground | SUNPOWER / landscape | margin10 | 15–20 | negative → positive |
| Neuchâtel / ground | JASOLAR / portrait | margin10 | 10–15 | negative → positive |
| Neuchâtel / ground | JASOLAR / portrait | manual_cleaning | 20–25 | zero → positive |
| Neuchâtel / ground | JASOLAR / landscape | margin10 | 20–25 | negative → positive |
| Neuchâtel / ground | LONGI / portrait | margin10 | 10–15 | negative → positive |
| Neuchâtel / ground | LONGI / portrait | manual_cleaning | 20–25 | zero → positive |
| Neuchâtel / ground | LONGI / landscape | margin10 | 15–20 | negative → positive |
| Neuchâtel / rooftop | YINGLI / portrait | margin10 | 5–10 | negative → positive |
| Neuchâtel / rooftop | YINGLI / landscape | margin10 | 10–15 | negative → positive |
| Neuchâtel / rooftop | SUNPOWER / portrait | manual_cleaning | 20–25 | zero → negative |
| Neuchâtel / rooftop | SUNPOWER / landscape | margin10 | 5–10 | negative → positive |
| Neuchâtel / rooftop | SUNPOWER / landscape | manual_cleaning | 20–25 | zero → negative |
| Neuchâtel / rooftop | JASOLAR / portrait | manual_cleaning | 20–25 | zero → negative |
| Neuchâtel / rooftop | JASOLAR / landscape | margin10 | 5–10 | negative → positive |
| Neuchâtel / rooftop | JASOLAR / landscape | manual_cleaning | 20–25 | zero → negative |
| Neuchâtel / rooftop | LONGI / portrait | manual_cleaning | 20–25 | zero → negative |
| Neuchâtel / rooftop | LONGI / landscape | margin10 | 5–10 | negative → positive |
| Neuchâtel / rooftop | LONGI / landscape | manual_cleaning | 20–25 | zero → negative |
| Almería / ground | YINGLI / landscape | margin10 | 5–10 | negative → positive |
| Doha / ground | YINGLI / landscape | margin10 | 10–15 | negative → positive |
| Doha / ground | YINGLI / landscape | margin10 | 15–20 | positive → negative |
| Doha / ground | JASOLAR / landscape | margin10 | 5–10 | negative → positive |

## Figure caption — ground

Sensitivity to the assumed washing threshold, ground market. Rectangular band, h = 60 mm, f = 0.1, initial added margin 0 mm; all other assumptions fixed. Columns: Neuchâtel, Almería and Doha. a–c, additional annual energy loss due to the band without manual cleaning: loss with band minus loss without band, as a percentage of identical clean production. d–f, net present value of an additional 10 mm margin over 25 years at a 5% discount rate, in EUR/kWp, including the extra-glass investment; costs with and without band are independently optimised at each margin, as in the existing economic analysis. g–i, annual net gain from manual cleaning for total soiling (uniform deposit plus band), as a percentage of clean revenue: total cost with no manual cleaning minus the minimum total cost of the evaluated manual policies, including cleaning charges. This selection considers economic-trigger and fixed-interval manual policies (14, 30, 60, 90, 120, 180, 270 and 365 days), excluding both no manual cleaning and the robot; negative gains are retained. Colours identify modules; solid lines/circles denote portrait and dashed lines/triangles landscape. The dotted vertical line marks R_H = 15 mm/event, an illustrative reference specified for the main analysis, not a threshold chosen by optimising these curves. Horizontal lines mark zero economic benefit. Panels a, g use a separate Neuchâtel scale; b–c and h–i share scales, and d–f use a common scale. Curves show individual configurations; no uncertainty bars or confidence intervals are used. Complete removal of uniform and band deposits is assumed when event precipitation P_e >= R_H. Events are separated by three consecutive dry hours. These washing thresholds are not calibrated for bottom-edge bands. The annual normalisation and available policy sets follow the saved simulations.

## Figure caption — rooftop

Sensitivity to the assumed washing threshold, rooftop market. Rectangular band, h = 60 mm, f = 0.1, initial added margin 0 mm; all other assumptions fixed. Columns: Neuchâtel, Almería and Doha. a–c, additional annual energy loss due to the band without manual cleaning: loss with band minus loss without band, as a percentage of identical clean production. d–f, net present value of an additional 10 mm margin over 25 years at a 5% discount rate, in EUR/kWp, including the extra-glass investment; costs with and without band are independently optimised at each margin, as in the existing economic analysis. g–i, annual net gain from manual cleaning for total soiling (uniform deposit plus band), as a percentage of clean revenue: total cost with no manual cleaning minus the minimum total cost of the evaluated manual policies, including cleaning charges. This selection considers the economic-trigger manual policy, the only manual policy evaluated for this market, excluding both no manual cleaning and the robot; negative gains are retained. Colours identify modules; solid lines/circles denote portrait and dashed lines/triangles landscape. The dotted vertical line marks R_H = 15 mm/event, an illustrative reference specified for the main analysis, not a threshold chosen by optimising these curves. Horizontal lines mark zero economic benefit. Neuchâtel uses a separate vertical scale in each row; Almería and Doha share a scale within each row. Curves show individual configurations; no uncertainty bars or confidence intervals are used. Complete removal of uniform and band deposits is assumed when event precipitation P_e >= R_H. Events are separated by three consecutive dry hours. These washing thresholds are not calibrated for bottom-edge bands. The annual normalisation and available policy sets follow the saved simulations.

[Trajectoires temporelles complémentaires](../../fig_band_growth_total_comparison.pdf) · [Conventions complètes des tableaux de synthèse](../threshold_sensitivity/README.md)
