# Sensibilité au seuil de lavage par événement

Calcul reproductible : `python chapters/summarize_economics.py`. Synthèse des simulations sauvegardées ; aucune simulation ni trajectoire existante n’est modifiée.

## Périmètre

- Seuils : 5, 10, 15, 20, 25 mm/événement ; référence illustrative R_H=15 mm/événement (choix utilisateur). P_e >= R_H suppose un lavage complet des deux dépôts ; sinon transfert si P_e >= 1 mm. Une mise à jour à la clôture après 3 heures consécutives sèches.
- Référence de la figure et des petits tableaux : bande rectangulaire de 60 mm, f=0,1, marge initiale 0 mm, 4 modules × 2 orientations × 3 sites. Les exports complets conservent toutes les hauteurs et les trois fractions déjà simulées, sans agréger ces paramètres entre seuils.
- Deux marchés séparés : ground (0,06 €/kWh ; 0,41 €/kWp/nettoyage), rooftop (0,28 €/kWh ; 3,3 €/kWp/nettoyage). Les données existantes proposent au sol : aucun nettoyage, déclenchement économique, intervalles 14/30/60/90/120/180/270/365 jours ; en toiture : aucun nettoyage ou déclenchement économique seulement. Le robot est exclu de l’optimisation manuelle, conformément aux figures existantes (son investissement n’est pas inclus dans le coût brut sauvegardé).

## Conventions et unités

Soit E0 la production propre et Eperdue l’énergie perdue, pour le même site, module, orientation, seuil et paramètres. Les comparaisons entre seuils gardent les mêmes données météo, E0, prix, coûts de nettoyage, dépôt, alpha, h, f et grille de stratégies ; input_audit.json en vérifie les invariants et les empreintes des fichiers.

- Énergie sans nettoyage manuel : L_total = Eperdue_avec_bande/E0 ; L_uniforme = Eperdue_sans_bande/E0 ; L_bande_supplémentaire = L_total − L_uniforme. Colonnes energy_none_* en % de la production propre. La bande n’est donc pas ajoutée deux fois à la perte uniforme. Les fractions du bilan 2016–2025 sont pondérées par la production propre, pas des moyennes arithmétiques des pourcentages annuels. annual_energy_reference.csv donne aussi les résultats année par année.
- Coût annuel total C = p Eperdue + c Nnettoyages ; revenu propre V0 = p E0. Les colonnes *_pct_clean_revenue donnent 100 C/V0 ou 100 ΔC/V0. Un écart entre ces pourcentages s’exprime en points de pourcentage de revenu propre.
- Sans nettoyage : coût supplémentaire de bande ΔC_none = C_avec,none − C_sans,none. Après optimisation : ΔC_opt = min_s C_avec,s − min_t C_sans,t ; s et t sont choisis indépendamment. Les politiques, intervalles et fréquences retenus sont exportés séparément. Le coût de référence sans bande reste celui du modèle uniforme aux mêmes seuils et météo. On ne soustrait pas le coût sans bande à la seule stratégie choisie avec bande.
- Gain manuel signé = C_avec,none − min_{s manuel} C_avec,s, coûts de nettoyage déjà inclus. Les politiques manuelles excluent none et robot ; un gain négatif ou nul est conservé. Le gain de la politique sélectionnée parmi none et les stratégies manuelles est exporté séparément (selected_gain_*), ainsi que le gain manuel du cas sans bande et la réduction du surcoût de bande due à l’optimisation. Un gain total de nettoyage ne doit pas être interprété comme un gain dû à la seule bande.
- Fréquence associée : opt_total_cleanings_per_year pour la politique finalement choisie avec bande ; manual_total_cleanings_per_year pour la meilleure politique manuelle, même lorsqu’elle n’est pas rentable. La politique sans bande possède sa propre fréquence. Une politique manuelle peut ne déclencher aucun nettoyage. Les fréquences minimales/maximales des politiques co-optimales sont également exportées ; le minimum exact est retenu, départagé par l’ordre des lignes sauvegardées en cas d’égalité (none, trigger, puis intervalles croissants).

### Marge de 10 mm

La convention existante est conservée : coût du module 100 €/kWp (0,10 €/Wp), investissement K10 = 100 × 10/L_pente en €/kWp, L_pente en mm. Le module et sa puissance active restent identiques ; seule la hauteur effectivement ombrée passe de h à max(h−10,0), sans changer la masse de bande ni son modèle de croissance. Chaque marge possède sa propre optimisation, avec et sans bande.
Économie annuelle S10 = ΔC_opt,marge0 − ΔC_opt,marge10. Facteur d’actualisation AF = (1−1,05^(−25))/0,05. Bénéfice net sur 25 ans = AF S10 − K10, en €/kWp ; bénéfice annuel équivalent = S10 − K10/AF, en €/kWp/an et en % du revenu propre. La figure utilise le bénéfice net sur 25 ans, investissement déduit, et non la seule valeur brute de l’énergie économisée. La rentabilité correspond au signe strictement positif de ce bénéfice net.

### Années et périodes censurées

Le nombre moyen d’événements est le compte sur la période 2016–2025 divisé par 10. Les coûts annuels, rendements annuels et fréquences manuelles conservent la normalisation durée couverte/365,25. Après correction des timestamps Open-Meteo, la dernière heure source manque : la couverture finit le 31 décembre 2025 à 23:00 UTC, soit 3652 jours et 23 heures, sans météo synthétique. Les colonnes study_start_utc, study_stop_utc et missing_trailing_source_hours explicitent cette couverture. Les colonnes *_per_calendar_year donnent aussi les fréquences manuelles ramenées à 10 années nominales. Les pertes et coûts en % ne dépendent pas de cette différence de normalisation. La définition économique de la marge reste identique.
Les intervalles sont mesurés entre clôtures d’événements lavants, sans remise à zéro au changement d’année. max_complete_washing_gap_days exclut les bornes censurées. Le début d’étude → premier lavage et le dernier lavage → fin d’étude sont exportés séparément avec dates et durées ; leur durée observée ne constitue pas un intervalle complet entre deux lavages. Tous les temps sont UTC et suivent le catalogue utilisé par les simulations principales.

## Robustesse

Chaque module/orientation est conservé. Les signes positif/négatif/zéro sont évalués aux cinq seuils (tolérance numérique 10⁻⁹ €/kWp/an ; orientation 10⁻⁸ point de pourcentage). Zéro est une indifférence numérique, pas une preuve de rentabilité. Ces tolérances ne sont pas des seuils de significativité physique ou statistique.
Le classement portrait/paysage est évalué séparément pour la perte supplémentaire de bande sans nettoyage et pour le surcoût de bande après optimisations indépendantes. Un écart portrait − paysage positif favorise le paysage, négatif favorise le portrait. Les passages par une égalité sont conservés comme changements de classement.
transitions_*.csv indique chaque intervalle de seuils adjacents où change un signe ou un classement, avec les deux valeurs et états. La position exacte de la bascule n’est pas interpolée : aucune simulation n’est disponible entre ces seuils. Les lignes de robustesse indiquent aussi tous les états dans l’ordre 5/10/15/20/25, les seuils rentables et si le signe/classement reste constant.

## Résultats à h=60 mm, f=0,1

| Site | Marché | Configurations marge rentable aux 5 seuils / 8 | Nettoyage avantageux aux 5 seuils / 8 | Classement énergétique stable / 4 modules | Classement économique stable / 4 modules |
|---|---|---|---|---|---|
| Neuchâtel | ground | 0 | 0 | 3 | 3 |
| Neuchâtel | rooftop | 3 | 0 | 3 | 3 |
| Almería | ground | 7 | 8 | 4 | 4 |
| Almería | rooftop | 8 | 8 | 4 | 4 |
| Doha | ground | 6 | 8 | 3 | 4 |
| Doha | rooftop | 8 | 8 | 3 | 4 |

Les changements individuels de signe/classement sont dans transitions_reference.csv ; ces comptes ne remplacent pas les résultats par configuration.

[Conclusions par module/orientation et liste des bascules](conclusions_reference.md) : seuils rentables, séquences de classement et valeurs de chaque changement. Une inversion de très faible amplitude ne démontre pas un avantage matériel ou une précision physique du modèle ; ses deux valeurs sont conservées.

## Fichiers

- summary_reference.csv : 240 lignes, 3 sites × 5 seuils × 8 configurations × 2 marchés, h=60 mm et f=0,1.
- summary_all_configurations.csv : mêmes métriques pour toutes les hauteurs et fractions sauvegardées, deux marchés ; margin10 compare toujours 0 et +10 mm.
- annual_energy_reference.csv : pertes totales, uniformes et supplémentaires de bande pour chaque année civile, sans nettoyage manuel, h=60 mm/f=0,1.
- washing_summary.csv, washing_intervals.csv, censored_washing_intervals.csv : comptes, intervalles complets et bornes censurées, communs aux configurations d’un site/seuil.
- robustness_reference.csv, robustness_all_configurations.csv : conclusions par module/orientation, avec classements énergétiques et économiques joints.
- orientation_reference.csv, orientation_all_configurations.csv : écarts portrait − paysage et rang à chaque seuil ; orientation_robustness_*.csv : stabilité du classement.
- transitions_reference.csv, transitions_all_configurations.csv : changements de signe ou rang entre seuils adjacents.
- site_reference_statistics.csv : médiane, minimum, maximum et nombre de configurations rentables ; il s’agit d’une dispersion entre configurations, jamais d’intervalles de confiance.
- input_audit.json, conventions.json, column_dictionary.csv : invariants vérifiés, paramètres/formules et unités des colonnes.

[Figure de synthèse au sol](../../fig_energy_economics_threshold_sensitivity.pdf) · [Figure en toiture](../../fig_energy_economics_threshold_sensitivity_rooftop.pdf)

Chaque courbe représente un module/orientation, sans intervalle de confiance. Les axes verticaux sont propres à chaque site.

Les trajectoires temporelles des cinq scénarios sont conservées comme figures complémentaires : [comparaison des trajectoires](../../fig_band_growth_total_comparison.pdf) et [index des figures par seuil](../../total_scenarios/README.md).

Le lavage complet reste une hypothèse et les seuils ne sont pas calibrés pour les bandes. Aucun seuil central n’est choisi selon la courbe ou le résultat économique.
