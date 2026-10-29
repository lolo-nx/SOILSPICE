# Résultats énergie et économie

Le modèle retenu prend une décision par événement, clos après trois heures consécutives sans précipitation.
Lavage complet supposé si le cumul atteint `R_H` ; sinon transfert à partir de 1 mm/événement.
La référence illustrative demandée est **R_H = 15 mm/événement** ; les seuils ne sont pas calibrés pour les bandes.

| Usage | Dossier |
|---|---|
| Figures principales, tableaux et CSV de référence | [total_scenarios/PH15mm](total_scenarios/PH15mm/) |
| Sensibilité à 5 mm/événement | [total_scenarios/PH5mm](total_scenarios/PH5mm/) |
| Sensibilité à 10 mm/événement | [total_scenarios/PH10mm](total_scenarios/PH10mm/) |
| Sensibilité à 20 mm/événement | [total_scenarios/PH20mm](total_scenarios/PH20mm/) |
| Sensibilité à 25 mm/événement | [total_scenarios/PH25mm](total_scenarios/PH25mm/) |

[Index des figures des cinq scénarios](total_scenarios/README.md).

Les comparaisons et figures de sensibilité sont conservées à la racine : `fig_rain_total_comparison`,
`fig_band_growth_total_comparison`, `fig_energy_economics_threshold_sensitivity*` et
`fig_washing_threshold_sensitivity_*` (PDF et PNG).
Leurs données et conventions sont dans `data/rain_total_*`, [data/threshold_sensitivity](data/threshold_sensitivity/README.md)
et [data/washing_threshold_sensitivity_3x3](data/washing_threshold_sensitivity_3x3/README.md).

```sh
python chapters/4_economics.py                    # référence 15 mm/événement
python chapters/4_economics.py --washing-total 10 # un autre seuil
python chapters/4_economics.py --all-thresholds   # les cinq seuils et leur comparaison
python chapters/summarize_economics.py            # synthèse depuis les CSV enregistrés
python chapters/fig_washing_threshold_sensitivity.py # figures 3 × 3 depuis les CSV
python chapters/rain_event_audit.py               # contrôle des événements et actions uniques
```

Le nettoyage du dépôt a supprimé les anciens scénarios en mm/h, les essais de comparaison MeteoSwiss/ERA5,
leurs diagnostics et leurs scripts/tests spécifiques. Les CSV et figures des cinq scénarios retenus,
les résultats des chapitres 1–3 et les données météorologiques sources sont conservés.
Les tests du modèle actuel restent dans `tests/`.

**Regénération terminée avec les timestamps corrigés.** Les cinq seuils, leurs figures PDF/PNG, tableaux,
CSV et synthèses de sensibilité ont été recalculés ensemble. La référence reste **R_H=15 mm/événement** ;
les paramètres physiques, économiques et les stratégies de nettoyage sont inchangés.

[Rapport de vérification](data/corrected_timestamp_regeneration.json) : **51 figures PDF/PNG**, **336 960 lignes de simulations**,
**8 640 lignes de synthèse économique**, **240 lignes de référence** et **30 tests réussis**. Les contrôles relient
les valeurs tracées aux simulations, vérifient les comptes d'événements et l'investissement robot ajouté une seule fois.

Les timestamps Open-Meteo de fin d'heure sont convertis en débuts d'intervalles après le calcul PV. La couverture
source se termine le 31 décembre 2025 à 23:00 UTC ; aucune météo synthétique n'a été ajoutée pour la dernière heure
manquante. Les conventions et périodes censurées sont explicites dans les tableaux de sensibilité.

[Revue historique du calage horaire](data/timing_review/README.md) : les colonnes avant/après documentent la correction ;
la référence corrigée de cette revue correspond maintenant à la référence publiée.

[Rapport historique du nettoyage et des vérifications](data/cleanup_report.json) : 297 fichiers obsolètes supprimés,
452 CSV/PDF/PNG conservés à l’identique, 22 tests réussis ; trajectoires reproduites aux cinq seuils
et 8 640 lignes de synthèse économique vérifiées contre les simulations enregistrées.
