# Graphiques : cumul par événement

**Référence illustrative : R_H = 15 mm/événement** (non calibrée pour les bandes). Les quatre autres seuils sont des scénarios de sensibilité.

Cinq scénarios : R_H = 5, 10, 15, 20 et 25 mm/événement. Les événements se terminent après 3 heures consécutives sans pluie. Le lavage complet est supposé lorsque le cumul atteint R_H ; sinon, transfert si le cumul atteint 1 mm. Les paramètres de dépôt et de géométrie sont identiques entre scénarios.

Chaque figure est exportée en PDF et PNG. Les CSV sont dans le sous-dossier data de chaque scénario ; les tableaux se trouvent à côté des figures.

| Figure | 5 mm/événement | 10 mm/événement | 15 mm/événement | 20 mm/événement | 25 mm/événement |
|---|---|---|---|---|---|
| Croissance de la bande | [PDF](PH5mm/fig_band_growth.pdf) · [PNG](PH5mm/fig_band_growth.png) | [PDF](PH10mm/fig_band_growth.pdf) · [PNG](PH10mm/fig_band_growth.png) | [PDF](PH15mm/fig_band_growth.pdf) · [PNG](PH15mm/fig_band_growth.png) | [PDF](PH20mm/fig_band_growth.pdf) · [PNG](PH20mm/fig_band_growth.png) | [PDF](PH25mm/fig_band_growth.pdf) · [PNG](PH25mm/fig_band_growth.png) |
| Cartes des pertes énergétiques | [PDF](PH5mm/fig_energy_maps.pdf) · [PNG](PH5mm/fig_energy_maps.png) | [PDF](PH10mm/fig_energy_maps.pdf) · [PNG](PH10mm/fig_energy_maps.png) | [PDF](PH15mm/fig_energy_maps.pdf) · [PNG](PH15mm/fig_energy_maps.png) | [PDF](PH20mm/fig_energy_maps.pdf) · [PNG](PH20mm/fig_energy_maps.png) | [PDF](PH25mm/fig_energy_maps.pdf) · [PNG](PH25mm/fig_energy_maps.png) |
| Leviers de conception | [PDF](PH5mm/fig_design_levers.pdf) · [PNG](PH5mm/fig_design_levers.png) | [PDF](PH10mm/fig_design_levers.pdf) · [PNG](PH10mm/fig_design_levers.png) | [PDF](PH15mm/fig_design_levers.pdf) · [PNG](PH15mm/fig_design_levers.png) | [PDF](PH20mm/fig_design_levers.pdf) · [PNG](PH20mm/fig_design_levers.png) | [PDF](PH25mm/fig_design_levers.pdf) · [PNG](PH25mm/fig_design_levers.png) |
| Hauteur et géométrie | [PDF](PH5mm/fig_design_height.pdf) · [PNG](PH5mm/fig_design_height.png) | [PDF](PH10mm/fig_design_height.pdf) · [PNG](PH10mm/fig_design_height.png) | [PDF](PH15mm/fig_design_height.pdf) · [PNG](PH15mm/fig_design_height.png) | [PDF](PH20mm/fig_design_height.pdf) · [PNG](PH20mm/fig_design_height.png) | [PDF](PH25mm/fig_design_height.pdf) · [PNG](PH25mm/fig_design_height.png) |
| Sensibilité | [PDF](PH5mm/fig_sensitivity.pdf) · [PNG](PH5mm/fig_sensitivity.png) | [PDF](PH10mm/fig_sensitivity.pdf) · [PNG](PH10mm/fig_sensitivity.png) | [PDF](PH15mm/fig_sensitivity.pdf) · [PNG](PH15mm/fig_sensitivity.png) | [PDF](PH20mm/fig_sensitivity.pdf) · [PNG](PH20mm/fig_sensitivity.png) | [PDF](PH25mm/fig_sensitivity.pdf) · [PNG](PH25mm/fig_sensitivity.png) |
| Nettoyage et coûts | [PDF](PH5mm/fig_cleaning.pdf) · [PNG](PH5mm/fig_cleaning.png) | [PDF](PH10mm/fig_cleaning.pdf) · [PNG](PH10mm/fig_cleaning.png) | [PDF](PH15mm/fig_cleaning.pdf) · [PNG](PH15mm/fig_cleaning.png) | [PDF](PH20mm/fig_cleaning.pdf) · [PNG](PH20mm/fig_cleaning.png) | [PDF](PH25mm/fig_cleaning.pdf) · [PNG](PH25mm/fig_cleaning.png) |
| Forme de bande | [PDF](PH5mm/fig_band_shape.pdf) · [PNG](PH5mm/fig_band_shape.png) | [PDF](PH10mm/fig_band_shape.pdf) · [PNG](PH10mm/fig_band_shape.png) | [PDF](PH15mm/fig_band_shape.pdf) · [PNG](PH15mm/fig_band_shape.png) | [PDF](PH20mm/fig_band_shape.pdf) · [PNG](PH20mm/fig_band_shape.png) | [PDF](PH25mm/fig_band_shape.pdf) · [PNG](PH25mm/fig_band_shape.png) |
| Diagnostic inverse | [PDF](PH5mm/fig_inverse_chain.pdf) · [PNG](PH5mm/fig_inverse_chain.png) | [PDF](PH10mm/fig_inverse_chain.pdf) · [PNG](PH10mm/fig_inverse_chain.png) | [PDF](PH15mm/fig_inverse_chain.pdf) · [PNG](PH15mm/fig_inverse_chain.png) | [PDF](PH20mm/fig_inverse_chain.pdf) · [PNG](PH20mm/fig_inverse_chain.png) | [PDF](PH25mm/fig_inverse_chain.pdf) · [PNG](PH25mm/fig_inverse_chain.png) |
| Cartes par cycle (SI) | [PDF](PH5mm/fig_pan_h_maps_SI.pdf) · [PNG](PH5mm/fig_pan_h_maps_SI.png) | [PDF](PH10mm/fig_pan_h_maps_SI.pdf) · [PNG](PH10mm/fig_pan_h_maps_SI.png) | [PDF](PH15mm/fig_pan_h_maps_SI.pdf) · [PNG](PH15mm/fig_pan_h_maps_SI.png) | [PDF](PH20mm/fig_pan_h_maps_SI.pdf) · [PNG](PH20mm/fig_pan_h_maps_SI.png) | [PDF](PH25mm/fig_pan_h_maps_SI.pdf) · [PNG](PH25mm/fig_pan_h_maps_SI.png) |

[Comparaison des indicateurs](../fig_rain_total_comparison.pdf) · [Comparaison de la croissance](../fig_band_growth_total_comparison.pdf)

Ces seuils ne sont pas calibrés pour les bandes. Aucun seuil central unique n’est sélectionné.

Les cinq dossiers ont été régénérés avec les timestamps Open-Meteo corrigés en débuts d’intervalles UTC.
[Source et contrôles de cohérence](../data/corrected_timestamp_regeneration.json).
