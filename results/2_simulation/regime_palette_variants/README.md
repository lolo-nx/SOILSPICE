# Variantes des couleurs de régime

Les couleurs sont exactement celles de la palette fournie. L’ordre est I / II / III.
A reprend la proposition actuelle ; les figures principales ne sont pas modifiées.

[Planche comparative](COMPARAISON_TRIOS.pdf) · [Six pages à taille complète](TOUS_LES_TRIOS.pdf)

| Variante | I : current limited | II : bypass clamped | III : mixed | Figure |
|---|---|---|---|---|
| A | cyan `#CCEEFF` | jaune `#EEEEBB` | violet `#EEBBDD` | [PDF](A_cyan_jaune_violet.pdf) |
| B | cyan `#CCEEFF` | vert `#BBDDBB` | rose `#FFBBCC` | [PDF](B_cyan_vert_rose.pdf) |
| C | gris `#DDDDDD` | cyan `#CCEEFF` | rose `#FFBBCC` | [PDF](C_gris_cyan_rose.pdf) |
| D | cyan `#CCEEFF` | jaune `#EEEEBB` | vert `#BBDDBB` | [PDF](D_cyan_jaune_vert.pdf) |
| E | bleu `#AACCEE` | jaune `#EEEEBB` | rose `#FFBBCC` | [PDF](E_bleu_jaune_rose.pdf) |
| F | gris `#DDDDDD` | vert `#BBDDBB` | violet `#EEBBDD` | [PDF](F_gris_vert_violet.pdf) |

Chaque variante conserve les 800 états, les 16 cartes, les limites d’axes et la mise en page.
Les PDF sont vectoriels et les PNG individuels sont exportés à 600 dpi.

Reproduction depuis SOILSPICE :

```sh
MPLCONFIGDIR=/tmp/soilspice-matplotlib MPLBACKEND=Agg python3 results/2_simulation/regime_palette_variants/build.py
```
