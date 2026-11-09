# Variantes des couleurs de régime

Les couleurs sont exactement celles de la palette fournie. L’ordre est I / II / III.
I et III utilisent des teintes proches ; II utilise une autre couleur. Les figures principales sont conservées.

[Planche comparative](COMPARAISON_TRIOS.pdf) · [Six pages à taille complète](TOUS_LES_TRIOS.pdf)

| Variante | I : current limited | II : bypass clamped | III : mixed | Figure |
|---|---|---|---|---|
| A | cyan `#CCEEFF` | jaune `#EEEEBB` | bleu `#AACCEE` | [PDF](A_cyan_jaune_bleu.pdf) |
| B | cyan `#CCEEFF` | rose `#FFBBCC` | bleu `#AACCEE` | [PDF](B_cyan_rose_bleu.pdf) |
| C | cyan `#CCEEFF` | vert `#BBDDBB` | bleu `#AACCEE` | [PDF](C_cyan_vert_bleu.pdf) |
| D | rose `#FFBBCC` | cyan `#CCEEFF` | violet `#EEBBDD` | [PDF](D_rose_cyan_violet.pdf) |
| E | rose `#FFBBCC` | vert `#BBDDBB` | violet `#EEBBDD` | [PDF](E_rose_vert_violet.pdf) |
| F | rose `#FFBBCC` | jaune `#EEEEBB` | violet `#EEBBDD` | [PDF](F_rose_jaune_violet.pdf) |

Chaque variante conserve les 800 états, les 16 cartes, les limites d’axes et la mise en page.
Les PDF sont vectoriels et les PNG individuels sont exportés à 600 dpi.

Reproduction depuis SOILSPICE :

```sh
MPLCONFIGDIR=/tmp/soilspice-matplotlib MPLBACKEND=Agg python3 results/2_simulation/regime_palette_variants/build.py --related
```
