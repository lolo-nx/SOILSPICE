# PH15 mm: robotic cleaning economic revision

## What changed

Previous ownership charge: investment / 10, without discounting, replacement timing or maintenance. New charge: discounted 25-year ownership cost converted to an equivalent annual cost, plus unchanged saved annual lost revenue.

Main figure: **linear remaining-life residual credit**, as selected by the user. The named alternative assumes zero residual value.

For initial investment C, r=0.05 and A=sum((1+r)^(-t), t=1..25):

    PV_no_residual = C [1 + (1+r)^(-12) + (1+r)^(-24)] + 0.05 C A
    PV_with_residual = PV_no_residual - (11/12) C (1+r)^(-25)
    EAC = PV_with_residual / A  (main case)
    annual total cost (% clean revenue) = 100 (saved energy-loss cost + EAC) / clean annual revenue

For manual strategies, multiplying annual costs and revenue by the same annuity factor leaves the plotted cost/revenue ratio unchanged. Comparing their existing annual costs to robot EAC is therefore consistent. The 2016–2025 meteorological record still supplies mean annual quantities; it is not a separate 10-year economic horizon.

## Exact plotted robot values

| Site | Residual assumption | Band low (%) | Band midpoint (%) | Band high (%) | No-band midpoint (%) |
|---|---|---:|---:|---:|---:|
| Neuchâtel | none | 2.446654899 | 5.347332081 | 8.248009264 | 5.347332081 |
| Almería | none | 1.881370782 | 3.979358067 | 6.077345351 | 3.979188788 |
| Doha | none | 1.758345689 | 3.696765339 | 5.635184988 | 3.696726954 |
| Neuchâtel | linear_remaining_life | 2.193965073 | 4.789308716 | 7.384652359 | 4.789308716 |
| Almería | linear_remaining_life | 1.698606549 | 3.575753718 | 5.452900887 | 3.575584439 |
| Doha | linear_remaining_life | 1.589482044 | 3.323858122 | 5.058234199 | 3.323819737 |

## Impact on other panels and conclusions

Only robot costs, their investment-range error bars and no-band robot markers in panel g change. Panels a–f do not use robot ownership costs: the manual optima explicitly exclude robots. The cleaning frequencies, physical energy losses, transmittances, band growth, threshold sensitivities and electrical calculations are unchanged. Other rainfall scenarios are not re-exported; a future full run uses the corrected robot ownership accounting.

A common ownership charge cancels exactly in the paired band-minus-no-band robot difference for each configuration. Site-level differences of medians need not equal the median of paired differences. All uncertainties here are investment bounds, separate from module/orientation variability.

- Almería: old robot total 2.246508857%; new robot total 3.575753718%; best displayed manual strategy Trigger rule: 5.434432695%; rain only 10.218918860%. Robot cheaper than the best displayed manual strategy at midpoint: True. At low and high robot investment respectively: True, False. Robot band-only effect (median paired difference) 0.000126936 percentage points.
- Doha: old robot total 2.095712281%; new robot total 3.323858122%; best displayed manual strategy Trigger rule: 8.547368741%; rain only 32.704909882%. Robot cheaper than the best displayed manual strategy at midpoint: True. At low and high robot investment respectively: True, True. Robot band-only effect (median paired difference) 0.000032372 percentage points.
- Neuchâtel: old robot total 2.951494795%; new robot total 4.789308716%; best displayed manual strategy Best fixed interval: 0.475287294%; rain only 0.475287294%. Robot cheaper than the best displayed manual strategy at midpoint: False. At low and high robot investment respectively: False, False. Robot band-only effect (median paired difference) 0.000000000 percentage points.

The robot almost removes the band-specific loss, but its total economic cost includes a substantial ownership charge. Avoid interpreting robot total soiling cost as a band-only loss or claiming one strategy is optimal across the whole investment range from midpoint bars alone.

See manuscript_text.md and the .tex snippets for ready-to-use Results, Methods and caption text; CSVs retain full precision. All physical inputs and saved model runs are protected by hashes; numerical cash-flow tests are independent of the figure implementation.
