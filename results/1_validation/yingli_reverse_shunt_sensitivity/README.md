# Yingli reverse-only shunt sensitivity

Diagnostic study, not a calibrated physical model. Reproduce with `MPLBACKEND=Agg python chapters/yingli_reverse_shunt_sensitivity.py`; add `--reuse-simulations` to reuse matching saved netlists/raw files.

## Model and scope

Only the shunt element in copies of the original constant-5.99-ohm Yingli netlists changes. Its current is Vj/5.99 for Vj >= 0 and Vj/Rreverse for Vj < 0. Current is continuous at zero; its derivative may have a kink. The law applies to every cell, so forward-biased cells retain the original conductance. Photocurrent, geometry, series resistance, diode parameters, bypass diodes and temperature handling remain unchanged. No Bishop term is added. The production module definition and all pre-existing results are preserved.

Rectangular 80 mm bands are compared at every available portrait/landscape campaign transmittance, including the clean control. The 0.442 landscape case was not measured.

Exact zero reverse conductance makes reverse-biased identical-cell voltage sharing numerically ill-conditioned. Negligible conduction is represented by 1 MΩ, explicitly checked against 10 MΩ. Their maximum module-current difference in the 5–25 V diagnostic window is 3.815 µA. The 1 MΩ case retains about 5 µA or less of reverse shunt current in sampled cell states; it is an approximation to the open-shunt limit, not a literal infinite-resistance simulation.

Numerical convergence uses reltol=1e-5, abstol=1e-10 and vntol=1e-8, with initial node guesses from the original short-circuit solution. These affect numerical solution only, not the circuit equations. The 10 MΩ convergence check uses the original solver defaults.

## Metric definitions

Both orientations use the existing clean-flash temperature corrections. All eight models of a configuration are evaluated on one common 0.1 V grid over their shared measured voltage support. Measured currents below 0.3 A are excluded. I–V nRMSE is normalized by measured Isc. The 5–25 V RMSE is in amperes and includes bypass transitions in landscape, so it is a fixed current window rather than an isolated fitted plateau. Each configuration has equal weight in the shaded-curve mean; clean curves are excluded from that mean. Minor differences from earlier pairwise metrics can result from using the shared support of all eight cases.

Maximum-power error means 100*(Pmax_sim/Pmax_measured−1), signed; positive means overprediction. Pmax is sampled on the original 0.05 V LTspice sweep without refitting.

## Portrait: I–V nRMSE [%]

| Reverse Rsh [Ω] | tau=1 | tau=0.706 | tau=0.536 | tau=0.442 | tau=0 |
|---|---:|---:|---:|---:|---:|
| 5.99 | 0.496670 | 3.647856 | 4.047921 | 4.160754 | 7.301676 |
| 10 | 0.496670 | 2.087682 | 2.129364 | 2.103146 | 4.627988 |
| 20 | 0.496670 | 0.845436 | 0.683974 | 0.582567 | 2.774211 |
| 50 | 0.496670 | 0.426071 | 0.472708 | 0.605471 | 1.951444 |
| 100 | 0.496670 | 0.553792 | 0.726382 | 0.898184 | 1.804357 |
| 200 | 0.496670 | 0.654371 | 0.868170 | 1.053196 | 1.766548 |
| 500 | 0.496670 | 0.721268 | 0.955647 | 1.147731 | 1.756443 |
| negligible | 0.496670 | 0.767832 | 1.014674 | 1.211213 | 1.755087 |

## Portrait: current RMSE in 5–25 V [A]

| Reverse Rsh [Ω] | tau=1 | tau=0.706 | tau=0.536 | tau=0.442 | tau=0 |
|---|---:|---:|---:|---:|---:|
| 5.99 | 0.015164 | 0.325907 | 0.327569 | 0.325112 | 0.470853 |
| 10 | 0.015164 | 0.186124 | 0.169286 | 0.161168 | 0.291029 |
| 20 | 0.015164 | 0.073318 | 0.047368 | 0.038466 | 0.157884 |
| 50 | 0.015164 | 0.012920 | 0.032668 | 0.048546 | 0.085726 |
| 100 | 0.015164 | 0.026065 | 0.056983 | 0.073182 | 0.067952 |
| 200 | 0.015164 | 0.037054 | 0.069501 | 0.085863 | 0.062013 |
| 500 | 0.015164 | 0.043943 | 0.077070 | 0.093538 | 0.059780 |
| negligible | 0.015164 | 0.048603 | 0.082131 | 0.098675 | 0.058924 |

## Landscape: I–V nRMSE [%]

| Reverse Rsh [Ω] | tau=1 | tau=0.706 | tau=0.536 | tau=0 |
|---|---:|---:|---:|---:|
| 5.99 | 0.496670 | 0.850635 | 0.731430 | 1.583385 |
| 10 | 0.496670 | 0.809166 | 0.707994 | 1.551785 |
| 20 | 0.496670 | 0.781815 | 0.697089 | 1.531702 |
| 50 | 0.496670 | 0.767131 | 0.693447 | 1.521173 |
| 100 | 0.496670 | 0.762550 | 0.692731 | 1.517924 |
| 200 | 0.496670 | 0.760321 | 0.692468 | 1.516349 |
| 500 | 0.496670 | 0.759004 | 0.692341 | 1.515421 |
| negligible | 0.496670 | 0.758134 | 0.692269 | 1.514808 |

## Landscape: current RMSE in 5–25 V [A]

| Reverse Rsh [Ω] | tau=1 | tau=0.706 | tau=0.536 | tau=0 |
|---|---:|---:|---:|---:|
| 5.99 | 0.015164 | 0.088845 | 0.084126 | 0.163362 |
| 10 | 0.015164 | 0.082886 | 0.081424 | 0.161363 |
| 20 | 0.015164 | 0.078930 | 0.080265 | 0.160284 |
| 50 | 0.015164 | 0.076803 | 0.079945 | 0.159812 |
| 100 | 0.015164 | 0.076140 | 0.079903 | 0.159684 |
| 200 | 0.015164 | 0.075818 | 0.079894 | 0.159625 |
| 500 | 0.015164 | 0.075628 | 0.079892 | 0.159592 |
| negligible | 0.015164 | 0.075502 | 0.079893 | 0.159571 |

## Maximum power: invariant across all eight resistances

The maximum-power point lies on the unchanged forward-bias branch in these cases. All eight sampled Pmax values are identical for each configuration; a lower plateau error therefore does not imply improved maximum-power agreement.

| Orientation | tau | Measured Pmax [W] | Simulated Pmax [W] | Error [W] | Error [%] |
|---|---:|---:|---:|---:|---:|
| portrait | 1 | 264.855000 | 263.553074 | -1.301926 | -0.491562 |
| portrait | 0.706 | 257.442076 | 256.008962 | -1.433114 | -0.556674 |
| portrait | 0.536 | 248.111792 | 246.052155 | -2.059636 | -0.830124 |
| portrait | 0.442 | 241.658028 | 239.561864 | -2.096164 | -0.867409 |
| portrait | 0 | 201.237777 | 204.073391 | +2.835614 | +1.409086 |
| landscape | 1 | 264.855000 | 263.553074 | -1.301926 | -0.491562 |
| landscape | 0.706 | 250.782015 | 250.692387 | -0.089628 | -0.035739 |
| landscape | 0.536 | 238.478101 | 237.082991 | -1.395110 | -0.585005 |
| landscape | 0 | 182.050943 | 183.892315 | +1.841372 | +1.011460 |

## Findings

Increasing reverse Rsh removes excess portrait plateau current, but agreement is not monotonic across all resistances. The translucent cases eventually undershoot the measured current when reverse leakage is suppressed too far. The opaque case continues to improve toward negligible conduction, while retaining a residual mismatch.

Lowest errors among the tested values are descriptive sensitivity results, not calibration or recommended physical parameters:

| Orientation | tau | Lowest I–V nRMSE at [Ω] | Lowest 5–25 V RMSE at [Ω] |
|---|---:|---:|---:|
| portrait | 0.706 | 50 | 50 |
| portrait | 0.536 | 50 | 50 |
| portrait | 0.442 | 20 | 20 |
| portrait | 0 | negligible | negligible |
| landscape | 0.706 | negligible | negligible |
| landscape | 0.536 | negligible | 500 |
| landscape | 0 | negligible | negligible |

Convergence of simulated current toward negligible leakage is quantified separately from improvement against measurements:

| Reverse Rsh [Ω] | Maximum distance to negligible case, portrait [mA] | Landscape [mA] |
|---:|---:|---:|
| 100 | 43.118477 | 5.902397 |
| 200 | 21.571636 | 2.952133 |
| 500 | 8.629799 | 1.180845 |

Portrait mean nRMSE is lowest at 50 Ω among the tested cases. Landscape mean nRMSE changes only slightly; none of its shaded cases has a higher whole-curve nRMSE than the original baseline. The clean curve is unchanged. These observations support reverse leakage as an important source of the portrait discrepancy, without establishing a unique or physically calibrated reverse shunt resistance.

## Artifacts and verification

- `fig_reverse_shunt_portrait.pdf/png`: full curves and 5–25 V detail for all portrait transmittances.
- `fig_reverse_shunt_landscape.pdf/png`: corresponding landscape checks.
- `fig_reverse_shunt_error_sensitivity.pdf/png`: nRMSE, current-window RMSE and signed Pmax error.
- `curve_metrics.csv`: full-precision errors for all 72 resistance/configuration combinations.
- `iv_curves.csv`: measured and simulated curve samples.
- `mean_shaded_curve_metrics.csv`: arithmetic means across shaded transmittances.
- `sampled_cell_states.csv`: direct internal branch-current checks at module 0, 5, 15, 25 and 35 V.
- `convergence_to_limit.csv` and `negligible_conduction_verification.csv`: endpoint checks.
- `verification.json`: unchanged-input hashes and numerical checks.

Maximum difference from the original 5.99 Ω module-current arrays: 0 A. Maximum KCL residual: 0.075279 mA. Maximum piecewise-shunt current-law residual: 0.609716 µA. Preserved pre-existing scientific/library files: 724.

The inherited `V1 N001 0 V` line generates the same missing-unswept-value warning as the original model; the explicit DC sweep is retained and the original baseline is reproduced. `runs/exploratory/` contains preliminary convergence probes excluded from all reported results.
