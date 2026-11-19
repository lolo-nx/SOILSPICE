# Annual irradiance-only sensitivity, PH15 mm

## Scope

Reference: rectangular band, h=60 mm, f=0.1, no added glass margin, rain only, all four modules in portrait and landscape, 2016–2025 weather. Original manuscript/economic outputs are not overwritten.

**No irradiance cutoff is used in the economic energy calculation.** The field-test power-ratio analysis uses 400 W/m²; field-test energy integration uses 20 W/m². Neither filter is applied to the economic model.

Hourly clean-production weights come from the unchanged energy.site_hours/pv_yield pipeline, retaining transposition, IAM and Faiman temperature weighting. Effective irradiance and cell temperature are recovered algebraically from the two existing clean-yield series with distinct temperature coefficients, and all four module yields are checked against the recovered inputs. No weather is invented.

The saved JA Solar landscape, 60 mm trajectory is reused exactly. For another orientation/module, tau is its power L_s/L_s,ref, as follows from the existing rain-only mass balance. Event closures, source hour boundaries and midnight intervals are unchanged. Uniform dust attenuation is reconstructed from the same saved event totals and deposition rate, and checked against the saved energy-weighted interval attenuation. All 24 original band-loss values must reproduce within 0.00002 pp before LTspice annual calculations proceed.

## Calculation

For each real hour i, Y_i is clean energy and U_i is uniform-layer transmission. The original additional band loss is 100 sum(Y_i U_i L_archive(tau_i))/sum(Y_i). The primary grid irradiance is G_i U_i, because uniform dust attenuates the irradiance reaching both covered and uncovered cells. An alternative using weather G_i before uniform attenuation is also exported.

At fixed electrical temperature 25 °C, paired clean/soiled LTspice calculations give ΔL(G,tau)=L_25(G,tau)−L_25(1000,tau). We apply L_corrected=L_archive(tau)+ΔL(G,tau), rather than replacing the archived baseline. This explicitly separates the irradiance effect from the previously identified 27 °C versus 25 °C nominal-reference discrepancy and STC LUT interpolation differences. The original tau quantisation (0.005) is retained. Linear interpolation is used in irradiance and tau; no statistical confidence interval is implied.

Electrical grid: G=[1, 5, 20, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200] W/m²; tau=[0.0, 0.1, 0.2, 0.3, 0.4, 0.44, 0.5, 0.6, 0.7, 0.71, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0]. All main-grid clean and soiled powers are paired at the same G and 25 °C. Previously computed identical netlists are reused where possible.

No hours or energy are silently dropped. Outside the solved G range, the original STC contribution is retained and a conservative 0–100% band-loss bound is exported. Values outside physical loss bounds are clipped only after baseline anchoring; counts and energy effect of any clipping are reported.

## Annual results

| Site | Original median (%) | Irradiance-corrected median (%) | Difference of medians (pp) | Paired differences, min–max (pp) |
|---|---:|---:|---:|---:|
| Almería | 3.768371 | 3.794740 | +0.026369 | -0.110049 to +0.123305 |
| Doha | 6.979457 | 6.939551 | -0.039907 | -0.189298 to +0.277446 |
| Neuchâtel | 0.136420 | 0.140246 | +0.003826 | +0.001512 to +0.005493 |

Minimum supported clean-energy fraction: 99.999466880%. Maximum unsupported-domain bound width: 0.000531 pp. Maximum clipping contribution: 0.000000000 pp.

The difference of site medians is not necessarily the median of paired configuration differences; both are exported. Annual differences by year and the clean-energy fractions below 20, 200 and 400 W/m² can be derived directly from the CSVs.

## Contribution of low irradiance

| Site | Clean energy below 20 W/m² (%) | Below 200 W/m² (%) | Below 400 W/m² (%) |
|---|---:|---:|---:|
| Neuchâtel | 0.308094 | 11.175442 | 31.131314 |
| Almería | 0.144287 | 5.096328 | 17.361415 |
| Doha | 0.101299 | 4.288609 | 13.225167 |

## Limits

This is an irradiance-only sensitivity of the annual energy estimator, not a validation or a complete replacement economic calculation. Cell temperature affects the existing clean-energy weights, but the relative band-loss correction is evaluated at fixed 25 °C because the original electrical temperature interface remains insufficient. Fixed Rsh and unvalidated low-irradiance parameters remain assumptions.

Cleaning is not reoptimised, and manual/robot costs, rain events and dust trajectories are not modified. These results must not be applied as a single percentage correction to every economic strategy or height.

## Reproduction

    python chapters/annual_irradiance_sensitivity.py --workers 4

See annual_results_per_configuration.csv, annual_site_summary.csv, annual_results_by_year.csv, irradiance_energy_distribution.csv, irradiance_tau_grid.csv, baseline_reproduction.csv and verification.json.

## Interpolation resolution diagnostic

Coarsening the existing grid to G=[1,20,100,400,800,1000,1200] and tau=[0,0.2,0.4,0.6,0.8,0.95,1] changes the primary annual result by at most 0.248322 pp across the 24 configurations. This is a resolution-sensitivity diagnostic, not an independent upper bound on interpolation error. Main results use the full grid.

Adding G=150,300,500,700,900 W/m² to the initial 11-level grid changes annual band loss by at most 0.052214 pp across the 24 configurations. Final main results use this refined 16-level grid. The comparison is retained in refinement_comparison.csv.

## Numerical solver handling

The Normal solver is limited to 20 seconds per new run. On timeout, its log is retained and the identical netlist is retried with the native LTspice Alternate solver (no parameter or tolerance change). Alternate was used in 0 completed run attempts. The paired SunPower clean/shaded solver comparison is saved in alternate_solver_comparison.csv. Completed original-grid runs are reused; interrupted calculations are resumed from per-netlist hash-checked metadata.
