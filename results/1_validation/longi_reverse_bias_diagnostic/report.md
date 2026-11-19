# LONGi portrait reverse-bias diagnostic

Existing 80 mm rectangular band; original parameters, 25 °C, original 0.2 V DC grid. Global MPP is the maximum over the entire first-quadrant sweep, with 0.2 V resolution. Voc is interpolated at the current zero crossing separately for each shaded state. The displayed validation value 0.44 is stored as tau_I=0.442 in experiments.py; this diagnostic preserves that exact campaign value.

## Geometry and signs

Shaded cells: XX1, XX18, XX19, XX36, XX37, XX54; f_cov = 0.5846 (unrounded geometric coverage 0.5845511482). All are in the bottom row, columns 1–6. The effective relative photocurrent is 1 − f_cov(1 − tau_I).

The three series substrings each contain two parallel 18-cell arms: XX1–18 / XX55–72, XX19–36 / XX73–90, XX37–54 / XX91–108. Two shaded cells occur in the lower arm of each substring; the parallel arm is unshaded.

V_terminal = V(+) − V(−); Vj = V(internal N) − V(−). Positive module output is I(V1). R2 is oriented internal N → positive terminal, so positive cell/arm output is +I(R2), not its negative. I(B1), I(R1), I(D1) are signed from junction to negative terminal. KCL: I_arm = Iph − I(R1) − I(D1) − I(B1). The reverse Bishop contribution to output is −I(B1). For each substring I_module = I_shaded_arm + I_unshaded_arm + I_bypass. Diode bypass current is positive in its forward conduction direction.

The unchanged Bishop source is I(B1) = (Vj/5.5438) × 0.002 × abs(1 − Vj/(−4.5))^(−3); LTspice pwr uses the absolute base. The ordinary shunt contribution −I(R1) is separate from Bishop breakdown.

## Short circuit and global MPP

Values below are the mean over six shaded cells; full individual values are in summary_per_cell.csv.

| tau_I | Point | Vmodule (V) | Imodule (A) | Vcell (V) | Vj (V) | Iarm (A) | Iparallel (A) | −IBishop (A) | Ishunt signed (A) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | SC | 0.00000 | 14.30285 | -4.17002 | -4.15618 | 7.10325 | 7.19960 | 3.36178 | -0.74970 |
| 0 | MPP | 35.20000 | 9.71950 | 0.07779 | 0.08359 | 2.97666 | 6.74284 | -0.00003 | 0.01508 |
| 0.442 | SC | 0.00000 | 14.30597 | -4.06971 | -4.05587 | 7.10637 | 7.19960 | 1.52201 | -0.73160 |
| 0.442 | MPP | 35.40000 | 11.49526 | 0.26167 | 0.27103 | 4.80378 | 6.69148 | -0.00008 | 0.04889 |

## Mechanism assessment

- tau_I=0: proposed SC-to-MPP mechanism SUPPORTED. All six junctions reverse biased at SC: True; all still reverse biased at MPP: False. Maximum shaded-cell spread in Vj over the full curve: 3.81e-06 V. SC-to-MPP Bishop nonincreasing (tolerance 1e-5 A): True. Maximum absolute bypass current over that interval: 4.24e-09 A.
- tau_I=0.442: proposed SC-to-MPP mechanism SUPPORTED. All six junctions reverse biased at SC: True; all still reverse biased at MPP: False. Maximum shaded-cell spread in Vj over the full curve: 3.81e-06 V. SC-to-MPP Bishop nonincreasing (tolerance 1e-5 A): True. Maximum absolute bypass current over that interval: 4.25e-09 A.

The decline must be interpreted through branch redistribution: module current sums the two parallel arm currents. A declining Bishop contribution is direct evidence for the proposed mechanism within this model, but does not establish it as the unique physical cause in the measurements. Ordinary shunt current also changes; no refitting or mechanism ablation was performed.

For tau_I=0, between V/Voc=0.4 and 0.8, module current decreases by 0.5999 A. The shaded arm decreases by 0.5093 A and its parallel arm by 0.0905 A. The Bishop contribution decreases by only 0.0341 A, whereas the reverse ordinary-shunt contribution decreases by 0.4753 A. Thus Bishop breakdown explains the early drop, but the gradual central decline is mainly the ordinary shunt contribution decreasing as reverse voltage relaxes.

For tau_I=0.442, between V/Voc=0.4 and 0.8, module current decreases by 0.5849 A. The shaded arm decreases by 0.4935 A and its parallel arm by 0.0913 A. The Bishop contribution decreases by only 0.0223 A, whereas the reverse ordinary-shunt contribution decreases by 0.4713 A. Thus Bishop breakdown explains the early drop, but the gradual central decline is mainly the ordinary shunt contribution decreasing as reverse voltage relaxes.

## Validation and reproducibility

Only copied netlists are run. Instrumentation adds a .save directive; the macOS library path is resolved to the installed Windows standard.dio. The bare source placeholder V1 N097 0 V is written as V1 N097 0 0 for Windows syntax; the original .dc sweep overrides this initial source value. No model values or original files are changed. The local diode library and source files are hashed in audit.json. Local versus previously saved validation curves may differ with LTspice/library versions; the discrepancies below quantify this. Instrumentation produces zero change in module current; the saved validation curves agree to within 1 microampere. With the original solver tolerances, the maximum full-sweep cell KCL residual is 4.17 mA and substring KCL residual is 2.62 mA; branch series currents agree exactly at RAW precision. These residuals are retained and reported rather than changing solver settings.

```json
{
  "checks": {
    "kcl_A": 0.004165902370004915,
    "series_ohm_V": 2.7088556289675547e-06,
    "shunt_ohm_A": 5.707818324096614e-07,
    "bishop_formula_A": 1.7222650423587282e-05,
    "parallel_kcl_A": 0.0026116327078362644,
    "series_branch_A": 0.0
  },
  "simulations": [
    {
      "tau_I": 0.0,
      "instrumented_IV_max_change_A": 0.0,
      "saved_validation_IV_max_change_A": 9.536743341698184e-07,
      "Voc_V": 40.08163859094664,
      "Vmp_V": 35.2,
      "sweep_step_V": 0.2
    },
    {
      "tau_I": 0.442,
      "instrumented_IV_max_change_A": 0.0,
      "saved_validation_IV_max_change_A": 1.8918200339612667e-13,
      "Voc_V": 40.11831461875795,
      "Vmp_V": 35.4,
      "sweep_step_V": 0.2
    }
  ]
}
```

Files: cell_traces.csv (all sweep points, including points beyond Voc); summary_per_cell.csv; summary.csv; fig_longi_reverse_bias_diagnostic.png/.pdf; audit.json; runs/ (baseline and instrumented netlists/raw).