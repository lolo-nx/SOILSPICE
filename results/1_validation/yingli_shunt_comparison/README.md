# Yingli shunt-resistance comparison

Reproduce: `MPLBACKEND=Agg python chapters/yingli_shunt_comparison.py`.

Only the cell shunt resistance changes. The template, diode, series resistance, bypass diodes, photocurrent, temperatures and band geometry are identical between the two models. No Bishop breakdown term is introduced. Constant Rsh is 5.99 ohm/cell. The PVsyst law uses R0=4*5.99, exponent 5.5 and g=1-f_cov*(1-tau_I), including partial cell coverage.

All 18 campaign states are plotted: 14 shaded measurements and the same clean flash repeated for four geometry/orientation comparisons. Transmittances are 1, .706, .536, .442 and 0; .442 was not measured in landscape. Both rectangular and triangular 80 mm bands are included. Measured shaded curves use the existing clean-flash temperature correction.

`fig_yingli_shunt_comparison` is the overview; the two `by_transmittance` figures use blue solid curves for PVsyst and grey dashed curves for constant Rsh. `fig_yingli_shunt_current_change` shows the small direct current change in mA. PDF and PNG versions are supplied.

`curve_errors.csv` compares both simulations on identical 0.1 V samples of their shared measured voltage support, excluding measured current below 0.3 A, as in the existing validation. IV and PV errors are normalized by measured Isc and Pmax. The 5-25 V current-window metric is predefined and includes transitions for triangular bands; it is not a fitted plateau region. `error_summary.csv` gives arithmetic means over shaded transmittances, excluding clean curves.

The improvement is modest and does not resolve the portrait plateau discrepancy. IV nRMSE decreases for all 14 shaded curves, but PV error and maximum-power error do not improve uniformly. The clean simulations are identical. No parameters were fitted.

Existing publication exports were preserved. They still describe the previous constant-Rsh simulations; this comparison does not regenerate chapter 2, energy or economic results.

Mean shaded-curve errors:

```text
    shape orientation    model  nRMSE_IV_pct  nRMSE_PV_pct  plateau_5_25V_RMSE_A
rectangle   landscape constant      1.055150      1.107691              0.112111
rectangle   landscape   pvsyst      1.054059      1.108463              0.111839
rectangle    portrait constant      4.796798      1.640220              0.363470
rectangle    portrait   pvsyst      4.646790      1.600101              0.352615
 triangle   landscape constant      2.926808      3.240014              0.248337
 triangle   landscape   pvsyst      2.875457      3.178599              0.244533
 triangle    portrait constant      4.325943      3.551139              0.453873
 triangle    portrait   pvsyst      4.292759      3.515565              0.450942
```
