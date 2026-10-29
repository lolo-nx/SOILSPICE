# Yingli portrait reverse-bias diagnostic

Reproduce with `MPLBACKEND=Agg python chapters/yingli_reverse_bias_diagnostic.py`.
The script copies the original constant-Rsh netlists from the shunt-comparison
directory and adds only `.save` requests. No electrical parameter, diode model,
geometry or transmittance is changed. No Bishop term is introduced. The module
I-V arrays are exactly identical to the saved baseline in all five cases.

The band is rectangular, 80 mm high, in portrait. The six bottom-row cells
(XX10, XX20, XX30, XX40, XX50, XX60) have coverage 0.3269 as written into LTspice.
The clean case probes the same positions as a control. The six junction voltages
agree within 3.82 microvolts; tables give their means. The results describe the
model's predicted internal state, not measured cell voltages.

`fig_yingli_reverse_bias_diagnostic.pdf/png` shows measured and constant-model
I-V curves, shaded-cell junction voltage, reverse shunt contribution, and current
residuals. `fig_yingli_reverse_bias_by_transmittance.pdf/png` puts the I-V and cell
voltage curves together for every transmittance. The grey region is 5-25 V.
Junction voltage is measured at the diode/shunt node relative to the cell minus
terminal; external terminal voltage also includes the series-resistance drop.

## Cell junction voltage [V]

| tau_I | Module 0 V | 5 V | 10 V | 15 V | 20 V | 25 V |
|---:|---:|---:|---:|---:|---:|---:|
| 1 (clean) | 0.045 | 0.128 | 0.211 | 0.294 | 0.378 | 0.461 |
| 0.706 | -4.166 | -3.524 | -2.813 | -2.068 | -1.305 | -0.529 |
| 0.536 | -4.614 | -3.843 | -3.062 | -2.274 | -1.481 | -0.684 |
| 0.442 | -4.729 | -3.941 | -3.148 | -2.351 | -1.551 | -0.749 |
| 0 | -5.032 | -4.222 | -3.410 | -2.598 | -1.785 | -0.971 |

All six shaded cells are reverse-biased throughout the 5-25 V window. They
return to forward bias at approximately 28.381, 29.278, 29.656 and 30.959 V of
module voltage, in descending shaded-transmittance order.

## Branch currents at module voltage 15 V [A]

| tau_I | Photocurrent | Reverse shunt contribution (-Ish) | Module current | Measured current |
|---:|---:|---:|---:|---:|
| 0.706 | 8.225657 | 0.345293 | 8.570951 | 8.280117 |
| 0.536 | 7.719928 | 0.379585 | 8.099513 | 7.793763 |
| 0.442 | 7.440289 | 0.392521 | 7.832810 | 7.531214 |
| 0 | 6.125393 | 0.433667 | 6.559060 | 6.131493 |

The cell diode current is about -0.961 nA at this point for the shaded cases.
The bypass diodes carry only approximately 3.08 nA in magnitude. Consequently
Icell_out = Iph - Id - Ish is approximately Iph + |Ish|. Branch currents come
directly from LTspice, not from a fitted reverse-bias law. Across all output
samples, KCL closes within 7.25 microamperes; junction/shunt and terminal/series
voltage identities close within 3.66 microvolts (raw traces use float32 storage).

For the opaque case at module short circuit, photocurrent is 6.125 A and reverse
shunt current contributes approximately 0.840 A to the simulated 6.965 A module
current. At module 15 V, the measured current is nearly equal to modeled
photocurrent, while the modeled reverse shunt adds 0.434 A. This is strong
evidence about the mechanism of excess current inside the present model. It
does not uniquely identify the true cell reverse characteristic.

## Breakdown context and interpretation

The largest modeled reverse bias is only -5.032 V. This is substantially less
negative than a typical crystalline-silicon junction-breakdown scale of about
-15 V discussed by [Pannebakker et al. (2017)](https://doi.org/10.1002/pip.2898).
For industrial acid-etched multicrystalline cells, [Breitenstein et al. (2011)](https://doi.org/10.1063/1.3562200)
report early breakdown in the -4 to -9 V range and avalanche breakdown beyond
about -13 V. The short-circuit predictions overlap the low-voltage end of the
early-breakdown range for some of those cells. These are material-specific
examples, not a measured threshold for this n-type Yingli module. Real-cell
mismatch could also concentrate reverse bias more strongly than this identical-
cell model predicts. Reverse breakdown cannot therefore be excluded in the
real module solely from these modeled voltages.

At fixed photocurrent, temperature and other parameters, adding a reverse
breakdown conduction path would tend to increase current supported by a shaded
cell at a given reverse voltage. It is not, by itself, a mechanism for reducing
the already excessive simulated plateau. These diagnostics do not justify
adding Bishop terms as the next correction.

The most directly motivated hypothesis is that the 5.99 ohm cell shunt,
inferred for the original electrical model, is too conductive when extrapolated
into the shaded reverse-bias regime. Its constant linear conductance is not
validated by agreement with the clean forward I-V curve. No replacement shunt
law or new parameter value is fitted here.

Other hypotheses remain possible: errors in effective coverage/photocurrent,
spatially nonuniform collection or leakage within a partly covered cell,
cell-to-cell mismatch, and flash-test transient or temperature effects. The
opaque-band discrepancy shows that filter transmittance calibration alone
cannot explain every case, but geometry/irradiance uncertainty remains. None
of these alternatives is identified uniquely by the current data. A direct
dark and illuminated reverse I-V measurement of the relevant cells, accompanied
by a check of physical band coverage and photocurrent, would discriminate them
more directly than adding an uncalibrated breakdown term.

## Numerical artifacts

- `diagnostic_summary.csv`: compact full-precision per-transmittance results.
- `diagnostic_points_mean.csv`: mean cell state at module 0, 0.1, 1, 5, 10, 15,
  20, 25, 30 and 35 V. Measured current at exact module zero is left empty when
  outside the measured support; it is not extrapolated.
- `diagnostic_points_by_cell.csv`: the same points for each of the six cells.
- `cell_states.csv`: directly saved cell voltages and currents across the full
  first-quadrant sweep, including terminal absorbed power.
- `verification.json`: source hashes, unchanged-module-I-V check, KCL and Ohm's
  law checks, and confirmation that all pre-existing results and library files
  were preserved.
- `runs/`: instrumented netlists, raw outputs and logs. The inherited source
  line `V1 N001 0 V` produces an LTspice warning that its unswept DC value is
  assumed zero; the explicit `.dc V1 0 39.5 0.05` sweep completes and exactly
  reproduces the original I-V arrays. The source line is deliberately retained.
