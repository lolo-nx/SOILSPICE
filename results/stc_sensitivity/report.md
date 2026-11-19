# Irradiance sensitivity of the STC loss approximation

## Scope and temperature limitation

280 requested loss conditions: 56 band configurations × five irradiances. Four modules; 80 mm rectangular/triangular bands in both orientations with exact tau_I=0.71, 0.44, 0; plus 60 mm rectangular bands at exact tau_I=0.44 in both orientations. These requested transmittances are not silently replaced by the campaign values 0.706/0.442.

The original pipeline writes t into a parameter that controls both photocurrent and diode TNOM. It does not set .temp. LTspice logs confirm operating temperature 27 °C in the archived-reproduction track, with module-specific nominal parameters 24.8 °C (Yingli), 25.4 °C (SunPower), 24.2 °C (JA Solar), 25 °C (LONGi). Therefore the original reference is not strictly a 1000 W/m², actual-25-°C simulation. This cannot be represented as an exact STC reproduction.

First, archived curves are reproduced using the original settings/grid. Separately, the irradiance sensitivity uses explicit .temp 25 AND t=25 for all G, without changing fitted IS/IS1, n, Rs, Rp, breakdown or bypass coefficients. This creates a consistently isothermal 25 °C reference, but its absolute powers can differ from the archive; preflight_explicit_25C_vs_archive.csv quantifies the difference. No archived result is overwritten. The requested L_STC is recomputed within this explicit-25-C track, never borrowed from another track.

The multi-temperature sweep is **not performed**: changing t would change the nominal diode calibration while changing photocurrent, instead of independently moving operating temperature with saturation-current and thermal-voltage responses referenced to a fixed TNOM. The native diode has a temperature law, but the current wrapper does not separate these controls. No arbitrary thermal law, datasheet power correction, or fitted temperature coefficient is added. Constants kv, q and k are declared but do not implement a separate cell-voltage temperature model.

## Model implementation

- Existing soilspice.spice.read_netlist, patch_cells, set_temperature, set_param, run_ltspice and iv_parameters are reused. A native binary RAW reader from the previous LONGi diagnostic supplies the saved traces when PyLTSpice is unavailable.
- Cell photocurrent is [Isc + Isc kc (t−25)] (Vi/1000) [1−f_cov(1−tau_I)]. Vi is set to each requested irradiance; geometric coverage uses the original geometry.fcov_grid/cell_fcov and the same four-decimal netlist precision.
- SunPower and LONGi contain the explicit Bishop branch I_B=(Vj/Rp) a pwr(1−Vj/Vb,−b). Yingli and JA Solar do not contain this branch in their current templates: reverse current is predominantly the ordinary shunt and any native diode leakage. All use their existing fixed shunt resistance.
- Module current sign is the configured +I(V1), checked against terminal voltage; MPP is the global maximum of V I over the entire first quadrant, including all local peaks. The sensitivity sweep is extended 15% and sampled at 0.02 V; original grids are retained for archived reproduction. Validity requires positive power, interior MPP and a Voc crossing.
- Bypass activity is measured directly from the three module bypass diode currents and terminal voltages. An MPP diode is active when its forward current exceeds max(0.1 mA, 1% of module MPP current) and its voltage is forward. The mask identifies which physical bypass diode is active; parallel branches are already represented in the original full-module netlist. These are MPP states, not claims about the entire I–V curve.

## Validation

preflight_archived_reference.csv checks all four clean modules and 48 campaign shading states (two shapes, two orientations, three exact available campaign transmittances). The unavailable 0.442 landscape campaign states use the saved 0.536 state for reproduction only; the main sensitivity always uses the requested exact 0.44. The original SunPower sweep can stop short of Voc, although its MPP is interior. That truncated sweep is preserved only for archived reproduction; every new sensitivity curve must bracket Voc. Archived powers are checked with a 0.02 W reproduction tolerance. The transparent-band control checks every unique geometry at tau_I=1 with loss tolerance 1e-5 percentage points. Every run and any failure is retained in run_manifest.csv and failed_simulations.csv. The worst overall, worst annual-reference and bypass-change cases are recomputed with a 0.01 V grid; grid_convergence.csv reports the resulting deviation changes.

Maximum archived power reproduction error: 1.33e-11 W. True-25-C powers differ from archived powers by 0.228865 to 1.117114%, including operating/nominal temperature settings and numerical grid refinement. This discrepancy is a reported finding, not a silently passed true-STC equivalence. Refining the grid from 0.02 V to 0.01 V for the global worst case, the annual-geometry worst case and the bypass-transition case changes ΔL by at most 0.000120 pp; see grid_convergence.csv.

Result rows: 280 successful / 280 requested; 0 failed loss conditions.

## Numerical findings

| Module | Minimum ΔL (pp) | Maximum ΔL (pp) | Maximum absolute ΔL (pp) |
|---|---:|---:|---:|
| JA Solar (TOPCon) | -3.280553 | 1.308373 | 3.280553 |
| LONGi (HIBC) | -5.234945 | 1.961131 | 5.234945 |
| SunPower (IBC) | -10.444600 | 0.771051 | 10.444600 |
| Yingli (n-Si) | -4.707685 | 2.457674 | 4.707685 |

Largest absolute deviation: **10.444600 pp**, SunPower (IBC), triangle, portrait, h=80 mm, tau_I=0, G=200 W/m². L_STC=24.622996%, L(G)=14.178396%.

Among 224 non-STC conditions, 186 are within 1 pp and 208 are within 2 pp. 16 exceed 2 pp.

1 conditions have a different bypass-diode MPP mask from their own STC state; see bypass_transitions.csv.

Shape maximum absolute deviations: rectangle: 7.385528 pp, triangle: 10.444600 pp.
Orientation maximum absolute deviations: landscape: 6.257792 pp, portrait: 10.444600 pp.
Transmittance maximum absolute deviations: 0.0: 10.444600 pp, 0.44: 6.923295 pp, 0.71: 1.937467 pp.
Band height maximum absolute deviations: 60.0: 0.780923 pp, 80.0: 10.444600 pp.

The 60 mm annual-reference geometries have a worst deviation of 0.780923 pp across the tested irradiances.

## Scientific interpretation and limitations

Changes with irradiance reflect the existing nonlinear circuit: photocurrent decreases with G, whereas fixed shunt conductance, series resistance, diode voltage and reverse current do not scale proportionally. MPP branches and bypass states can therefore change. Differences between modules cannot be attributed solely to architecture: fitted cell parameters and breakdown models also differ.
Geometry and transmittance can change which cells are limiting and the preferred MPP. Use configuration_summary.csv rather than extrapolating one geometry to all cases.
The one MPP bypass transition is SunPower, triangular opaque band, 80 mm, landscape: bypass diode D1 conducts at the STC MPP but not at the 200 W/m² MPP. The largest global exception occurs instead in portrait without an MPP bypass-mask change; reverse/shunt response and competing MPPs also matter.
This is a numerical sensitivity test, **not independent experimental validation**. The discrepancy between the archived nominal reference and true 25 °C, the insufficient temperature interface, constant Rsh and unvalidated low-irradiance/reverse-bias laws limit inference. Small errors for some configurations do not validate the annual STC approximation; large exceptions are explicitly retained.
No annual energy-error estimate is made. Raw weather data exist, but a defensible hourly estimate would require an independently controlled temperature model and the time-varying band transmittance, including irradiance below 200 W/m² and above 1000 W/m². The fixed-tau grid here must not be substituted for that integration.

The annual model also interpolates an archived transmittance table, whereas this test directly simulates tau_I=0.44. The 60 mm results assess that representative geometry, not the full error of the hourly annual pipeline.

## Reproduction

    python chapters/stc_sensitivity.py --workers 4
Netlists, RAW outputs, logs, and per-run metadata are retained in runs/. Successful runs with identical netlist hashes are reused. Windows source-placeholder and .lib path adaptations are confined to copies. The source/input/manuscript hash check is in verification.json.
