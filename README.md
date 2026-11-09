# soilspice — bottom-edge soiling bands on PV modules

Code of the paper *Impact of bottom-edge soiling on photovoltaic module performance: experimental validation and
cell-level modeling across module architectures*. One script per chapter of the paper; every paper figure is written,
in the paper style, into the results folder of its chapter.

```
python chapters/run_all.py              # all paper figures/tables, reference R_H=15 mm/event
python chapters/1_validation.py --full  # a chapter from scratch (LTspice for chapters 1 and 2)
```

## Layout

```
soilspice/                 the library (all the functions)
    config.py              paths: the only machine-specific file (flash data, LTspice executable)
    modules.py             the five modules: geometry, netlist, cell numbering, bypass diodes
    geometry.py            band -> covered fraction f_cov of every cell, A_cov, band height for a target A_cov
    spice.py               netlist patching, LTspice batch runs, I-V parameters
    measurements.py        flash-test .asc files, temperature correction
    experiments.py         the flash-test campaign (which file is which band and filter)
    validation.py          SPICE against the flash tests
    simulation.py          loss maps at equal A_cov and at equal band height
    analysis.py            losses, error metrics, electrical regimes
    resilience.py          resilience metrics of the loss maps (S, M, M0, S*, P_an,S)
    uncertainty.py         model uncertainty from the validation residuals
    bypass.py              ideal-bypass model (the middle of the three model levels)
    field.py               roof field test
    energy.py              climate (ERA5, pvlib), band growth model, cleaning strategies, costs
    rain.py                hourly rain events and the event/midnight simulation grid
    rain_scenarios.py       event-total thresholds (reference 15 mm/event), metrics and comparisons
    style.py               paper style shared by all figures
    plots.py               diagnostic figures
    fig_validation.py      figures of chapter 1
    fig_simulation.py      figures of chapter 2
    fig_energy.py          figures of chapter 4
    schematics.py          electrical layout of the modules (Methods)
    tables.py              LaTeX tables
chapters/
    1_validation.py        experimental validation of the model
    2_simulation.py        electrical resilience across modules
    3_field_test.py        roof field test
    4_economics.py         energy yield and cleaning costs
    run_all.py             the four chapters in order
data/
    netlists/              LTspice netlists (.net)
    climate/               ERA5 (Open-Meteo) at the three sites and Kaifeng, MeteoSwiss station NEU
    field_test/            IV-tracer export of the roof test; photos/ = the two photos of fig_field_photos
    inverse_example/       SPICE I-V curves used as a synthetic measurement (inverse diagnosis)
results/
    1_validation/          paper figures and tables of chapter 1; data/ = everything behind them
    2_simulation/          id. chapter 2 (data/ holds the loss maps used by chapters 3 and 4)
    3_field_test/          id. chapter 3
    4_economics/           reference total_scenarios/PH15mm/; sensitivity PH5mm, PH10mm, PH20mm, PH25mm
```

## Paper figures

| Paper | File |
|---|---|
| iv_grid_rectangle, errors_grid, fig_model_hierarchy, pv_module_schematics, Table uncertainty | results/1_validation/ |
| optical_loss_modules_grid, regime_map, regime_atlas, average_power_windows_horizontal, initial_amplification | results/2_simulation/ |
| Table field (summary.json), fig_field_test (SI), fig_field_photos | results/3_field_test/ |
| fig_band_growth, fig_energy_maps, fig_design_height, fig_design_levers, fig_cleaning, fig_inverse_chain, ed_massband_sensitivity (= fig_sensitivity), fig_pan_h_maps (SI), Tables sites and clean | results/4_economics/total_scenarios/PH15mm/ (reference); PH5mm, PH10mm, PH20mm, PH25mm for sensitivity |
| fig_rain_total_comparison, fig_band_growth_total_comparison, energy/economics threshold sensitivity | results/4_economics/ |

Made by hand: the annotated Soiling-images source montage, bishop-circuit,
simulation_framework and spectral_current_final. The soiling montage is reused
without retouching by `soilspice/fig_soiling.py` in the combined figure.

`python Nature-builder-panel.py` assembles soiling photographs, electrical layouts
(without matrix panel E) and the simulation framework into
`results/1_validation/fig_soiling_layout_framework.pdf` and `.png`.
The existing photographic views form panel a at upper left; electrical layouts
form panel b at upper right, and the framework spans the bottom (panel c).
The original annotated photographs are clipped from `data/figure_assets/Soiling-images.png`;
no photograph content is synthesised or retouched. Use `--without-photos` to
reproduce the previous two-panel layout.
The source framework is kept in `data/figure_assets/`; the electrical layouts are regenerated
by `soilspice/schematics.py`. PDF panels remain vector in the assembled figure (183 × 148 mm).

`python Nature-builder-panel.py --framework-only` exports the simulation framework alone
to `results/1_validation/fig_simulation_framework.pdf` (vector) and `.png` (600 dpi),
at 183 mm wide with 3 mm margins and no outer panel letter.

`python Nature-builder-validation.py` creates `results/1_validation/fig_validation_panel.pdf` and `.png`:
rectangle I-V curves across the top (a), curve errors at bottom left (b), and model comparison at bottom right (c).
The saved results are redrawn in compact formats to fit a 183 × 170 mm figure,
with matching bottom panels, 11 pt bold panel letters and vector plots.

`python Nature-builder-regimes.py` creates `results/2_simulation/fig_regimes_panel.pdf` and `.png`
at 183 × 170 mm: an integrated response map above (a), shared legends, and an atlas below (b)
with modules in rows and band configurations in columns. The drawing code is `soilspice/fig_regimes.py`.
It uses the saved `resilience_states.csv`; the standalone map and atlas retain their existing layout.

`optical_loss_modules_grid` shows positive power losses, mean curves in 1 % optical-loss
bins and min–max envelopes. The Y-axis label includes the definition
`100 × (Pmp,clean − Pmp,soiled) / Pmp,clean [%]`, without an `L_P` symbol. The positive loss
scale runs from 0 % at the bottom to 100 % at the top; no extra line is drawn at zero.
Dotted horizontal guides
labelled `Bypass loss reference` mark the power share of
the smallest bypass-protected substring: 1/3 (33.3 %) for Yingli, JA Solar and LONGi,
and 1/4 (25 %) for SunPower's unbalanced 2/4/2-column layout. These are idealised loss
levels associated with bypassing a substring, not exact activation thresholds or fitted
plateau values; the simulated losses can lie above these references.

Chapter 2 exports `average_power_windows_horizontal.pdf` and `.png`: the retained
average-power values and scale limits, with the 16 module/configuration combinations
in columns and the eight optical-loss windows in rows, with a vertical colour bar
on the right, in a wide format (183 × 65 mm). It uses the original grey–gold–green
palette, with all numbers in black. The colour scale remains 0.5–1.0.

`python chapters/fig_resilience.py --horizontal-power-only` regenerates this figure
from the saved metrics and checks all 128 values against the source data. This is the
only retained standalone average-power-window layout; trial layouts, palettes and
comparison exports have been removed.

## Energy and economics: retained rain-event model

The illustrative reference is **R_H = 15 mm/event**. The sensitivity scenarios are **5, 10, 20 and 25 mm/event**,
with identical weather, deposition, geometry and economic assumptions. No threshold has been calibrated for bands.

```
python chapters/4_economics.py                      # full reference scenario, R_H=15 mm/event
python chapters/4_economics.py --washing-total 10   # one sensitivity scenario
python chapters/4_economics.py --all-thresholds     # all five scenarios and comparisons
python chapters/4_economics.py --comparison-only    # all five comparison outputs
python chapters/4_economics.py --band-growth-only   # reference temporal trajectories
```

Events end after **three consecutive dry hours**. One or two dry hours, midnight and the change of year do not
split an event. Each completed event causes one update, at closure:

1. P_e >= R_H: assumed complete washing, resetting both uniform and band deposits.
2. Otherwise, P_e >= 1 mm: the existing mass balance transfers dust to the band.
3. Otherwise, rain has no effect. Deposition continues between event updates.

The criterion uses **event total precipitation**, not hourly peak intensity. Source hourly peaks remain recorded
in the event catalogue as descriptive data. Incomplete trailing events cause no artificial update; missing records
cannot count as dry. The Open-Meteo precipitation column is used as supplied. All event times are UTC.

Reference figures/tables and source CSVs: `results/4_economics/total_scenarios/PH15mm/`.
The four other scenarios have identical folder structure under `total_scenarios/`. Its README links all five sets.
Shared supplementary figures: `fig_rain_total_comparison` and `fig_band_growth_total_comparison` (PDF/PNG).
`data/rain_total_*.csv` stores event counts, explicitly censored washing intervals, annual and per-configuration
energy/transmission metrics and temporal trajectories for all five thresholds. `data/rain_events.csv` stores
source events; `data/rain_events_total_scenarios.csv` stores their classifications.

`python chapters/rain_event_audit.py` independently partitions positive hourly rainfall, checks mass conservation,
and verifies one action per completed event for each of the five thresholds (`data/rain_total_event_audit.json`).
The current automated tests cover event separation, mass balance, integration, independent economic optimisation,
signed manual-cleaning gains and threshold robustness (`python -m unittest discover -s tests`).

Retired hourly-intensity scenarios, daily washing trials and MeteoSwiss/ERA5 source/threshold-comparison experiments
have been removed, together with their dedicated scripts and tests. Raw climate data are retained. The cleanup
leaves numerical outputs of the five event-total scenarios and the other three chapters unchanged.

`python chapters/summarize_economics.py` summarises the five saved total-rain scenarios without rerunning the model.
It exports `results/4_economics/data/threshold_sensitivity/` (reference and full-grid energy/cost tables, independent
with/without-band optima, signed manual-cleaning gains, +10 mm margin net value, censored washing intervals,
per-configuration robustness and adjacent-threshold reversals). `README.md`, `conclusions_reference.md` and
`conventions.json` in that folder define units, annualisation, counterfactuals and policy sets. The reference is
h=60 mm, f=0.1; both markets and all saved heights/fractions remain in the tables. The two figures
`fig_energy_economics_threshold_sensitivity.pdf` and `_rooftop.pdf` show every module/orientation, without confidence
intervals; R_H=15 mm/event is the illustrative reference. Existing five-scenario temporal trajectories remain supplementary figures.

`python chapters/fig_washing_threshold_sensitivity.py` creates **additional** 3 × 3 figures for the ground and rooftop
markets: `fig_washing_threshold_sensitivity_ground.pdf`/`.png` and `_rooftop.pdf`/`.png` (183 × 170 mm).
The third row shows the signed manual-cleaning gain for total soiling, excluding none/robot from the manual-policy
selection. A vertical line marks the user-specified illustrative reference R_H=15 mm/event. New plotted CSVs,
captions, sign-change brackets and checks are in `data/washing_threshold_sensitivity_3x3/`. Values reconcile with
the saved summary and raw runs; all retained figures, including the main economics figures, are protected.

## Requirements

### Yingli irradiance-dependent shunt comparison

`python chapters/yingli_shunt_comparison.py` compares every available Yingli flash curve
(rectangle/triangle, portrait/landscape) with constant and PVsyst exponential cell shunt resistance.
All modules use constant shunt resistance by default, including Yingli (5.99 ohm/cell).
Only this separate comparison script explicitly enables `PVsystShuntLaw(5.99)` for its diagnostic variant.
The law uses the same per-cell effective irradiance as the photocurrent, including partial coverage.
Plots, paired curve metrics, cell-level resistances and provenance are isolated in
`results/1_validation/yingli_shunt_comparison/`. Existing publication exports are preserved and
are not regenerated by this command; the PVsyst variant is not used by the main simulation pipeline.

Python >= 3.9, `pip install -r requirements.txt`, and LTspice (only with `--full`). Set the path of the flash-test
files and of LTspice in `soilspice/config.py`.

The flash-test campaign (300 MB of .asc files) is not in the repository: `DATA_DIR` in `config.py` points to it.
Without `--full`, chapters 1 and 2 start from the saved simulations in `results/1_validation/data/` and
`results/2_simulation/data/`.

After changing the campaign transmittances in `soilspice/experiments.py`, run
`python chapters/1_validation.py --full` to recompute validation simulations, errors and figures.
Chapter 1 checks the saved configurations and transmittances before reusing them and asks for `--full`
if they differ from the current campaign.
