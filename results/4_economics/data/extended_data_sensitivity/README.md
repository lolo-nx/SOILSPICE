# Extended Data sensitivity figures

Reproduce: `MPLBACKEND=Agg python chapters/extended_data_sensitivity.py`.
This script reads existing corrected simulations; it runs no new physical model.
Figures: `../../fig_ed1_parameter_sensitivity.pdf` and
`../../fig_ed2_washing_threshold_sensitivity.pdf`, with matching 600 dpi PNGs.
Both vector PDFs have exact dimensions 183 × 125 mm, editable embedded fonts and six panels.
The earlier figures, source results and manuscript are unchanged.

## Source data and conventions

- `ed1_parameter_cases.csv`: all 33 site/case medians used for Fig. 1, including
  reference values, actual inputs and units. The five parameters vary one at a time
  at RH=15 mm/event. Alpha changes with D held fixed; a change of soiling rate
  changes D according to the existing mass conversion. Blue/orange bars terminate
  at the lower/higher case median; both extend from the same black reference line.
- `ed2_configurations.csv`: all 120 ground-market configuration/threshold rows,
  including the full original manual-cleaning results and signed gains, whether
  or not plotted. No policy is reselected by this figure script.
- `ed2_plot_data.csv`: all 15 site/threshold groups, each with exactly eight
  configurations. Minimum, median and maximum are calculated from the individual
  configuration metrics, not from differences of separately aggregated quantities.
  The manual-cleaning aggregates remain in this file but are not plotted.
- `reference_checks.csv`: independent comparisons of Fig. 1, Fig. 2 and raw policy
  runs at RH=15. The energy metric is comparable; the two economic metrics have
  different definitions and units and must not be equated.
- `caption_ed1.tex` and `caption_ed2.tex`: standalone LaTeX captions.
- `verification.json`: source hashes, data reconciliation and output checks.

Reference geometry: 60 mm rectangular band, f=0.1, no initial added glass margin.
Ground-market prices: 0.06 EUR/kWh and 0.41 EUR/kWp per manual intervention.
Energy percentages are extra lost energy with band minus without band under
natural rainfall only, divided by the identical clean production. They are not
total-soiling losses or percentages of soiled production.
Fig. 1 cost percentages compare the two independently simulated trigger policies
with and without band, including cleaning charges, divided by clean revenue.
They do not minimise over fixed-interval strategies or include a robot.
Fig. 2 margin NPV uses the existing independently optimised costs with/without
band at each margin. It equals the 25-year annuity factor at 5% times the annual
reduction in additional band cost, less the added-glass investment
(100 EUR/kWp × 10 mm / module length along the slope). The existing choice among
no manual cleaning, trigger and fixed intervals is preserved; robot excluded.
The negative values of net benefit or manual-cleaning gain are never clipped.

## Timing and washing

The corrected regeneration manifest is checked against the current rain/energy
code, all weather and SPICE inputs, and the source exports before plotting.
Open-Meteo hour-end labels are shifted once to interval starts after PV yield
calculation. Three consecutive dry hours close an event, including across
midnight/year boundaries; one rain action occurs at closure. Decisions use only
the completed event total, never intensity: total >= RH assumes full removal,
otherwise total >= RL transfers dust; smaller totals have no rain action.
UTC daily manual decisions and all cleaning assumptions remain unchanged.
The physical source coverage ends at 2025-12-31 23:00 UTC: the last source hour
is absent and no synthetic weather is added. Existing annualisation by actual
duration / 365.25 and clean-yield ratios are retained.

RH=15 mm/event is an illustrative reference, not selected from these curves.
Complete washing is an assumption; none of the washing thresholds is calibrated
for actual bottom-edge bands. Shaded ranges describe variation across the eight
configurations, not confidence intervals. Lines connect only the five simulated
thresholds; no additional scenario or exact crossing threshold is inferred.
Fig. 1 horizontal scales vary by site. Fig. 2 panel a has a separate energy scale;
panels b–c share a scale, and d–f share a scale.

All comparable RH=15 references and all raw-run/synthesis-table identities agree
within floating-point tolerance (1e-11 percentage points). No inconsistency found.
