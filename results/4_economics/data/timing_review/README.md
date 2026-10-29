# Rain timing and cleaning review

This review changes one calculation: Open-Meteo's hourly precipitation and radiation rows label the **end** of
the preceding hour. `site_hours()` now shifts the joined climate/yield frame back one hour, after `pv_yield()`
has evaluated the sun at the original hour midpoint. Rain and energy therefore share physical interval starts.
The solar calculation, physical parameters, SPICE model, event-total thresholds and cleaning strategies are unchanged.
[Source timestamp definitions](https://open-meteo.com/en/docs/historical-weather-api#hourly-parameter-definition).

## Checks and correction

- Events require three consecutive zero-precipitation hours to close. One/two dry hours, midnight and year
  boundaries do not split them. Incomplete trailing events cause no action. The reference remains 15 mm/event;
  sensitivity thresholds remain 5, 10, 20 and 25 mm/event. Hourly intensity is not used for washing.
- Every completed event maps to exactly one simulation action, at the end of its third dry hour. This was
  checked against all three source series and the five thresholds; see `event_audit.csv`.
- Maintenance decisions remain at UTC midnight, once per day. The preceding interval's rain action occurs
  before that decision. Natural washing at midnight can therefore cancel an otherwise due manual intervention.
  The source CSV headers all specify GMT with a zero UTC offset. Local daily boundaries were not introduced.
- Manual cleaning and complete natural washing reset both deposits, `cum`, `cumb` and the fixed-interval age.
  Intermediate transfer removes the uniform deposit and adds the trapped fraction to the band, retaining its age.
- `cum` accumulates total losses since removal, adjusted after transfer to retain `cumb`, the accumulated
  **additional band loss**. Repeated transfers retain that band counter correctly. The historical `lost` energy
  is never erased. No trigger-rule change was needed.
- The raw robot rows deliberately contain lost-revenue costs without investment. `fig_cleaning` adds the
  existing annual investment range **1.104–3.772 EUR/kWp/year** once, to both band and uniform-only costs.
  Its centre is 2.438 EUR/kWp/year. Recalculation reproduces all three previously saved robot comparison values;
  see `robot_cost_check.csv`. Manual-policy optimisation continues to exclude the robot.
- The 29 passing tests include seven focused additions in `tests/test_energy_timing.py`: timestamp conversion,
  year-crossing rain, midnight ordering, trapped-mass conservation, repeated-transfer trigger counters,
  natural/manual clock resets, and daily robot decisions despite extra event splits. The existing hourly versus
  grouped integration tests cover all four strategies. The analytic trigger test also checks historical energy loss.

## Before/after reference results

Both complete R_H=15 reference runs contain **67,392 rows**, spanning the existing geometry/parameter/policy grids
at all three sites and both markets. The freshly run baseline matches the pre-correction saved reference.
The simulation, cost, solar-yield and event-detection function bodies were also checked to be unchanged apart
from documentation. Only the input-time alignment changed.

For the rectangular reference band **h=60 mm, f=0.1, initial added margin 0**:

- Maximum absolute change in additional band energy loss without manual cleaning: **0.004428 pp of clean yield**.
- Maximum absolute change in manual-cleaning gain: **0.049560 pp of clean revenue**.
- Maximum absolute change in the 25-year NPV of +10 mm margin: **0.859244 EUR/kWp**.
- Intervention counts change by at most **one intervention over the decade** for an individual reference
  policy/configuration. The maximum annualised frequency change is **0.100119 interventions/year**.
- No manual-cleaning gain changes sign. One margin NPV does: **Doha / ground / Yingli / landscape**, from
  **−0.284047 to +0.176863 EUR/kWp**. This remains a small, timing-sensitive economic benefit.

Across the entire original grid (all heights, fractions, margins and policies), the largest total-energy-loss
change is **0.078559 pp**. `full_grid_change_summary.csv` retains results by site, market and strategy.
These are numerical differences between two timings, not uncertainty or confidence intervals.

## Files and publication status

- `reference_model_runs_after.csv`: corrected full reference, with the existing raw cost convention for robot rows.
- `reference_growth_series_after.csv`, `reference_site_table_after.csv`: corrected reference states and site statistics.
- `reference_energy_cleaning_comparison.csv`: before/after losses, costs and actual/annualised cleaning counts
  for every reference module/orientation/policy, both markets.
- `reference_economic_comparison.csv`: independent optimisation, signed manual gains and +10 mm margin NPV,
  before and after, for all 48 site/market/module/orientation combinations.
- `coverage.csv`, `event_audit.csv`, `robot_cost_check.csv`, `verification.json`: detailed verification.

**Publication status: all five scenarios and their comparison figures have now been regenerated with the corrected timestamp handling.** The published R_H=15 model runs reproduce the corrected reference in this review. The before/after columns in this folder remain historical review evidence; current figures and numerical tables are in `total_scenarios/` and the two sensitivity folders. See the [regeneration verification](../corrected_timestamp_regeneration.json).

The complete regeneration commands are:

```sh
python chapters/4_economics.py --all-thresholds
python chapters/summarize_economics.py
python chapters/fig_washing_threshold_sensitivity.py
python chapters/rain_event_audit.py
```

## Limitations to state in the paper

The trigger is a retrospective, rain-adjusted loss heuristic, not an optimisation of forecast recoverable revenue.
Resetting its uniform component after transfer does not recover historical lost energy and does not reconstruct
the monetary history of the fraction of dust retained in the band. Fixed intervals and triggers are checked once
per UTC day; decisions are not scheduled from a rain forecast. Daily robotic cleaning remains idealised.

Rain acts once at event closure, including its three dry hours; deposition continues between updates, including
during rainfall. Intermediate rainfall transfers dust and complete washing removes both deposits under the existing
assumptions. Complete washing and the event-total thresholds remain uncalibrated for actual bottom-edge bands.
The input column includes total precipitation; snow processes are not modelled separately. Initial deposits are clean.

The source files do not contain the hour ending at 2026-01-01 00:00 UTC. After correcting labels, coverage ends
at **2025-12-31 23:00 UTC**, one hour before the nominal ten-year endpoint. No weather is invented to fill it.
The total clean energy is unchanged because the excluded old boundary hour had zero production at all sites.
Economic annualisation still uses covered duration / 365.25, so its denominator changes from 3653 days to
3652 days 23 hours. Site-table annual yield remains total / ten calendar years; see `coverage.csv`.

No Methods text was rewritten.
