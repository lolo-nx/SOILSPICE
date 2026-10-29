"""Audit hourly rain partition and single event actions for the five event-total scenarios.

    python chapters/rain_event_audit.py
"""
import json
import sys
from math import fsum
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from soilspice import config, energy, rain

START = pd.Timestamp(f"{energy.YEARS[0]}-01-01", tz="UTC")
STOP = pd.Timestamp(f"{energy.YEARS[1] + 1}-01-01", tz="UTC")
YEARS = list(range(energy.YEARS[0], energy.YEARS[1] + 1))


def audit_catalogue(saved, growth, threshold):
    """Independently partition positive-rain hours, then match source, catalogue and simulation actions."""
    reports = []
    for site, filename in energy.FILES.items():
        hourly, *_ = energy.era5(filename)
        # Audit physical intervals: Open-Meteo's source labels are hour ends.
        hourly.index = hourly.index - pd.Timedelta(hours=1)
        detected = rain.events(hourly.rain)
        wet = np.flatnonzero(hourly.rain.to_numpy() > 0)
        # Three zero records between positive records means a separation of >= 4 hours.
        groups = np.split(wet, np.flatnonzero(np.diff(wet) > rain.DRY_HOURS) + 1) if len(wet) else []
        assert len(groups) == len(detected), site
        for positions, event in zip(groups, detected.itertuples()):
            amounts = hourly.rain.iloc[positions]
            assert event.start == amounts.index[0]
            assert event.last_rain == amounts.index[-1]
            assert event.wet_hours == len(positions)
            np.testing.assert_allclose(event.rain, fsum(amounts), rtol=0, atol=1e-12)
            np.testing.assert_allclose(event.rain_max, amounts.max(), rtol=0, atol=1e-12)
            if event.complete:
                end = event.last_rain + pd.Timedelta(hours=1 + rain.DRY_HOURS)
                assert event.end == end
                assert (hourly.rain.loc[event.last_rain + pd.Timedelta(hours=1):
                                        event.end - pd.Timedelta(hours=1)] == 0).all()
                assert event.duration_hours == (end - event.start).total_seconds() / 3600
        assert sum(map(len, groups)) == len(wet)
        np.testing.assert_allclose(fsum(detected.rain), fsum(hourly.rain), rtol=0, atol=1e-9)
        expected = detected[(detected.start < STOP) & (detected.end.isna() | (detected.end > START))]
        actual = saved[saved.site == site].set_index("event_id").sort_index()
        assert actual.index.is_unique and not actual.start.duplicated().any()
        assert not actual.loc[actual.complete, "end"].duplicated().any()
        pd.testing.assert_frame_equal(actual[expected.columns], expected, check_dtype=False, check_freq=False)
        complete = expected[expected.complete & (expected.end > START) & (expected.end <= STOP)]
        study = hourly[(hourly.index >= START) & (hourly.index < STOP)].copy()
        study["Y"] = 0.0  # Only the event/midnight grid is needed for this audit.
        steps = rain.steps(study, expected)
        actions = steps[steps.event_id >= 0].set_index("event_id").sort_index()
        assert actions.index.is_unique
        np.testing.assert_array_equal(actions.index, complete.index)
        np.testing.assert_allclose(actions.rain, complete.rain, rtol=0, atol=1e-12)
        pd.testing.assert_series_equal(actions.end_time, complete.end, check_names=False)
        stored = growth[(growth.site == site) & (growth.event_id >= 0)].set_index("event_id").sort_index()
        assert stored.index.is_unique
        np.testing.assert_array_equal(stored.index, complete.index)
        np.testing.assert_allclose(stored.rain, complete.rain, rtol=0, atol=1e-12)
        np.testing.assert_array_equal(pd.to_datetime(stored.end_time, utc=True), complete.end)
        expected_washing = complete.rain >= threshold
        np.testing.assert_array_equal(stored.washing, expected_washing)
        reports.append(dict(site=site, family="total", threshold=threshold, catalogue_events=len(actual), complete_actions=len(actions),
                            duplicate_event_ids=int(actual.index.duplicated().sum()),
                            duplicate_simulation_actions=int(stored.index.duplicated().sum()),
                            independent_hourly_partition_matches=True, rain_mass_conserved=True,
                            washing_actions=int(stored.washing.sum())))
    return reports


def main():
    from soilspice.rain_scenarios import scenario_tag
    events = pd.read_csv(config.ECON_DIR / "rain_events.csv")
    for column in ("start", "last_rain", "end"):
        events[column] = pd.to_datetime(events[column], utc=True)
    reports = []
    thresholds = energy.TOTAL_THRESHOLDS
    for threshold in thresholds:
        growth = pd.read_csv(config.CH4_DIR / "total_scenarios" / scenario_tag(threshold)
                             / "data" / "growth_series_reference.csv")
        reports.extend(audit_catalogue(events, growth, threshold))
    target = config.ECON_DIR / "rain_total_event_audit.json"
    target.write_text(json.dumps(reports, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(reports, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
