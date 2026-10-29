"""Rain events from hourly records and a simulation grid ending at event closures and midnight."""
from math import fsum

import numpy as np
import pandas as pd

DRY_HOURS = 3
RAIN_INTERVAL_HOURS = 1.0  # validated consecutive hourly records


def events(hourly_rain, dry_hours=DRY_HOURS):
    """A positive hourly accumulation opens an event; dry_hours consecutive zero hours close it.
    Each timestamp labels the START of its one-hour rain interval. Closure is after the third dry interval.
    Open-Meteo's end-labelled inputs are converted once in energy.site_hours, before this function is called.
    rain_max is max(hourly accumulation / 1 h), in mm/h, not the event mean.
    One/two dry hours and midnight do not split an event. A trailing open event stays incomplete.
    Missing/negative rain or missing hours cannot count as dry weather.
    """
    if not isinstance(dry_hours, int) or dry_hours < 1:
        raise ValueError("dry_hours must be a positive integer")
    index = hourly_rain.index
    if not isinstance(index, pd.DatetimeIndex) or not index.is_monotonic_increasing or index.has_duplicates:
        raise ValueError("hourly rain must have an ordered, unique DatetimeIndex")
    if len(index) > 1 and not np.all(np.diff(index.asi8) == pd.Timedelta(hours=1).value):
        raise ValueError("rain events require consecutive hourly records")
    values = hourly_rain.to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("hourly rain must be finite and nonnegative")
    rows, active, dry, amounts = [], None, 0, []
    for time, amount in zip(index, values):
        if amount > 0:
            if active is None:
                active = dict(start=time, last_rain=time, rain=0.0, rain_max=0.0, wet_hours=0)
                amounts = []
            active["last_rain"] = time
            active["rain"] += amount
            amounts.append(amount)
            hourly_intensity = amount / RAIN_INTERVAL_HOURS  # mm / h
            active["rain_max"] = max(active["rain_max"], hourly_intensity)
            active["wet_hours"] += 1
            dry = 0
        elif active is not None:
            dry += 1
            if dry == dry_hours:
                end = time + pd.Timedelta(hours=1)
                active["rain"] = fsum(amounts)
                rows.append({**active, "end": end, "complete": True,
                             "duration_hours": (end - active["start"]) / pd.Timedelta(hours=1)})
                active, dry = None, 0
    if active is not None:
        active["rain"] = fsum(amounts)
        rows.append({**active, "end": pd.NaT, "complete": False,
                     "duration_hours": (index[-1] + pd.Timedelta(hours=1) - active["start"]) / pd.Timedelta(hours=1)})
    columns = ["start", "last_rain", "rain", "rain_max", "wet_hours", "end", "complete", "duration_hours"]
    result = pd.DataFrame(rows, columns=columns)
    result.index.name = "event_id"
    for column in ("start", "last_rain", "end"):
        result[column] = pd.to_datetime(result[column], utc=index.tz is not None)
    result["complete"] = result.complete.astype(bool)
    return result


def steps(hourly, event_table):
    """Combine consecutive hours without an event closure or a daily cleaning decision inside them.
    Sum hourly yields in each interval; event totals/peaks are present only on the interval closing that event.
    The hourly index labels interval starts. At a midnight closure, the rain action belongs to the preceding
    interval, followed by one manual decision at the next interval's midnight start.
    This grouping leaves all event and maintenance times intact.
    """
    index = hourly.index
    if not len(index):
        raise ValueError("no hourly climate data in the simulation period")
    stop = index[-1] + pd.Timedelta(hours=1)
    closed = event_table[event_table.complete & (event_table.end > index[0]) & (event_table.end <= stop)]
    event_stops = np.searchsorted(index.asi8, pd.DatetimeIndex(closed.end).asi8)
    midnight = np.flatnonzero(index.hour == 0)
    stops = np.unique(np.r_[midnight[midnight > 0], event_stops, len(index)]).astype(int)
    starts = np.r_[0, stops[:-1]]
    columns = [c for c in hourly if c == "Y" or c.startswith("Y_")]
    result = pd.DataFrame(np.add.reduceat(hourly[columns].to_numpy(), starts, axis=0),
                          index=index[starts], columns=columns)
    result["hour_start"], result["hour_stop"] = starts, stops
    result["duration_days"] = (stops - starts) / 24
    result["clean_start"] = index[starts].hour == 0
    result["end_time"] = index[stops - 1] + pd.Timedelta(hours=1)
    result["rain"], result["rain_max"], result["event_id"] = 0.0, 0.0, -1
    if len(closed):
        positions = np.searchsorted(stops, event_stops)
        result.iloc[positions, result.columns.get_loc("rain")] = closed.rain.to_numpy()
        result.iloc[positions, result.columns.get_loc("rain_max")] = closed.rain_max.to_numpy()
        result.iloc[positions, result.columns.get_loc("event_id")] = closed.index.to_numpy()
    return result
