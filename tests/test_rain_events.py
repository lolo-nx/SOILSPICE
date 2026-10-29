"""Event separation, cross-midnight totals and equivalence to an hourly mass/energy calculation."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from soilspice import energy, rain


def hourly(values, start="2020-01-01"):
    return pd.Series(values, index=pd.date_range(start, periods=len(values), freq="h", tz="UTC"))


class RainEventTests(unittest.TestCase):
    def test_one_or_two_dry_hours_do_not_split_an_event(self):
        r = hourly([1, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0])
        e = rain.events(r)
        self.assertEqual(len(e), 2)
        self.assertEqual(e.rain.tolist(), [3, 3])
        self.assertEqual(e.end.iloc[0], r.index[7])
        self.assertEqual(e.start.iloc[1], r.index[7])
        self.assertTrue(e.complete.all())

    def test_midnight_does_not_split_a_washing_event(self):
        r = hourly([4, 4, 4, 4, 4, 4, 1, 0, 0, 0], start="2020-01-01 20:00")
        e = rain.events(r)
        self.assertEqual(len(e), 1)
        self.assertEqual(e.rain.iloc[0], 25)
        self.assertEqual(e.rain_max.iloc[0], 4)
        self.assertTrue(energy.washing_flags(e.rain)[0])
        d = rain.steps(pd.DataFrame({"rain": r, "Y": np.ones(len(r))}), e)
        self.assertEqual(d.rain.sum(), 25)
        self.assertEqual((d.event_id >= 0).sum(), 1)
        self.assertTrue((d.index.hour == 0).any())
        self.assertEqual(d.end_time.iloc[-1], r.index[-1] + pd.Timedelta(hours=1))

    def test_incomplete_trailing_event_is_not_artificially_closed(self):
        r = hourly([0, 2, 0, 0])
        e = rain.events(r)
        self.assertFalse(e.complete.iloc[0])
        self.assertTrue(pd.isna(e.end.iloc[0]))
        d = rain.steps(pd.DataFrame({"rain": r, "Y": np.ones(len(r))}), e)
        self.assertEqual(d.rain.sum(), 0)
        self.assertTrue((d.event_id == -1).all())

    def test_dry_period_has_no_events(self):
        r = hourly(np.zeros(48))
        e = rain.events(r)
        self.assertTrue(e.empty)
        d = rain.steps(pd.DataFrame({"rain": r, "Y": np.ones(48)}), e)
        np.testing.assert_allclose(d.duration_days, [1, 1])
        np.testing.assert_allclose(d.Y, [24, 24])

    def test_fractional_hourly_rain_reaches_exact_event_thresholds(self):
        for amount, hours, total in ((0.1, 10, 1.0), (0.1, 250, 25.0)):
            with self.subTest(total=total):
                e = rain.events(hourly([amount] * hours + [0, 0, 0]))
                self.assertEqual(e.rain.iloc[0], total)

    def test_missing_rain_or_hours_cannot_count_as_dry(self):
        for values in ([1, np.nan, 0, 0], [1, -0.1, 0, 0]):
            with self.assertRaises(ValueError):
                rain.events(hourly(values))
        with self.assertRaisesRegex(ValueError, "consecutive hourly"):
            rain.events(hourly([1, 0, 0, 0]).drop(pd.Timestamp("2020-01-01 01:00", tz="UTC")))

    def test_grouped_integration_matches_hourly_simulation(self):
        values = np.zeros(48)
        values[[1, 3]] = 1
        values[21:28] = 4
        values[40] = 0.5
        values[47] = 2  # right-censored; no fabricated end-of-event action
        r = hourly(values)
        h = pd.DataFrame({"rain": r, "Y_JASOLAR": 0.1 + np.sin(np.arange(48)) ** 2})
        e = rain.events(r)
        d = rain.steps(h, e)
        D, alpha = 0.4, 0.2
        gain = np.array([2.0])
        loss_table = (0.3 * (1 - energy.TAU))[None, :]
        event_rain = np.zeros(48)
        for event in e[e.complete].itertuples():
            last = h.index.searchsorted(event.end) - 1
            event_rain[last] = event.rain
        energy.attenuated_yields.cache_clear()
        with patch.object(energy, "site_hours", return_value=h), patch.object(energy, "site_steps", return_value=d):
            dust = energy.attenuated_yields("test", D, alpha)
        energy.attenuated_yields.cache_clear()
        y = h[["Y_JASOLAR"]].to_numpy()
        for mode in ("none", "trigger", "interval", "robot"):
            with self.subTest(mode=mode):
                fine = energy.simulate(y, event_rain, D, gain, loss_table, mode, k=1, price=1, cost=0.01,
                                       washing_total=15.0,
                                       alpha=alpha, duration_days=np.full(48, 1 / 24),
                                       dust_y=y * np.exp(-alpha * D / 24), clean_start=h.index.hour == 0, track=0)
                grouped = energy.simulate(d[["Y_JASOLAR"]].to_numpy(), d.rain.to_numpy(), D, gain, loss_table,
                                          mode, k=1, price=1, cost=0.01, alpha=alpha,
                                          washing_total=15.0,
                                          duration_days=d.duration_days, dust_y=dust.to_numpy(),
                                          clean_start=d.clean_start, track=0)
                for actual, expected in zip(grouped[:3], fine[:3]):
                    np.testing.assert_allclose(actual, expected, atol=1e-12)
                np.testing.assert_allclose(np.array(grouped[3])[:, :2],
                                           np.array(fine[3])[d.hour_start.to_numpy(), :2], atol=1e-12)


if __name__ == "__main__":
    unittest.main()
