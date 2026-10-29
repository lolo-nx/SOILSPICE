"""Physical hour labels, rainfall/maintenance ordering and trigger bookkeeping."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from soilspice import energy, rain


class EnergyTimingTests(unittest.TestCase):
    def test_source_hour_ends_are_shifted_once_after_pv_calculation(self):
        raw = pd.DataFrame({'rain': [15., 0., 0., 0.]},
                           index=pd.date_range('2020-12-31 21:00', periods=4, freq='h', tz='UTC'))
        production = pd.DataFrame({'Y': [.1, .2, .3, .4]}, index=raw.index)
        energy.site_hours.cache_clear()
        self.addCleanup(energy.site_hours.cache_clear)
        with patch.object(energy, 'era5', return_value=(raw, 0, 0, 0)), \
             patch.object(energy, 'pv_yield', return_value=production) as calculate:
            hours = energy.site_hours('Doha')
        # Solar calculation must still receive original end labels; energy values do not change.
        pd.testing.assert_index_equal(calculate.call_args.args[0].index, raw.index)
        pd.testing.assert_index_equal(hours.index, raw.index - pd.Timedelta(hours=1))
        np.testing.assert_array_equal(hours.Y, production.Y)
        event = rain.events(hours.rain).iloc[0]
        self.assertEqual(event.start, pd.Timestamp('2020-12-31 20:00', tz='UTC'))
        self.assertEqual(event.end, pd.Timestamp('2021-01-01 00:00', tz='UTC'))

    def test_year_boundary_does_not_split_one_event(self):
        amounts = pd.Series([4., 0., 0., 11., 0., 0., 0.],
                            index=pd.date_range('2020-12-31 22:00', periods=7, freq='h', tz='UTC'))
        events = rain.events(amounts)
        self.assertEqual(len(events), 1)
        self.assertEqual(events.rain.iloc[0], 15.)
        self.assertEqual(events.end.iloc[0], pd.Timestamp('2021-01-01 05:00', tz='UTC'))

    def test_midnight_washing_precedes_due_fixed_interval_cleaning(self):
        hours = pd.DataFrame({'rain': np.zeros(72), 'Y': np.ones(72)},
                             index=pd.date_range('2020-12-30', periods=72, freq='h', tz='UTC'))
        hours.iloc[44, hours.columns.get_loc('rain')] = 15.  # Dec 31, 20:00--21:00
        events = rain.events(hours.rain)
        steps = rain.steps(hours, events)
        self.assertEqual(events.end.iloc[0], pd.Timestamp('2021-01-01', tz='UTC'))
        np.testing.assert_array_equal(steps.rain, [0, 15, 0])
        args = (steps.Y.to_numpy(), steps.rain.to_numpy(), 1., np.ones(1),
                np.zeros((1, len(energy.TAU))))
        kwargs = dict(mode='interval', k=2, duration_days=steps.duration_days,
                      clean_start=steps.clean_start)
        # A wash at the due midnight resets the age before the manual decision.
        self.assertEqual(energy.simulate(*args, **kwargs)[1][0], 0)
        self.assertEqual(energy.simulate(args[0], np.zeros(3), *args[2:], **kwargs)[1][0], 1)

    def test_transfer_mass_conserves_the_trapped_fraction_and_wash_resets_both_layers(self):
        fraction, length, height, deposition = .1, 1000., 50., 2.
        gain = fraction * length / height
        intervals = np.array([.25, .75, .5, .5, .5])
        _, _, _, states = energy.simulate(
            np.ones(5), np.array([0, 1, 1, 15, 0]), deposition, np.array([gain]),
            np.zeros((1, len(energy.TAU))), duration_days=intervals,
            clean_start=np.zeros(5, dtype=bool), track=0)
        tau, _, uniform_loss = np.asarray(states).T
        band_mass = -np.log(tau) / energy.ALPHA
        np.testing.assert_allclose(band_mass, [0, 0, 4, 6, 0], atol=1e-12)
        # Mass per total module area in the band = trapped fraction of the transferred uniform mass.
        self.assertAlmostEqual(band_mass[2] * height / length, fraction * deposition * intervals[:2].sum())
        self.assertAlmostEqual((band_mass[3] - band_mass[2]) * height / length,
                               fraction * deposition * intervals[2])
        self.assertAlmostEqual(uniform_loss[-1], 1 - np.exp(-energy.ALPHA * deposition * intervals[-1]))

    def test_repeated_transfers_retain_band_counter_without_recovering_historical_energy(self):
        # alpha=ln(2) halves transmission for one unit of deposited mass; all LUT points here are exact.
        lost, cleanings, _, states = energy.simulate(
            np.ones(4), np.array([1, 1, 1, 0]), 1., np.ones(1), (1 - energy.TAU)[None, :],
            mode='trigger', price=1., cost=.6, alpha=np.log(2), track=0)
        # After consecutive transfers: cum = cumb = 0, .25, .625, hence only the fourth step cleans.
        np.testing.assert_allclose(np.asarray(states)[:, 0], [1, .5, .25, 1])
        self.assertEqual(cleanings[0], 1)
        self.assertAlmostEqual(lost[0], .5 + .75 + .875 + .5)

    def test_natural_and_manual_resets_restart_trigger_and_interval_clocks(self):
        for mode in ('trigger', 'interval'):
            with self.subTest(mode=mode):
                lost, cleanings, _, states = energy.simulate(
                    np.ones(6), np.array([0, 15, 0, 0, 0, 0]), 1., np.ones(1),
                    np.zeros((1, len(energy.TAU))), mode=mode, k=2,
                    price=1., cost=.9, alpha=np.log(2), track=0)
                self.assertEqual(cleanings[0], 1)
                np.testing.assert_allclose(np.asarray(states)[:, 2], [.5, .75, .5, .75, .5, .75])
                self.assertAlmostEqual(lost[0], 3.75)

    def test_multiple_event_splits_do_not_add_daily_robot_decisions(self):
        amounts = np.zeros(72)
        amounts[[3, 12, 27, 39, 52, 65]] = 1.
        hours = pd.DataFrame({'rain': amounts, 'Y': np.ones(72)},
                             index=pd.date_range('2020-12-31', periods=72, freq='h', tz='UTC'))
        steps = rain.steps(hours, rain.events(hours.rain))
        self.assertGreater(len(steps), 3)
        self.assertEqual(steps.clean_start.sum(), 3)
        result = energy.simulate(steps.Y.to_numpy(), steps.rain.to_numpy(), 1., np.ones(1),
                                 np.zeros((1, len(energy.TAU))), mode='robot',
                                 duration_days=steps.duration_days, clean_start=steps.clean_start)
        self.assertEqual(result[1][0], 3)


if __name__ == '__main__':
    unittest.main()
