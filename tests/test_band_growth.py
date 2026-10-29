"""Event-total washing, mass balance, and consistency of energy calculations."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from soilspice import energy


class BandGrowthTests(unittest.TestCase):
    def test_event_total_threshold_resets_both_deposits_after_transfer(self):
        rain = np.array([0.999, 1.0, 14.9, 15.0, 0.0])
        _, _, _, states = energy.simulate(
            np.ones(5), rain, 1.0, np.array([2.0]), np.zeros((1, len(energy.TAU))), track=0)
        tau, _, uniform_loss = np.array(states).T
        np.testing.assert_allclose(tau, np.exp(-energy.ALPHA * np.array([0, 0, 4, 6, 0])))
        np.testing.assert_allclose(uniform_loss, 1 - np.exp(-energy.ALPHA * np.array([1, 2, 1, 1, 1])))

    def test_reference_and_sensitivity_thresholds_use_event_totals(self):
        self.assertEqual(energy.REFERENCE_TOTAL, 15)
        self.assertEqual(energy.TOTAL_THRESHOLDS, (5, 10, 15, 20, 25))
        totals = np.array([0.999, 1, 4.999, 5, 9.999, 10, 14.999, 15, 19.999, 20, 24.999, 25])
        for threshold in energy.TOTAL_THRESHOLDS:
            np.testing.assert_array_equal(energy.washing_flags(totals, washing_total=threshold), totals >= threshold)
        np.testing.assert_array_equal(energy.washing_flags(totals), totals >= 15)
        for threshold in (0.0, -1.0, np.inf, np.nan):
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                energy.washing_flags([1.0], washing_total=threshold)
        for total in (np.nan, np.inf, -1):
            with self.subTest(total=total), self.assertRaises(ValueError):
                energy.washing_flags([total])

    def test_series_classification_uses_event_total_and_resets_both_deposits(self):
        days = self.weather()
        with patch.object(energy, 'SITES', ('Neuchâtel',)), \
                patch.object(energy, 'site_steps', return_value=days), \
                patch.object(energy, 'attenuated_yields', side_effect=lambda site, D, alpha: days[['Y_JASOLAR']] * np.exp(-alpha * D)), \
                patch.object(energy, 'lut', return_value=np.zeros((1, len(energy.TAU)))):
            series = energy.growth_series()
        self.assertEqual(series.rain_regime.tolist(),
                         ['negligible', 'transfer', 'washing', 'transfer', 'transfer', 'no_event_end'])
        self.assertLess(series.tau.iloc[2], 1)
        self.assertEqual(series.tau.iloc[3], 1)
        self.assertAlmostEqual(series.L_u.iloc[3], series.L_u.iloc[0])

    @staticmethod
    def weather():
        days = pd.DataFrame({
            'rain': [0.999, 1.0, 30.0, 5.0, 1.0, 0.0],
            'rain_max': [0.5, 1.0, 0.5, 5.0, 1.0, 0.0],
            'Y_JASOLAR': [0.5, 2.0, 1.0, 0.1, 3.0, 1.0],
        }, index=pd.date_range('2020-01-01', periods=6, tz='UTC'))
        days['duration_days'] = 1.0
        days['clean_start'] = True
        days['end_time'] = days.index + pd.Timedelta(days=1)
        days['event_id'] = [0, 1, 2, 3, 4, -1]
        return days

    def test_washing_mask_must_match_intervals(self):
        with self.assertRaisesRegex(ValueError, 'one flag per simulated interval'):
            energy.simulate(np.ones(2), np.zeros(2), 1.0, np.ones(1),
                            np.zeros((1, len(energy.TAU))), washing_days=[True])

    def test_maps_cost_runs_and_series_agree_with_energy_weighting(self):
        days = self.weather()
        loss_table = (0.2 * (1 - energy.TAU))[None, :]
        with patch.object(energy, 'SITES', ('Neuchâtel',)), \
                patch.object(energy, 'site_steps', return_value=days), \
                patch.object(energy, 'attenuated_yields', side_effect=lambda site, D, alpha: days[['Y_JASOLAR']] * np.exp(-alpha * D)), \
                patch.object(energy, 'ORDER', ('JASOLAR',)), \
                patch.object(energy, 'ORIENTATIONS', ('landscape',)), \
                patch.object(energy, 'F3', [0.1]), \
                patch.object(energy, 'MARGINS', [0.0]), \
                patch.object(energy, 'INTERVALS', []), \
                patch.object(energy, 'MARKETS', {'ground': (0.06, 0.41)}), \
                patch.object(energy, 'spice_table', return_value=(np.array([0, 0.9]), np.array([60.0]), {})), \
                patch.object(energy, 'lut', return_value=loss_table):
            series = energy.growth_series(washing_total=15)
            grid = energy.tau_grid(fs=[0.1], washing_total=15)
            runs = energy.model_runs(washing_total=15)
        loss = np.average(1 - (1 - series.L_u) * (1 - series.L_band), weights=series.Y)
        tau = np.average(series.tau, weights=series.Y)
        self.assertAlmostEqual(grid.loss.iloc[0], loss)
        self.assertAlmostEqual(grid.tau_mean.iloc[0], tau)
        reference = runs[runs['mode'] == 'none'].iloc[0]
        self.assertAlmostEqual(reference.loss, loss)
        self.assertAlmostEqual(reference.tau_mean, tau)

    def test_cleaning_bounds_and_cycles_use_event_totals(self):
        days = pd.DataFrame({'rain': [1, 5, 30, 0, 1, 0], 'rain_max': [1, 0.5, 5, 0, 1, 0], 'Y': np.ones(6)},
                            index=pd.date_range('2020-01-01', periods=6, tz='UTC'))
        days['clean_start'] = True
        with patch.object(energy, 'SITES', ('Neuchâtel',)), patch.object(energy, 'site_steps', return_value=days):
            bound = energy.cleaning_bound(washing_total=15)
            cycles = energy.washing_cycles(washing_total=15)
        np.testing.assert_allclose(bound.washing_per_yr, 1.0)
        np.testing.assert_allclose(bound.E_next, 2.0)
        self.assertEqual(cycles['Neuchâtel'], 6.0)


if __name__ == '__main__':
    unittest.main()
