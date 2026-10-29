"""Inter-washing intervals must span years and distinguish censored boundary waits."""
import unittest

import numpy as np
import pandas as pd

from soilspice.rain_scenarios import event_intervals


class ScenarioIntervalTests(unittest.TestCase):
    def setUp(self):
        self.start = pd.Timestamp('2020-01-01', tz='UTC')
        self.stop = pd.Timestamp('2022-01-01', tz='UTC')
        self.events = pd.DataFrame({
            'end': pd.to_datetime(['2020-12-20', '2021-01-10', '2021-02-10', None], utc=True),
            'rain_max': [3.0, 3.0, 0.5, 50.0], 'rain': [9.0, 12.0, 30.0, 100.0],
            'complete': [True, True, True, False],
        })

    def test_gap_spans_new_year_and_incomplete_event_is_not_a_washing(self):
        metrics, intervals = event_intervals(self.events, self.start, self.stop, 5.0)
        self.assertEqual(metrics['washing_events'], 3)
        self.assertEqual(metrics['max_between_washings_days'], 31)
        self.assertGreater(metrics['max_observed_without_washing_days'], 300)
        self.assertEqual(intervals.left_censored.sum(), 1)
        self.assertEqual(intervals.right_censored.sum(), 1)
        self.assertEqual(intervals.days.sum(), 731)

    def test_no_washing_has_no_complete_gap_and_whole_period_is_censored(self):
        metrics, intervals = event_intervals(self.events, self.start, self.stop, 50.0)
        self.assertEqual(metrics['washing_events'], 0)
        self.assertTrue(np.isnan(metrics['max_between_washings_days']))
        self.assertEqual(metrics['max_observed_without_washing_days'], 731)
        self.assertTrue(intervals.left_censored.iloc[0] and intervals.right_censored.iloc[0])

    def test_total_rule_selects_only_total_without_peak_override(self):
        metrics, _ = event_intervals(self.events, self.start, self.stop, 10)
        self.assertEqual(metrics['washing_events'], 2)
        self.assertEqual(metrics['max_between_washings_days'], 31)


if __name__ == '__main__':
    unittest.main()
