"""Counterfactual optimisation, signed profits and threshold/ranking reversals."""
import unittest

import numpy as np
import pandas as pd

from soilspice import energy
from soilspice import energy_sensitivity as S


def saved_runs():
    rows = []
    for margin, costs in ((0, (20, 6, 9)), (10, (18, 5, 8))):
        for mode, k, cost, cost_u, frequency, frequency_u in (
            ('none', 0, costs[0], 10, 0, 0),
            ('trigger', 0, costs[1], 5, 2, 2),
            ('interval', 30, costs[2], 2, 4, 1),
        ):
            rows.append(dict(site='Neuchâtel', market='ground', module='JASOLAR', orientation='portrait',
                             h=60., f=.1, margin=float(margin), mode=mode, k=k, Y=1000., revenue=100.,
                             loss=.2 if mode == 'none' else .03, loss_u=.1 if mode == 'none' else .01,
                             cost=float(cost), cost_u=float(cost_u), cleans=float(frequency), cleans_u=float(frequency_u)))
    return pd.DataFrame(rows)


class EnergySensitivityTests(unittest.TestCase):
    def test_independent_optima_have_different_policies_and_frequencies(self):
        q = S.summarize_runs(saved_runs(), 10, slope_length=lambda m, o: 1000).iloc[0]
        self.assertEqual(q.opt_total_policy, 'trigger')
        self.assertEqual(q.opt_uniform_policy, 'interval')
        self.assertEqual(q.opt_total_cleanings_per_year, 2)
        self.assertEqual(q.opt_uniform_cleanings_per_year, 1)
        self.assertEqual(q.cost_opt_band_extra_eur_kwp_year, 4)  # 6 - 2; not 6 - 5.
        self.assertAlmostEqual(q.energy_none_band_extra_loss_pct, 10)
        self.assertEqual(q.cost_none_band_extra_eur_kwp_year, 10)

    def test_margin_deducts_investment_after_independent_reoptimisation(self):
        q = S.summarize_runs(saved_runs(), 10, slope_length=lambda m, o: 1000).iloc[0]
        self.assertEqual(q.margin10_capex_eur_kwp, 1)
        self.assertEqual(q.margin10_savings_eur_kwp_year, 1)
        self.assertAlmostEqual(q.margin10_npv25_eur_kwp, energy.AF - 1)
        self.assertAlmostEqual(q.margin10_net_eur_kwp_year * energy.AF, q.margin10_npv25_eur_kwp)
        self.assertEqual(q.margin10_state, 'positive')

    def test_unprofitable_manual_option_is_not_hidden_by_including_none(self):
        frame = saved_runs()
        frame.loc[frame['mode'] == 'none', 'cost'] = 2.
        frame.loc[frame['mode'] == 'trigger', 'cost'] = 3.
        frame.loc[frame['mode'] == 'interval', 'cost'] = 4.
        q = S.summarize_runs(frame, 10, slope_length=lambda m, o: 1000).iloc[0]
        self.assertEqual(q.manual_gain_eur_kwp_year, -1)
        self.assertEqual(q.manual_gain_state, 'negative')
        self.assertEqual(q.manual_total_policy, 'trigger')
        self.assertEqual(q.selected_gain_eur_kwp_year, 0)
        self.assertEqual(q.opt_total_policy, 'none')
        self.assertEqual(q.opt_total_cleanings_per_year, 0)
        self.assertEqual(q.margin10_state, 'negative')

    def test_tied_policies_export_frequency_range(self):
        frame = saved_runs()
        frame.loc[(frame.margin == 0) & (frame['mode'] == 'interval'), 'cost'] = 6.
        q = S.summarize_runs(frame, 10, slope_length=lambda m, o: 1000).iloc[0]
        self.assertEqual(q.opt_total_policy, 'trigger')
        self.assertEqual(q.opt_total_cooptimal_policy_count, 2)
        self.assertEqual(q.opt_total_cooptimal_frequency_min, 2)
        self.assertEqual(q.opt_total_cooptimal_frequency_max, 4)

    def test_robustness_retains_individual_signs_ties_and_all_adjacent_reversals(self):
        rows = []
        for j, t in enumerate(S.THRESHOLDS):
            for o in ('portrait', 'landscape'):
                rows.append(dict(site='Neuchâtel', market='ground', module='JASOLAR', orientation=o,
                                 h_mm=60., f=.1, R_H_mm_event=t,
                                 margin10_net_eur_kwp_year=[-1, 0, 1, -1, 1][j] if o == 'portrait' else 1,
                                 manual_gain_eur_kwp_year=[-1, 0, 0, 1, 1][j] if o == 'portrait' else 1,
                                 energy_none_band_extra_loss_pct=[2, 2, 1, 0, 2][j] if o == 'portrait' else 1,
                                 cost_opt_band_extra_pct_clean_revenue=[2, 1, 1, 0, 2][j] if o == 'portrait' else 1))
        table = pd.DataFrame(rows)
        orientations = S.orientation_comparison(table)
        robust, rank, changes = S.classify_across_thresholds(table, orientations)
        p = robust[robust.orientation == 'portrait'].iloc[0]
        l = robust[robust.orientation == 'landscape'].iloc[0]
        self.assertFalse(p.margin10_all_profitable)
        self.assertTrue(l.margin10_all_profitable)
        self.assertEqual(p.margin10_states_5_10_15_20_25, 'negative;zero;positive;negative;positive')
        self.assertFalse(rank.energy_ranking_stable.iloc[0])
        energy_changes = changes[changes.conclusion == 'orientation_energy']
        self.assertEqual(list(zip(energy_changes.R_H_low_mm, energy_changes.R_H_high_mm)), [(10, 15), (15, 20), (20, 25)])
        self.assertEqual(energy_changes.state_low.tolist(), ['landscape', 'tie', 'portrait'])
        with self.assertRaisesRegex(ValueError, 'all five thresholds'):
            S.classify_across_thresholds(table[table.R_H_mm_event != 15], orientations)

    def test_complete_washing_gap_crosses_year_and_excludes_long_censored_wait(self):
        dates = pd.to_datetime(['2016-01-01', '2016-12-20', '2017-01-10', '2026-01-01'], utc=True)
        intervals = pd.DataFrame(dict(site='Neuchâtel', P_H_mm=10,
                                      start=dates[:-1], end=dates[1:],
                                      days=np.diff(dates.asi8) / pd.Timedelta(days=1).value,
                                      left_censored=[True, False, False], right_censored=[False, False, True]))
        q, _ = S.washing_summary(intervals)
        self.assertEqual(q.max_complete_washing_gap_days.iloc[0], 21)
        self.assertEqual(q.washing_events.iloc[0], 2)
        self.assertEqual(q.washing_events_per_calendar_year.iloc[0], .2)
        self.assertGreater(q.final_censored_days.iloc[0], 3000)
        self.assertEqual(q.initial_censored_start.iloc[0], dates[0])

    def test_missing_final_source_hour_is_not_added_to_washing_intervals(self):
        dates = pd.to_datetime(['2016-01-01', '2020-06-01', '2025-12-31 23:00'], utc=True, format='mixed')
        intervals = pd.DataFrame(dict(site='Neuchâtel', P_H_mm=15,
                                      start=dates[:-1], end=dates[1:],
                                      days=np.diff(dates.asi8) / pd.Timedelta(days=1).value,
                                      left_censored=[True, False], right_censored=[False, True]))
        q, _ = S.washing_summary(intervals)
        self.assertEqual(q.washing_events.iloc[0], 1)
        self.assertEqual(q.missing_trailing_source_hours.iloc[0], 1)
        self.assertEqual(q.study_stop_utc.iloc[0], dates[-1])
        self.assertAlmostEqual(q.economic_years_365p25.iloc[0], (3653 - 1/24) / 365.25)


if __name__ == '__main__':
    unittest.main()
