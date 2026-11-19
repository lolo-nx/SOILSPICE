"""Independent discounted cash-flow and saved-input accounting checks."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from soilspice.robot_economics import (annuity_factor, robot_cashflows,
    robot_ownership_cost, robot_strategy_costs, ROBOT_INVESTMENT_EUR_KWP)


class RobotEconomicsTests(unittest.TestCase):
    def test_purchase_timing_and_yearly_maintenance(self):
        flows = robot_cashflows(100, residual_mode="none")
        self.assertEqual(flows.loc[flows.investment_eur_kwp > 0, "year"].tolist(), [0, 12, 24])
        self.assertEqual(flows.maintenance_eur_kwp.iloc[0], 0)
        np.testing.assert_array_equal(flows.maintenance_eur_kwp.iloc[1:], np.full(25, 5.))
        expected = sum(100 / 1.05**y for y in (0, 12, 24)) + sum(5 / 1.05**y for y in range(1, 26))
        self.assertAlmostEqual(flows.discounted_net_cost_eur_kwp.sum(), expected, places=12)
        self.assertAlmostEqual(robot_ownership_cost(100, residual_mode="none")["ownership_eac_eur_kwp_year"] * annuity_factor(), expected)

    def test_residual_is_explicit_discounted_and_not_an_avoided_replacement(self):
        no = robot_cashflows(120, residual_mode="none")
        yes = robot_cashflows(120, residual_mode="linear_remaining_life")
        self.assertEqual(yes.residual_credit_eur_kwp.iloc[-1], 110.)
        self.assertEqual(yes.investment_eur_kwp.iloc[24], 120.)
        self.assertAlmostEqual(no.discounted_net_cost_eur_kwp.sum()-yes.discounted_net_cost_eur_kwp.sum(), 110/1.05**25)
        with self.assertRaises(TypeError):
            robot_cashflows(120)
        with self.assertRaises(ValueError):
            robot_cashflows(120, residual_mode="resale_guaranteed")

    def test_zero_rate_and_exact_lifetime_boundary(self):
        self.assertEqual(annuity_factor(25, 0), 25)
        flows = robot_cashflows(120, residual_mode="linear_remaining_life", years=24, rate=0)
        self.assertEqual(flows.investment_eur_kwp.sum(), 240)
        self.assertEqual(flows.residual_credit_eur_kwp.sum(), 0)
        self.assertEqual(flows.net_cost_eur_kwp.sum(), 240 + 24*6)

    def test_ratio_uses_each_modules_clean_revenue_and_adds_cost_once(self):
        rows = pd.DataFrame(dict(mode=["robot"]*3, revenue=[50., 100., 200.],
                                 loss=[.01]*3, loss_u=[.005]*3,
                                 cost=[.5, 1., 2.], cost_u=[.25, .5, 1.]))
        before=rows.copy(deep=True)
        result=robot_strategy_costs(rows, residual_mode="none")
        pd.testing.assert_frame_equal(rows, before)
        for bound, investment in zip(("low", "mid", "high"),
                (ROBOT_INVESTMENT_EUR_KWP[0], np.mean(ROBOT_INVESTMENT_EUR_KWP), ROBOT_INVESTMENT_EUR_KWP[1])):
            annual=robot_ownership_cost(investment, residual_mode="none")["ownership_eac_eur_kwp_year"]
            np.testing.assert_allclose(result[f"band_{bound}_pct_clean_revenue"], 1+100*annual/rows.revenue)
            np.testing.assert_allclose(result[f"band_{bound}_pct_clean_revenue"]-result[f"noband_{bound}_pct_clean_revenue"], .5)
            np.testing.assert_allclose(result[f"band_total_{bound}_pv_eur_kwp"]/(rows.revenue*annuity_factor())*100,
                                       result[f"band_{bound}_pct_clean_revenue"])
        bad=rows.copy()
        bad.loc[0,'cost']+=1
        with self.assertRaises(AssertionError):
            robot_strategy_costs(bad, residual_mode="none")

    def test_figure_bars_error_bounds_and_noband_markers_use_discounted_totals(self):
        # Integration with actual saved physical results; never rerun physics.
        from soilspice import config, fig_energy
        import matplotlib.pyplot as plt
        from matplotlib.container import BarContainer
        scenario=config.CH4_DIR/'total_scenarios/PH15mm'
        if not (scenario/'data/model_runs.csv').exists():
            self.skipTest('Saved PH15mm results not available')
        runs=pd.read_csv(scenario/'data/model_runs.csv')
        best=pd.read_csv(scenario/'data/best_strategy.csv')
        figures=[]
        for mode in ('none','linear_remaining_life'):
            with patch.object(fig_energy,'save',side_effect=lambda fig,path: figures.append(fig)):
                fig_energy.cleaning_figures(runs,best,scenario,robot_residual_mode=mode)
        self.addCleanup(lambda: [plt.close(fig) for fig in figures])
        for fig,mode in zip(figures,('none','linear_remaining_life')):
            axis=fig.axes[2]
            bars=[container for container in axis.containers if isinstance(container,BarContainer)]
            self.assertEqual(len(bars),12)
            for i,site in enumerate(fig_energy.SITE):
                rows=runs[(runs.site==site)&(runs.market=='ground')&(runs['mode']=='robot')
                          &np.isclose(runs.h,60)&np.isclose(runs.f,.1)&np.isclose(runs.margin,0)]
                factor=(sum(1/1.05**year for year in (0,12,24))/annuity_factor()+.05)
                if mode!='none':
                    factor-=(11/12)/1.05**25/annuity_factor()
                lo,hi=ROBOT_INVESTMENT_EUR_KWP
                expected=[np.median(100*(rows.cost+C*factor)/rows.revenue) for C in (lo,(lo+hi)/2,hi)]
                robot=bars[i*4+3]
                self.assertAlmostEqual(robot.patches[0].get_height(),expected[1])
                segment=robot.errorbar.lines[2][0].get_segments()[0]
                np.testing.assert_allclose(segment[:,1],[expected[0],expected[2]])
                x=3+(i-1)*.26
                markers=[line for line in axis.lines if len(line.get_xdata())==2
                         and np.allclose(line.get_xdata(),[x-.26*.4,x+.26*.4])]
                self.assertEqual(len(markers),1)
                expected_noband=np.median(100*(rows.cost_u+(lo+hi)/2*factor)/rows.revenue)
                np.testing.assert_allclose(markers[0].get_ydata(),expected_noband)
        # Rain/manual cost curves and trigger frequencies are scenario-independent.
        for index in (0,1):
            a,b=figures[0].axes[index],figures[1].axes[index]
            self.assertEqual(len(a.lines),len(b.lines))
            for left,right in zip(a.lines,b.lines):
                np.testing.assert_array_equal(left.get_xdata(),right.get_xdata())
                np.testing.assert_array_equal(left.get_ydata(),right.get_ydata())


if __name__ == "__main__":
    unittest.main()
