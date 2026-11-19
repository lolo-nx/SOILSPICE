"""Hourly alignment, weather recovery and energy-weighted baseline anchoring."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from chapters import annual_irradiance_sensitivity as study


class AnnualIrradianceTests(unittest.TestCase):
    def test_effective_irradiance_and_temperature_recover_existing_yield_inputs(self):
        G=np.array([0.,10.,400.,1000.])
        T=np.array([20.,15.,40.,65.])
        frame=pd.DataFrame({f'Y_{name}':G/1000*(1+gamma*(T-25))
                            for name,gamma in study.energy.GAMMA_T.items()})
        actual_G,actual_T=study.effective_weather(frame)
        np.testing.assert_allclose(actual_G,G,atol=1e-9)
        np.testing.assert_allclose(actual_T[1:],T[1:],atol=1e-9)
        self.assertTrue(np.isnan(actual_T[0]))

    def test_saved_event_intervals_cover_actual_hours_without_inventing_boundaries(self):
        index=pd.date_range('2016-01-01',periods=5,freq='h',tz='UTC')
        hours=pd.DataFrame(index=index)
        series=pd.DataFrame(dict(date=[index[0],index[2]],
            end_time=[index[2],index[-1]+pd.Timedelta(hours=1)],duration_days=[2/24,3/24]))
        result=study.align_saved_series(hours,series)
        np.testing.assert_array_equal(result['interval'],[0,0,1,1,1])
        np.testing.assert_array_equal(result['age_hours'],[1,2,1,2,3])
        broken=series.copy()
        broken.loc[1,'date']=index[3]
        with self.assertRaises(AssertionError):
            study.align_saved_series(hours,broken)

    def test_native_delta_grid_removes_absolute_STC_baseline_offset(self):
        Gs=(100,1000)
        runner=SimpleNamespace(results={})
        for name in study.ORDER:
            for G in Gs:
                P=G*.4
                runner.results[study.electrical.State(name,G=G)]=dict(status='ok',Pmax_W=P)
                for orientation in study.ORIENTATIONS:
                    for tau in study.TAUS:
                        loss=(1-tau)*(.1+.02*(1000-G)/900)
                        runner.results[study.electrical.State(name,'rectangle',orientation,60,tau,G)]=dict(status='ok',Pmax_W=P*(1-loss))
        with tempfile.TemporaryDirectory() as directory:
            runner.out=Path(directory)
            lookup=study.native_grid(runner,Gs)
            for interpolator in lookup.values():
                delta=interpolator([[1000,.5],[100,.5]])
                np.testing.assert_allclose(delta,[0,.01],atol=1e-12)
                # Original archived loss can be 0.30; corrected loss is 0.31,
                # not native absolute 0.06 and not a new temperature baseline.
                np.testing.assert_allclose(.30+delta,[.30,.31])

    def test_annual_changes_are_energy_weighted_not_averaged_by_hour(self):
        hours=pd.DataFrame({f'Y_{name}':[1.,9.] for name in study.ORDER},
                           index=pd.date_range('2016-01-01',periods=2,freq='h',tz='UTC'))
        weather=dict(hours=hours,G=np.array([1000.,100.]),uniform=np.ones(2),
                     tau_ref=np.array([1.,.5]),annual_years=1)
        def interpolator(query):
            return np.where(query[:,0]<1000,.1*(1-query[:,1]),0.)
        lookup={(name,ori):interpolator for name in study.ORDER for ori in study.ORIENTATIONS}
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(study.energy,'slope_length',return_value=1000), \
                patch.object(study.energy,'curve',return_value=np.array([.5,0.])), \
                patch.object(study.energy,'spice_table',return_value=(np.array([0.,1.]),None,None)):
            result,_=study.aggregate_sites({'synthetic':weather},lookup,(1,1000),Path(directory))
            outside=dict(weather,G=np.array([1000.,.1]))
            unsupported,_=study.aggregate_sites({'synthetic':outside},lookup,(1,1000),Path(directory))
        np.testing.assert_allclose(result.original_band_loss_pct,22.5)
        np.testing.assert_allclose(result.corrected_band_loss_pct,27.)
        np.testing.assert_allclose(result.annual_difference_pp,4.5)
        np.testing.assert_allclose(unsupported.corrected_band_loss_pct,22.5)
        np.testing.assert_allclose(unsupported.supported_clean_energy_fraction,.1)
        np.testing.assert_allclose(unsupported.corrected_lower_bound_pct,0)
        np.testing.assert_allclose(unsupported.corrected_upper_bound_pct,90)
        # A simple mean over hours would instead give 2.5 pp: the 9x
        # production weight of the affected hour must be retained.


if __name__=='__main__':
    unittest.main()
