"""Controls/normalisation and global-MPP checks for the isolated sensitivity."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from chapters import stc_sensitivity as study


class STCSensitivityTests(unittest.TestCase):
    def test_requested_configuration_grid_is_complete_and_exact(self):
        grid=study.cases()
        self.assertEqual(len(grid),56)
        self.assertEqual(len(set(grid)),56)
        self.assertEqual(sum(row[3]==80 for row in grid),48)
        self.assertEqual(sum(row[3]==60 for row in grid),8)
        self.assertEqual({row[4] for row in grid}, {.71,.44,0})

    def test_explicit_operating_temperature_and_irradiance_do_not_refit_cells(self):
        for name in study.ORDER:
            with self.subTest(module=name):
                baseline=study.spice.read_netlist(study.modules.get(name).netlist_path)
                text,diodes,nodes=study.netlist(study.State(name,'rectangle','portrait',80,.44,400),Path('standard.dio'))
                self.assertIn('.temp 25',text)
                self.assertRegex(text,r'(?im)^\.param\s+Vi\s*=\s*400\b')
                self.assertRegex(text,r'(?im)^\.param\s+t\s*=\s*25\b')
                self.assertEqual(len(diodes),3)
                # Exact fitted coefficients: only operating controls t/Vi may differ.
                import re
                for key in ('IS','IS1','n1','rp','rs','Vb','a','b','kc'):
                    regex=rf'\b{key}\s*=\s*([-+\d.eE]+)'
                    original=re.findall(regex,baseline,re.I)
                    actual=re.findall(regex,text,re.I)
                    self.assertEqual(original,actual,key)
                legacy,_,_=study.netlist(study.State(name,track='legacy'),Path('standard.dio'))
                self.assertNotIn('.temp ',legacy)

    def test_global_mpp_includes_competing_high_voltage_peak(self):
        voltage=np.linspace(0,5,51)
        current=np.interp(voltage,np.arange(6),[10,8,3,4,2,-1])
        raw={'v1':voltage,'i(v1)':current,'v(out)':voltage}
        for k in (1,2,3):
            raw[f'i(d{k})']=np.where(np.isclose(voltage,3),.5 if k==2 else 0,0)
            raw[f'v(d{k})']=np.full(len(voltage),.5)
        with patch.object(study,'read_raw',return_value=raw):
            result=study.measure(Path('unused.raw'),[(f'd{k}','0') for k in (1,2,3)],('out','0'))
        self.assertEqual(result['Pmax_W'],12)
        self.assertEqual(result['Vmp_V'],3)
        self.assertEqual(result['bypass_MPP_mask'],2)
        self.assertTrue(result['Voc_bracketed'])

    def test_missing_voc_allowed_only_for_archived_reproduction(self):
        voltage=np.linspace(0,4,41)
        raw={'v1':voltage,'i(v1)':np.interp(voltage,np.arange(5),[10,8,3,1,.1]),'v(out)':voltage}
        for k in (1,2,3):
            raw[f'i(d{k})']=np.zeros(len(voltage))
            raw[f'v(d{k})']=np.zeros(len(voltage))
        with patch.object(study,'read_raw',return_value=raw):
            with self.assertRaises(ValueError):
                study.measure(Path('unused.raw'),[(f'd{k}','0') for k in (1,2,3)],('out','0'))
            result=study.measure(Path('unused.raw'),[(f'd{k}','0') for k in (1,2,3)],('out','0'),require_voc_bracket=False)
            self.assertFalse(result['Voc_bracketed'])

    def test_clean_and_soiled_are_normalised_at_same_irradiance(self):
        class MockRunner:
            results={}
        runner=MockRunner()
        for name in study.ORDER:
            for G in study.GS:
                runner.results[study.State(name,G=G)]=dict(status='ok',Pmax_W=.4*G,bypass_MPP_mask=0)
        for name,shape,ori,h,tau in study.cases():
            for G in study.GS:
                runner.results[study.State(name,shape,ori,h,tau,G)]=dict(status='ok',Pmax_W=.4*G*.8,bypass_MPP_mask=0)
        with tempfile.TemporaryDirectory() as temp:
            runner.out=Path(temp)
            frame,summary=study.compile_results(runner)
        self.assertEqual(len(frame),280)
        np.testing.assert_allclose(frame.relative_power_loss_percent,20)
        np.testing.assert_allclose(frame.Delta_L_pp,0,atol=1e-12)
        np.testing.assert_array_equal(summary.valid_conditions,5)


if __name__=='__main__':
    unittest.main()
