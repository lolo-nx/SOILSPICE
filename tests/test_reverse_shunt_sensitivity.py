"""Checks for the isolated continuous reverse-only diagnostic model."""
import unittest

import numpy as np

from chapters.yingli_reverse_shunt_sensitivity import RESISTANCES, patch_shunt, shunt_current
from soilspice import modules
from soilspice.spice import read_netlist


class ReverseShuntSensitivityTests(unittest.TestCase):
    def test_forward_branch_and_zero_continuity(self):
        for resistance in RESISTANCES:
            np.testing.assert_allclose(shunt_current([0., .1, .6], resistance), np.array([0., .1, .6]) / 5.99)
            self.assertEqual(float(shunt_current(0., resistance)), 0.)
            self.assertLess(abs(float(shunt_current(-1e-9, resistance))), 1e-9)
            self.assertLess(abs(float(shunt_current(1e-9, resistance))), 1e-9)

    def test_reverse_branch(self):
        for resistance in RESISTANCES:
            self.assertAlmostEqual(float(shunt_current(-5., resistance)), -5. / resistance)
        self.assertAlmostEqual(float(shunt_current(-5., 1e6)), -5e-6)

    def test_netlist_only_replaces_the_shunt_element(self):
        text = read_netlist(modules.get("YINGLI").netlist_path)
        for resistance in RESISTANCES:
            changed = patch_shunt(text, resistance)
            a, b = text.splitlines(), changed.splitlines()
            differences = [(before, after) for before, after in zip(a, b) if before != after]
            self.assertEqual(len(a), len(b))
            self.assertEqual(len(differences), 1)
            self.assertEqual(differences[0][0], "Rp N001 - {rp}")
            self.assertTrue(differences[0][1].startswith("Bsh N001 - I=if(V(N001,-)>=0,V(N001,-)/{rp},"))
            self.assertNotIn("pwr(", changed.lower())


if __name__ == "__main__":
    unittest.main()
