"""Physical limits and cell-by-cell netlist isolation of the Yingli shunt law."""
from dataclasses import replace
from pathlib import Path
from math import exp
import re
import unittest

import numpy as np

from soilspice import modules, spice


class YingliShuntTests(unittest.TestCase):
    def setUp(self):
        self.module = replace(modules.get("YINGLI"), rp_law=modules.PVsystShuntLaw(r_stc=5.99))

    def test_endpoints_and_partial_coverage(self):
        self.assertAlmostEqual(spice.shaded_rp(self.module, 0., 0.), 5.99)
        self.assertAlmostEqual(spice.shaded_rp(self.module, 1., 1.), 5.99)
        self.assertAlmostEqual(spice.shaded_rp(self.module, 1., 0.), 23.96)
        # Same photocurrent factor g=.5 must give the same resistance, including an opaque partial band.
        self.assertAlmostEqual(spice.shaded_rp(self.module, .5, 0.),
                               spice.shaded_rp(self.module, 1., .5))
        base = (5.99 - 23.96 * exp(-5.5)) / (1 - exp(-5.5))
        self.assertAlmostEqual(spice.shaded_rp(self.module, .7, .1),
                               base + (23.96 - base) * exp(-5.5 * .37))
        values = [spice.shaded_rp(self.module, 1., tau) for tau in np.linspace(0, 1, 101)]
        self.assertTrue(np.all(np.diff(values) < 0))

    def test_all_modules_use_constant_shunts_by_default(self):
        for name, module in modules.MODULES.items():
            self.assertIsNone(module.rp_law, name)

    def test_netlist_changes_only_cell_shunt_parameters(self):
        template = spice.read_netlist(self.module.netlist_path)
        coverage = {name: (.314159 if k % 2 else 0.)
                    for k, name in enumerate(self.module.instances.values())}
        variable = spice.patch_cells(template, self.module, coverage, .442)
        constant = spice.patch_cells(template, replace(self.module, rp_law=None), coverage, .442)
        stripped = "\n".join(re.sub(r" rp=[^\s]+", "", line)
                             if line.startswith(self.module.cell_prefix) else line
                             for line in variable.splitlines()) + "\n"
        self.assertEqual(stripped, constant)
        for line in variable.splitlines():
            if not line.startswith(self.module.cell_prefix):
                continue
            f = float(re.search(r"f_cov=(\S+)", line)[1])
            rp = float(re.search(r"rp=(\S+)", line)[1])
            self.assertAlmostEqual(rp, spice.shaded_rp(self.module, f, .442), places=9)
        self.assertNotIn("pwr(", variable.lower())

    def test_simulator_cache_separates_shunt_models(self):
        from tempfile import TemporaryDirectory
        from unittest.mock import patch
        with TemporaryDirectory() as temp, patch.object(spice.config, "RUNS_DIR", Path(temp)):
            old = dict(spice._SIMULATORS)
            try:
                spice._SIMULATORS.clear()
                variable = spice.Simulator.for_module(self.module)
                constant = spice.Simulator.for_module(replace(self.module, rp_law=None))
                self.assertIsNot(variable, constant)
                self.assertIs(variable, spice.Simulator.for_module(self.module))
            finally:
                spice._SIMULATORS.clear()
                spice._SIMULATORS.update(old)


if __name__ == "__main__":
    unittest.main()
