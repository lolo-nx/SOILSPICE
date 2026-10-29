"""SPICE-only loss maps.

For every module, band shape, orientation and target A_cov, the band height is
chosen so that all modules have exactly the same covered fraction, then one
simulation is run per tau_I. The clean reference is simulated once per module.

The output table (one row per state) is what the energy model reads:
    module, shape, orientation, target_A_cov, A_cov, height_mm,
    tau_I, tau_eff, Pmp_loss_%, Vmp_loss_%, Imp_loss_%
"""
import pandas as pd

from . import geometry
from .analysis import losses
from .spice import Simulator, run_tag

CONFIGS = (("rectangle", "portrait"), ("rectangle", "landscape"),
           ("triangle", "portrait"), ("triangle", "landscape"))
ACOV_TARGETS = (0.02, 0.04, 0.06, 0.08, 0.10)
TAU_VALUES = (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)


def loss_map(module, configs=CONFIGS, acov_targets=ACOV_TARGETS, taus=TAU_VALUES, apex="right"):
    """Loss map of one module (DataFrame, one row per simulated state)."""
    sim = Simulator.for_module(module)
    clean = sim.clean().params
    rows = []
    for shape, orientation in configs:
        for target in acov_targets:
            height, achieved = geometry.height_for_acov(module, shape, orientation, target, apex)
            fcov = geometry.fcov_grid(module, shape, orientation, height, apex)
            cells = geometry.cell_fcov(module, fcov, orientation)
            print(f"  {module.name} {shape:<9} {orientation:<9} A_cov={target:.0%} -> h={height:6.1f} mm")
            for tau in taus:
                tag = run_tag(module, shape, orientation, achieved, tau, "matched")
                try:
                    result = sim.run(cells, tau, tag)
                except RuntimeError as err:
                    print(f"    [!] {tag}: {err}")
                    continue
                rows.append({"module": module.name, "shape": shape, "orientation": orientation,
                             "target_A_cov": target, "A_cov": achieved, "height_mm": height,
                             "tau_I": tau, "tau_eff": geometry.tau_eff(achieved, tau),
                             **losses(result.params, clean)})
    return pd.DataFrame(rows)


def loss_maps(modules, **kwargs):
    """Loss maps of several modules, stacked in one table."""
    tables = []
    for module in modules:
        try:
            tables.append(loss_map(module, **kwargs))
        except FileNotFoundError as err:
            print(f"  [!] {module.name} skipped: {err}")
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


# --- equal band height ------------------------------------------------------------------
# The loss maps above compare modules at the same covered fraction of the cells
# (A_cov), which removes the geometry. Here every module gets the same band
# height h, measured from the lower edge of the module: the margin below the
# cells, the cell height and the module size now count, as they do for a real
# soiling band. Heights up to 10 cm cover what is observed for bottom-edge
# soiling on tilted modules.
HEIGHTS_MM = (10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80, 100)


def height_map(module, configs=CONFIGS, heights=HEIGHTS_MM, taus=TAU_VALUES, apex="right"):
    """Losses of one module for a band of fixed height h [mm]. A band that
    stays inside the margin shades no cell: zero loss, no simulation."""
    sim = Simulator.for_module(module)
    clean = sim.clean().params
    zero = {"Pmp_loss_%": 0.0, "Vmp_loss_%": 0.0, "Imp_loss_%": 0.0}
    rows = []
    for shape, orientation in configs:
        for h in heights:
            fcov = geometry.fcov_grid(module, shape, orientation, h, apex)
            a = geometry.acov(fcov)
            print(f"  {module.name} {shape:<9} {orientation:<9} h={h:5.0f} mm -> A_cov={a:6.2%}")
            cells = geometry.cell_fcov(module, fcov, orientation) if a > 0 else None
            for tau in taus:
                base = {"module": module.name, "shape": shape, "orientation": orientation,
                        "height_mm": h, "A_cov": a, "tau_I": tau, "tau_eff": geometry.tau_eff(a, tau)}
                if a <= 0:
                    rows.append({**base, **zero})
                    continue
                tag = run_tag(module, shape, orientation, a, tau, f"h{h:g}")
                try:
                    result = sim.run(cells, tau, tag)
                except RuntimeError as err:
                    print(f"    [!] {tag}: {err}")
                    continue
                rows.append({**base, **losses(result.params, clean)})
    return pd.DataFrame(rows)


def height_maps(modules, **kwargs):
    tables = []
    for module in modules:
        try:
            tables.append(height_map(module, **kwargs))
        except FileNotFoundError as err:
            print(f"  [!] {module.name} skipped: {err}")
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
