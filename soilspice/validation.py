"""Comparison of SPICE with the flash-test campaign.

The band of the experiments (80 mm, see experiments.py) is reproduced in
SPICE for every measured tau_I, and both are compared on
  - the I-V parameters (Isc, Voc, Vmp, Imp, Pmax, FF),
  - the losses relative to the clean flash / clean simulation,
  - the whole curve (nRMSE, nMAE, nMBE).
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import experiments, geometry, measurements
from .analysis import curve_errors, error_summary, losses
from .spice import Simulator, run_tag

PARAMETERS = (("Isc_A", "Isc"), ("Voc_V", "Voc"), ("Vmp_V", "Vmp"),
              ("Imp_A", "Imp"), ("Pmax_W", "Pmp"), ("FF_percent", "FF"))


def check_saved_campaign(table, shape=None, include_clean=True):
    """Reject saved results whose configurations/tau_I differ from experiments.py."""
    keys = ["module", "shape", "orientation", "tau_I"]
    if shape is not None:
        table = table.assign(shape=shape)
    expected = {
        (*key, entry["tau_I"])
        for key, entries in experiments.CAMPAIGN.items()
        if shape is None or key[1] == shape
        for entry in entries
        if include_clean or entry["tau_I"] != 1.0
    }
    actual = set(table[keys].itertuples(index=False, name=None))
    if actual != expected:
        raise ValueError(
            "Saved validation results do not match the configurations/tau_I in experiments.py. "
            "Recompute simulations and figures with: python chapters/1_validation.py --full"
        )


@dataclass
class ConfigResult:
    """Measured and simulated curves of one module / shape / orientation."""
    module: object
    shape: str
    orientation: str
    fcov: np.ndarray
    measured: list = field(default_factory=list)      # measurements.Measurement, clean first
    simulated: dict = field(default_factory=dict)     # tau_I -> spice.IVResult

    @property
    def acov(self):
        return geometry.acov(self.fcov)


def band_fcov(module, shape, orientation):
    return geometry.fcov_grid(module, shape, orientation, experiments.BAND_HEIGHT_MM, experiments.TRIANGLE_APEX)


def load_measurements(module, shape, orientation):
    """Flashes of one configuration, clean first. Shaded flashes are corrected
    to the temperature of the clean flash."""
    out = []
    t_clean = None
    for entry in experiments.flashes(module.name, shape, orientation):   # sorted by tau_I, clean first
        try:
            m = measurements.load(entry, module, t_target=t_clean)
            if entry["tau_I"] == 1.0:
                t_clean = m.params["Tmodule_C"]
            elif t_clean is None:
                print(f"  [!] no clean flash for {module.name} {shape} {orientation}: "
                      f"{entry['path']} not temperature-corrected")
            out.append(m)
        except (OSError, ValueError) as err:
            print(f"  [!] {entry['path']}: {err}")
    return out


def run_configuration(module, shape, orientation, sim=None, with_spice=True):
    """Measurements of one configuration and, optionally, the matching simulations."""
    res = ConfigResult(module, shape, orientation, band_fcov(module, shape, orientation),
                       load_measurements(module, shape, orientation))
    if not with_spice or not res.measured:
        return res
    sim = sim or Simulator.for_module(module)
    cells = geometry.cell_fcov(module, res.fcov, orientation)
    for m in res.measured:
        tag = run_tag(module, shape, orientation, res.acov, m.tau_i, "valid")
        try:
            res.simulated[m.tau_i] = sim.clean() if m.tau_i == 1.0 else sim.run(cells, m.tau_i, tag)
        except RuntimeError as err:
            print(f"  [!] {tag}: {err}")
    return res


def parameter_errors(results):
    """Long table: one row per configuration, tau_I and I-V parameter."""
    rows = []
    for r in results:
        for m in r.measured:
            sim = r.simulated.get(m.tau_i)
            if sim is None:
                continue
            for key, short in PARAMETERS:
                meas, spice = m.params.get(key), sim.params.get(key)
                if meas in (None, 0) or spice is None:
                    continue
                rows.append({"module": r.module.name, "shape": r.shape, "orientation": r.orientation,
                             "tau_I": m.tau_i, "parameter": short, "measured": meas, "spice": spice,
                             "rel_error_%": (spice - meas) / meas * 100.0})
    return pd.DataFrame(rows)


def loss_comparison(results):
    """Measured vs simulated losses, one row per shaded flash."""
    rows = []
    for r in results:
        clean_meas = next((m for m in r.measured if m.tau_i == 1.0), None)
        clean_sim = r.simulated.get(1.0)
        if clean_meas is None:
            continue
        for m in r.measured:
            if m.tau_i == 1.0:
                continue
            row = {"module": r.module.name, "shape": r.shape, "orientation": r.orientation,
                   "label": m.label, "tau_I": m.tau_i}
            row.update({k.replace("_loss", "_loss_exp"): v for k, v in losses(m.params, clean_meas.params).items()})
            sim = r.simulated.get(m.tau_i)
            if sim is not None and clean_sim is not None:
                row.update({k.replace("_loss", "_loss_sim"): v for k, v in losses(sim.params, clean_sim.params).items()})
            rows.append(row)
    return pd.DataFrame(rows)


def whole_curve_errors(results):
    rows = []
    for r in results:
        for m in r.measured:
            sim = r.simulated.get(m.tau_i)
            if sim is None:
                continue
            err = curve_errors(m.v, m.i, sim.v, sim.i, m.params.get("Isc_A"), m.params.get("Pmax_W"))
            if err:
                rows.append({"module": r.module.name, "shape": r.shape, "orientation": r.orientation,
                             "tau_I": m.tau_i, **err})
    return pd.DataFrame(rows)


def summarize(param_errors, tol_pct=5.0):
    """MAE / RMSE / max error per configuration and parameter, plus all parameters together."""
    keys = ["module", "shape", "orientation"]
    by_param = (param_errors.groupby(keys + ["parameter"])["rel_error_%"]
                .apply(error_summary, tol_pct=tol_pct).unstack().reset_index())
    overall = (param_errors.groupby(keys)["rel_error_%"]
               .apply(error_summary, tol_pct=tol_pct).unstack().reset_index())
    overall["parameter"] = "all"
    return pd.concat([by_param, overall], ignore_index=True)


# --- measured losses of all modules on a common x axis ---------------------------------

def measured_losses(modules, shape, orientation, with_spice=False):
    """Measured (and optionally simulated) losses of several modules for the
    same 80 mm band, with the severity indicators that make different module
    geometries comparable (tau_eff, see geometry.band_severity)."""
    rows = []
    for module in modules:
        res = run_configuration(module, shape, orientation, with_spice=with_spice)
        if not any(m.tau_i == 1.0 for m in res.measured):
            print(f"  [!] {module.name} {shape} {orientation}: no clean flash, skipped")
            continue
        severity, group = geometry.band_severity(module, res.fcov, shape, orientation)
        table = loss_comparison([res])
        for _, row in table.iterrows():
            rows.append({
                "module": module.name, "label": row["label"], "tau_I": row["tau_I"], "A_cov": res.acov,
                "f_cov_cell": severity if shape == "rectangle" else np.nan,
                "A_cov_substring": severity if shape == "triangle" else np.nan,
                "bpd_group_idx": group, "tau_eff": geometry.tau_eff(severity, row["tau_I"]),
                **{f"{m}_loss_%": row[f"{m}_loss_exp_%"] for m in ("Pmp", "Vmp", "Imp")},
                **{f"{m}_loss_sim_%": row.get(f"{m}_loss_sim_%", np.nan) for m in ("Pmp", "Vmp", "Imp")},
            })
    return pd.DataFrame(rows)


def power_ratio(module, configs):
    """Measured Pmax / Pmax_clean for each configuration of the campaign."""
    rows = []
    for shape, orientation in configs:
        meas = load_measurements(module, shape, orientation)
        clean = next((m for m in meas if m.tau_i == 1.0), None)
        if clean is None:
            continue
        for m in meas:
            rows.append({"module": module.name, "shape": shape, "orientation": orientation, "tau_I": m.tau_i,
                         "Pmax_exp_W": m.params["Pmax_W"], "Pclean_exp_W": clean.params["Pmax_W"],
                         "power_ratio_%": 100.0 * m.params["Pmax_W"] / clean.params["Pmax_W"]})
    return pd.DataFrame(rows)
