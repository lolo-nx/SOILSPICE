"""Reading the flash-test files (.asc export of the I-V tracer).

The file starts with metadata written as pairs of lines (a header line, then
a value line), followed by the measurement table. The separator and the
decimal mark depend on the export settings, so both are detected.
"""
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import config
from .spice import iv_parameters

_V_KEYS = ("umeas", "ucorr", "voltage", "tension", "u[v]", "v[v]")
_I_KEYS = ("imeas", "icorr", "current", "courant", "i[a]", "i[ma]")
_P_KEYS = ("pcorr", "power", "puissance", "p[w]")


def _read_lines(path):
    for enc in ("cp1252", "latin-1", "utf-8"):
        try:
            return path.read_text(encoding=enc).splitlines()
        except UnicodeDecodeError:
            continue
    return path.read_text(encoding="utf-8", errors="replace").splitlines()


def _separator(lines):
    counts = {d: sum(line.count(d) for line in lines[:80]) for d in (",", ";", "\t")}
    best = max(counts, key=counts.get)
    return best if counts[best] else ","


def _table_start(lines):
    """Index of the header line of the measurement table."""
    for k, line in enumerate(lines):
        low = line.lower()
        if ("time" in low or "temps" in low) and any(s in low for s in _V_KEYS):
            return k
    for k, line in enumerate(lines):
        low = line.lower()
        if any(s in low for s in _V_KEYS) and any(s in low for s in _I_KEYS):
            return k
    return None


def _column(columns, keys, prefer="corr"):
    found = [c for c in columns if any(k in c.lower() for k in keys)]
    preferred = [c for c in found if prefer in c.lower()]
    return (preferred or found or [None])[0]


def read_asc(path):
    """Returns (metadata dict, measurement DataFrame)."""
    path = Path(path)
    lines = _read_lines(path)
    if not lines:
        raise ValueError(f"empty file: {path}")
    sep = _separator(lines)
    decimal_comma = sep in (";", "\t") and any(re.search(r"\d,\d", line) for line in lines[:80])
    start = _table_start(lines)
    if start is None:
        raise ValueError(f"no measurement table found in {path} (separator {sep!r})")

    meta = {}
    k = 0
    while k < start:
        line = lines[k]
        if line.strip().lower() == "observations":
            meta["Observations"] = lines[k + 1].strip() if k + 1 < start else ""
            k += 2
        elif line.strip() and sep in line and k + 1 < start:
            meta.update({h.strip(): v.strip() for h, v in zip(line.split(sep), lines[k + 1].split(sep))
                         if h.strip()})
            k += 2
        else:
            k += 1

    columns = [c.strip() for c in lines[start].split(sep)]
    rows = []
    for line in lines[start + 1:]:
        parts = line.split(sep)
        if len(parts) != len(columns):
            continue
        try:
            rows.append([float(p.strip().replace(",", ".") if decimal_comma else p.strip()) for p in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"no usable measurement row in {path}")
    return meta, pd.DataFrame(rows, columns=columns)


def _meta_float(meta, *keys):
    for key in keys:
        try:
            return float(meta[key])
        except (KeyError, ValueError):
            continue
    return None


def iv_from_table(df):
    """Voltage [V] and current [A] columns, first quadrant only."""
    v_col, i_col = _column(df.columns, _V_KEYS), _column(df.columns, _I_KEYS)
    if v_col is None or i_col is None:
        raise ValueError(f"voltage/current columns not found in {list(df.columns)}")
    v = df[v_col].to_numpy(dtype=float)
    i = df[i_col].to_numpy(dtype=float)
    unit = i_col.lower().replace("[", "").replace("]", "")
    if "ma" in unit:
        i = i / 1000.0
    elif not ("[a]" in i_col.lower() or i_col.lower().strip().endswith((" a", "(a)"))):
        i = i / 1000.0 if np.abs(i).max() > 100 else i       # unit not stated: guess from magnitude
    keep = (v >= 0) & (i >= 0)
    return v[keep], i[keep]


def params_from_meta(meta):
    """I-V parameters as computed by the tracer software."""
    isc_ma, imp_ma = _meta_float(meta, "Module current Isc [mA]", "Isc [mA]"), _meta_float(meta, "Imp [mA]")
    g = _meta_float(meta, "Mean irradiance [kW/m2]", "Mean irradiance [kW/m²]")
    return {
        "Voc_V": _meta_float(meta, "Module voltage Voc [V]", "Voc [V]"),
        "Isc_A": isc_ma / 1000.0 if isc_ma is not None else None,
        "Vmp_V": _meta_float(meta, "Vmp [V]"),
        "Imp_A": imp_ma / 1000.0 if imp_ma is not None else None,
        "Pmax_W": _meta_float(meta, "Pmax [W]"),
        "FF_percent": _meta_float(meta, "Fill factor [%]"),
        "Irradiance_Wm2": g * 1000.0 if g is not None else None,
        "Tmodule_C": _meta_float(meta, "Module temperature [C]", "Module temperature [°C]"),
    }


def correct_temperature(v, i, t_meas, t_ref, kv, kc):
    """Linear translation of the curve to t_ref with relative coefficients."""
    dt = t_ref - t_meas
    return v * (1.0 + kv * dt), i * (1.0 + kc * dt)


@dataclass
class Measurement:
    label: str
    tau_i: float
    v: np.ndarray
    i: np.ndarray
    params: dict
    path: Path


def load(entry, module, t_target=None):
    """One flash from an experiments.py entry.

    The clean flash (tau_I = 1) is kept as measured. Shaded flashes are moved
    to t_target, the module temperature of the clean flash as written by the
    flasher, so that every flash of a configuration is at the same temperature
    as the clean one. Their parameters are recomputed from the corrected curve.

    module.t_ref is not used here: it is the temperature of the SPICE model
    (fitted on the clean flash), not a measured temperature."""
    path = Path(entry["path"])
    meta, table = read_asc(path)
    v, i = iv_from_table(table)
    params = params_from_meta(meta)
    t_meas = params["Tmodule_C"]
    if t_meas is None:
        raise ValueError(f"no module temperature in {path}")

    if config.NORMALIZE_TO_REF_TEMPERATURE and entry["tau_I"] != 1.0 and t_target is not None:
        v, i = correct_temperature(v, i, t_meas, t_target, module.kv, module.kc)
        order = np.argsort(v)
        v, i = v[order], i[order]
        params.update(iv_parameters(v, i))
        params["Tmodule_C"] = t_target
    else:
        order = np.argsort(v)
        v, i = v[order], i[order]
    return Measurement(entry.get("label") or path.stem, entry["tau_I"], v, i, params, path)
