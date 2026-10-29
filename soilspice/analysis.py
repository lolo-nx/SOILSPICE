"""Losses, error metrics and resilience indicators.

Sign convention used everywhere: loss_% = (X - X_clean) / X_clean * 100,
negative for a loss, positive for a gain.
"""
import numpy as np
import pandas as pd

METRICS = (("Pmax_W", "Pmp"), ("Vmp_V", "Vmp"), ("Imp_A", "Imp"))


def relative_change(value, clean):
    if value is None or clean in (None, 0):
        return np.nan
    return (value - clean) / clean * 100.0


def losses(params, clean):
    """{'Pmp_loss_%': ..., 'Vmp_loss_%': ..., 'Imp_loss_%': ...}"""
    return {f"{short}_loss_%": relative_change(params.get(key), clean.get(key)) for key, short in METRICS}


# --- comparing a simulated curve with a measured one -------------------------------

def curve_errors(v_meas, i_meas, v_sim, i_sim, isc_ref, pmax_ref, v_step=0.1, i_min=0.3):
    """Whole-curve errors, normalised by the measured Isc (I-V) and Pmax (P-V).

    Both curves are resampled on a common voltage grid; points where the
    measured current is below i_min (noisy tail near Voc) are left out.
    nMBE keeps the sign: > 0 means SPICE over-predicts."""
    if len(v_meas) < 2 or len(v_sim) < 2 or not isc_ref or not pmax_ref:
        return None
    om, os_ = np.argsort(v_meas), np.argsort(v_sim)
    v_meas, i_meas, v_sim, i_sim = v_meas[om], i_meas[om], v_sim[os_], i_sim[os_]
    lo, hi = max(v_meas[0], v_sim[0]), min(v_meas[-1], v_sim[-1])
    grid = np.arange(lo, hi, v_step)
    if len(grid) < 2:
        return None

    im, is_ = np.interp(grid, v_meas, i_meas), np.interp(grid, v_sim, i_sim)
    mask = im >= i_min
    if mask.sum() < 2:
        mask[:] = True
    out = {"n_grid_pts": int(mask.sum())}
    for tag, d, ref in (("IV", is_ - im, isc_ref), ("PV", grid * (is_ - im), pmax_ref)):
        d = d[mask]
        out[f"nRMSE_{tag}_%"] = float(np.sqrt(np.mean(d ** 2)) / ref * 100)
        out[f"nMAE_{tag}_%"] = float(np.mean(np.abs(d)) / ref * 100)
        out[f"nMBE_{tag}_%"] = float(np.mean(d) / ref * 100)
    return out


def error_summary(errors, tol_pct=5.0):
    """MAE, RMSE, max and share of points within tolerance of a series of relative errors [%]."""
    e = np.asarray(errors, dtype=float)
    return pd.Series({"n_pts": len(e), "MAE_%": np.mean(np.abs(e)), "RMSE_%": np.sqrt(np.mean(e ** 2)),
                      "max_abs_%": np.max(np.abs(e)), "within_tol_%": 100.0 * np.mean(np.abs(e) <= tol_pct)})


# --- resilience indicators ---------------------------------------------------------

def linear_fit(x, y):
    """Slope, intercept and R2 of y = a x + b on the finite points (nan if < 2 points)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 2 or np.ptp(x) == 0:
        return np.nan, np.nan, np.nan
    a, b = np.polyfit(x, y, 1)
    ss_res, ss_tot = np.sum((y - (a * x + b)) ** 2), np.sum((y - y.mean()) ** 2)
    return float(a), float(b), float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan


def gamma_table(df, x_col, group_cols=("shape", "orientation", "module"),
                metrics=("Pmp_loss_%", "Vmp_loss_%", "Imp_loss_%")):
    """Resilience slope gamma = d(loss %)/dx, given per 0.1 of x.
    A larger |gamma| means the loss grows faster with shading, i.e. a less
    resilient module."""
    rows = []
    for keys, sub in df.groupby(list(group_cols)):
        row = dict(zip(group_cols, keys))
        for metric in metrics:
            slope, _b, r2 = linear_fit(sub[x_col], sub[metric])
            short = metric.replace("_loss_%", "")
            row[f"gamma_{short}_per_0.1"] = slope * 0.1
            row[f"R2_{short}"] = r2
        rows.append(row)
    return pd.DataFrame(rows)


# Electrical regimes, from the change of Vmp and Imp at the MPP
#   I   current-limited : dVmp > -15 %          (MPP stays on the high-voltage peak)
#   II  clamped bypass  : dVmp <= -15 %, |dImp| < 3 %
#   III mixed           : dVmp <= -15 %, |dImp| >= 3 %
REGIME_DV_PCT = -15.0
REGIME_DI_PCT = 3.0
REGIME_NAMES = {0: "I- current-limited", 1: "II- clamped bypass", 2: "III- Mixed"}


def regime(d_vmp, d_imp, dv=REGIME_DV_PCT, di=REGIME_DI_PCT):
    if d_vmp is None or d_imp is None or not (np.isfinite(d_vmp) and np.isfinite(d_imp)):
        return np.nan
    if d_vmp > dv:
        return 0.0
    return 1.0 if abs(d_imp) < di else 2.0


def add_regime(df):
    out = df.copy()
    out["regime"] = [regime(v, i) for v, i in zip(out["Vmp_loss_%"], out["Imp_loss_%"])]
    return out


def regime_census(df):
    """Number and share of states in each regime."""
    counts = df["regime"].value_counts().reindex([0.0, 1.0, 2.0], fill_value=0)
    total = int(df["regime"].notna().sum())
    return pd.DataFrame({"regime": [REGIME_NAMES[int(k)] for k in counts.index], "n": counts.values,
                         "share_%": 100.0 * counts.values / total if total else np.nan})


def resilience_summary(sweep):
    """One line per module and band configuration, summarising its loss map:
    mean and worst Pmp loss, resilience slope against tau_eff, and the share
    of states in each electrical regime."""
    sweep = sweep if "regime" in sweep else add_regime(sweep)
    rows = []
    for (module, shape, orientation), sub in sweep.groupby(["module", "shape", "orientation"]):
        slope, _b, r2 = linear_fit(sub["tau_eff"], sub["Pmp_loss_%"])
        n = sub["regime"].notna().sum()
        rows.append({
            "module": module, "shape": shape, "orientation": orientation,
            "mean_Pmp_loss_%": sub["Pmp_loss_%"].mean(),
            "worst_Pmp_loss_%": sub["Pmp_loss_%"].min(),
            "gamma_Pmp_per_0.1_tau_eff": slope * 0.1, "R2": r2,
            **{f"share_{REGIME_NAMES[k].split('-')[0]}_%": 100.0 * (sub["regime"] == k).sum() / n if n else np.nan
               for k in (0, 1, 2)},
        })
    return pd.DataFrame(rows)


# --- multi-peak P-V curves ------------------------------------------------------------

def pv_peaks(v, i, min_prominence=0.01, min_height=0.02, n_grid=2000):
    """Local maxima of the P-V curve under partial shading.

    P_highV: peak at the highest voltage (all sub-strings active).
    P_lowV : largest peak at lower voltage (bypass diode conducting).
    R_peaks = P_lowV / P_highV, 0 if there is a single peak.
    Two neighbouring maxima are merged when the valley between them is
    shallower than min_prominence * Pmax."""
    v, i = np.asarray(v, dtype=float), np.asarray(i, dtype=float)
    ok = np.isfinite(v) & np.isfinite(i)
    v, i = v[ok], i[ok]
    empty = dict(n_peaks=0, P_peak_highV_W=np.nan, V_peak_highV_V=np.nan,
                 P_peak_lowV_W=np.nan, V_peak_lowV_V=np.nan, R_peaks=np.nan)
    if len(v) < 5:
        return empty
    order = np.argsort(v)
    vg = np.linspace(v[order][0], v[order][-1], n_grid)
    pg = vg * np.interp(vg, v[order], i[order])
    pmax = pg.max()
    if pmax <= 0:
        return empty

    peaks = [k for k in range(1, n_grid - 1)
             if pg[k] > pg[k - 1] and pg[k] >= pg[k + 1] and pg[k] >= min_height * pmax] or [int(np.argmax(pg))]
    merged = True
    while merged and len(peaks) > 1:
        merged = False
        for a, b in zip(peaks[:-1], peaks[1:]):
            if min(pg[a], pg[b]) - pg[a:b + 1].min() < min_prominence * pmax:
                peaks.remove(a if pg[a] < pg[b] else b)
                merged = True
                break

    high = max(peaks, key=lambda k: vg[k])
    out = dict(n_peaks=len(peaks), P_peak_highV_W=float(pg[high]), V_peak_highV_V=float(vg[high]),
               P_peak_lowV_W=np.nan, V_peak_lowV_V=np.nan, R_peaks=0.0)
    others = [k for k in peaks if k != high]
    if others:
        low = max(others, key=lambda k: pg[k])
        out.update(P_peak_lowV_W=float(pg[low]), V_peak_lowV_V=float(vg[low]), R_peaks=float(pg[low] / pg[high]))
    return out
