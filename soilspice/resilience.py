"""Soiling resilience metrics, computed on the SPICE loss maps.

All metrics use the optical loss of the band,

    S = A_cov (1 - tau_I)                           [% of the module light]

which is the share of the light lost by the module: a band twice as large
and half as opaque gives the same S. With uniform soiling the power loss
would simply equal S, so every metric compares the real loss with S.

  M     loss amplification, M = |dPmp / Pmp,clean| / S. Percent of power lost
        per percent of light blocked. 1 = uniform soiling, >> 1 = the series
        mismatch amplifies the loss.
  M0    initial amplification: slope of |dPmp| against S for light soiling
        (S <= 2 %), least squares through the origin. What happens at the
        first stage of soiling.
  S*    critical optical loss: smallest S at which a bypass diode conducts
        (dVmp <= -15 %, regime II or III). NaN: no bypass in the simulated range.
  max loss  worst power loss of the map (A_cov = 10 %, tau_I = 0).
  P_an,S    average normalized power under soiling (Klasen et al. 2022,
        Fauzan et al. 2023) with S instead of the shaded-area fraction:
            P_an,S = int P/Pclean dS / int (1 - S) dS   over [0, S_ref]
        P/Pclean is first averaged in 1 % bins of S, as Klasen et al. average
        over all shading scenarios with the same shaded fraction. 1 = ideal
        (loss = light blocked), lower = less resilient. Always computed on the
        same window for every configuration, otherwise the values cannot be
        compared (P_an,S tends to 1 on a short window).
"""
import numpy as np
import pandas as pd

from .analysis import regime

_trapz = getattr(np, "trapezoid", None) or np.trapz

LOW_S_MAX = 2.0                       # [%] window of M0
PAN_BIN = 1.0                         # [%] bin width of the mean power curve
PAN_NODES = 201                       # integration nodes of P_an,S
S_TOL = 0.005                         # relative tolerance on a window edge (A_cov from bisection)
PAN_WINDOWS = (1, 2, 3, 4, 5, 6, 8, 10)   # [%] windows [0, S_ref] of the P_an,S table
KEYS = ["module", "shape", "orientation"]


def add_optical_loss(sweep):
    """Columns S_% (optical loss, from the achieved A_cov), M and regime."""
    out = sweep.copy()
    out["S_%"] = 100.0 * out["A_cov"] * (1.0 - out["tau_I"])
    with np.errstate(divide="ignore", invalid="ignore"):
        out["M"] = np.where(out["S_%"] > 0, -out["Pmp_loss_%"] / out["S_%"], np.nan)
    if "regime" not in out:
        out["regime"] = [regime(v, i) for v, i in zip(out["Vmp_loss_%"], out["Imp_loss_%"])]
    return out


def initial_amplification(S, loss, s_max=LOW_S_MAX):
    """M0: slope of -loss against S through the origin, for S <= s_max."""
    low = (S <= s_max) & np.isfinite(loss) & (S > 0)
    if not low.any():
        return np.nan
    return float(np.sum(S[low] * -loss[low]) / np.sum(S[low] ** 2))


def power_curve(S, loss, s_max=None, bin_width=PAN_BIN):
    """Mean normalized power P/Pclean in bins of S, starting at (0, 1)."""
    S = np.round(np.asarray(S, dtype=float), 2)       # same bin edges for every module
    P = 1.0 + np.asarray(loss, dtype=float) / 100.0
    ok = np.isfinite(S) & np.isfinite(P)
    if s_max is not None:
        ok &= S <= s_max * (1 + S_TOL) + 1e-9
    S, P = S[ok], P[ok]
    if S.size == 0:
        return np.array([0.0]), np.array([1.0])
    edges = np.arange(0.0, S.max() + bin_width + 1e-9, bin_width)
    idx = np.clip(np.digitize(S, edges) - 1, 0, len(edges) - 2)
    xs, ys = [0.0], [1.0]
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.any():
            xs.append(S[sel].mean())
            ys.append(P[sel].mean())
    return np.array(xs), np.array(ys)


def average_power(S, loss, s_ref):
    """P_an,S on [0, s_ref]. NaN if the configuration does not reach s_ref."""
    S = np.asarray(S, dtype=float)
    if not np.isfinite(S).any() or s_ref > np.nanmax(S) * (1 + S_TOL):
        return np.nan
    x, y = power_curve(S, loss, s_max=s_ref)
    if x.size < 2:
        return np.nan
    if x[-1] < s_ref:                       # last bin centre below s_ref: extend the last segment
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2]) if x[-1] > x[-2] else 0.0
        y = np.r_[y, max(0.0, y[-1] + slope * (s_ref - x[-1]))]
        x = np.r_[x, s_ref]
    nodes = np.linspace(0.0, s_ref, PAN_NODES)
    yn = np.interp(nodes, x, y)
    xf = nodes / 100.0
    return float(_trapz(yn, xf) / _trapz(1.0 - xf, xf))


def summary(sweep, windows=PAN_WINDOWS):
    """One row per module and band configuration: M0, S*, max loss, mean and
    largest M, share of bypass states, P_an,S on the full common window and on
    the windows [0, S_ref]."""
    df = add_optical_loss(sweep)
    s_full = min(np.nanmax(g["S_%"]) for _, g in df.groupby(KEYS))     # common window for P_an,S
    rows = []
    for (module, shape, orientation), g in df.groupby(KEYS, sort=False):
        S, loss = g["S_%"].to_numpy(float), g["Pmp_loss_%"].to_numpy(float)
        bypass = g["regime"].to_numpy(float) > 0
        rows.append({
            "module": module, "shape": shape, "orientation": orientation,
            "M0": initial_amplification(S, loss),
            "S_star_%": float(S[bypass].min()) if bypass.any() else np.nan,
            "max_loss_%": float(-np.nanmin(loss)),
            "M_mean": float(np.nanmean(g["M"])), "M_max": float(np.nanmax(g["M"])),
            "bypass_share_%": 100.0 * bypass.mean(),
            "P_an_S": average_power(S, loss, s_full),
            **{f"P_an_S_0to{w:g}": average_power(S, loss, w) for w in windows},
        })
    out = pd.DataFrame(rows)
    out.attrs["s_full"] = s_full
    return out


# --- metrics at equal band height -------------------------------------------------------
# Same band height h on every module (margins, cell height and module size
# included), from scripts/run_equal_height.py:
#   h1, h5   band height [mm] at which the Pmp loss reaches 1 % and 5 %
#            (linear interpolation between the simulated heights). h5 / growth
#            rate of the band gives the time before a 5 % loss.
#   P_an,h   mean normalized power over the optical band height
#            h_opt = h (1 - tau_I), all h and tau_I together (as S for P_an,S):
#            P_an,h = (1/h_max) int_0^h_max P(h_opt)/P_clean dh_opt   (1 = no loss)
H_MAX = 100.0
LOSS_THRESHOLDS = (1.0, 5.0)


def height_at_loss(h, loss, threshold):
    """First band height where -loss reaches `threshold` [%]; NaN if never
    reached in the simulated range."""
    h, lost = np.r_[0.0, h], np.r_[0.0, -np.asarray(loss, dtype=float)]
    lost = np.maximum.accumulate(lost)                 # first crossing
    if lost[-1] < threshold:
        return np.nan
    k = int(np.argmax(lost >= threshold))
    return float(np.interp(threshold, lost[k - 1:k + 1], h[k - 1:k + 1]))


def average_power_height(h, loss, h_max=H_MAX):
    """P_an,h on [0, h_max]; NaN if the heights do not reach h_max."""
    h = np.r_[0.0, np.asarray(h, dtype=float)]
    p = np.r_[1.0, 1.0 + np.asarray(loss, dtype=float) / 100.0]
    if h.max() < h_max - 1e-9:
        return np.nan
    nodes = np.linspace(0.0, h_max, 401)
    return float(_trapz(np.interp(nodes, h, p), nodes) / h_max)


OPT_BIN_MM = 5.0     # bin width of the mean power curve against the optical band height


def optical_height(h, tau):
    """Optical band height h_opt = h (1 - tau_I) [mm]: the height of an opaque
    band that blocks the same light as a band of height h and transmittance
    tau_I. It is the equal-height counterpart of S = A_cov (1 - tau_I): S
    counts the blocked light on the cells (module-dependent), h_opt counts the
    soil on the module edge (the same for every module)."""
    return np.asarray(h, dtype=float) * (1.0 - np.asarray(tau, dtype=float))


def average_power_optical_height(h_opt, loss, h_max=H_MAX, bin_mm=OPT_BIN_MM):
    """P_an,h: mean normalized power over h_opt in [0, h_max]. P/Pclean is first
    averaged in bins of h_opt (all band heights and tau_I together, as for
    P_an,S). NaN if h_opt does not reach h_max."""
    h_opt = np.asarray(h_opt, dtype=float)
    if np.nanmax(h_opt) < h_max * (1 - S_TOL):
        return np.nan
    x, y = power_curve(h_opt, loss, s_max=h_max, bin_width=bin_mm)
    if x[-1] < h_max:
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2]) if len(x) > 1 and x[-1] > x[-2] else 0.0
        y = np.r_[y, max(0.0, y[-1] + slope * (h_max - x[-1]))]
        x = np.r_[x, h_max]
    nodes = np.linspace(0.0, h_max, 401)
    return float(_trapz(np.interp(nodes, x, y), nodes) / h_max)


def height_summary(df, taus=(0.5, 0.0), h_max=H_MAX, thresholds=LOSS_THRESHOLDS):
    """One row per module, configuration and tau_I in `taus`:
      h1, h5     band height [mm] for 1 % and 5 % loss at that tau_I
      P_an_h     mean normalized power over the optical band height
                 h_opt = h (1 - tau_I) in [0, h_max], all tau_I together
                 (same value on every row of a configuration)
      P_an_h_tau mean normalized power over h in [0, h_max] at that tau_I only"""
    rows = []
    for (module, shape, orientation), g in df.groupby(KEYS, sort=False):
        pan = average_power_optical_height(optical_height(g["height_mm"], g["tau_I"]),
                                           g["Pmp_loss_%"].to_numpy(float), h_max)
        for tau in taus:
            gt = g[np.isclose(g["tau_I"], tau)].sort_values("height_mm")
            if gt.empty:
                continue
            h, loss = gt["height_mm"].to_numpy(float), gt["Pmp_loss_%"].to_numpy(float)
            rows.append({"module": module, "shape": shape, "orientation": orientation, "tau_I": tau,
                         **{f"h{t:g}_mm": height_at_loss(h, loss, t) for t in thresholds},
                         "P_an_h": pan,
                         "P_an_h_tau": average_power_height(h, loss, h_max),
                         f"loss_at_{h_max:g}mm_%": float(np.interp(h_max, h, loss)) if h.max() >= h_max else np.nan})
    out = pd.DataFrame(rows)
    out.attrs["h_max"] = h_max
    return out
