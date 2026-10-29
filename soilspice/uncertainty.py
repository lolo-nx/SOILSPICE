"""Model uncertainty of the SPICE losses, estimated on the flash-test campaign.

For every shaded flash i the residual of the loss is

    e_i = dX_SPICE,i - dX_EXP,i          [% points, X = Pmp, Vmp or Imp]

and the residuals are described as a bias plus a random part,

    e_i = b + eps_i,   b = MBE = mean(e_i),   sigma = std(e_i - b).

Since dX_EXP = dX_SPICE - e, a simulated loss is turned into an estimate of
the real loss with

    dX_real = dX_SPICE - b  +/-  z * sigma          (z = 1.96 for 95 %)

The bias is SUBTRACTED: a positive b means SPICE predicts a smaller loss
(less negative) than measured, so the estimate must be moved down.

b and sigma can be global or depend on a factor (electrical regime, module,
...). homogeneity() says whether a single global sigma is justified.
"""
import numpy as np
import pandas as pd

from . import config, experiments, geometry
from .analysis import REGIME_NAMES, linear_fit, regime
from .validation import loss_comparison

METRICS = ("Pmp", "Vmp", "Imp")
Z95 = 1.96
MEASURED_TAUS = tuple(sorted({entry["tau_I"] for entries in experiments.CAMPAIGN.values()
                              for entry in entries if entry["tau_I"] < 1.0}, reverse=True))


def tau_class(tau):
    """Nearest shaded transmittance configured in experiments.py."""
    return f"{min(MEASURED_TAUS, key=lambda t: abs(t - tau)):g}"


def regime_name(code):
    return REGIME_NAMES[int(code)] if np.isfinite(code) else "n/a"


# --- residuals ----------------------------------------------------------------------

def residuals(results):
    """One row per shaded flash with a simulation: losses, residuals and the
    variables the residuals are checked against. The regime is computed from
    the SPICE losses, as it is for the loss maps the model is applied to."""
    rows = []
    for r in results:
        if not r.simulated:
            continue
        severity, _ = geometry.band_severity(r.module, r.fcov, r.shape, r.orientation)
        for _, row in loss_comparison([r]).iterrows():
            if not np.isfinite(row.get("Pmp_loss_sim_%", np.nan)):
                continue
            out = {"module": r.module.name, "shape": r.shape, "orientation": r.orientation,
                   "tau_I": row["tau_I"], "tau_class": tau_class(row["tau_I"]), "A_cov": r.acov,
                   "tau_eff": geometry.tau_eff(severity, row["tau_I"])}
            for m in METRICS:
                out[f"{m}_loss_exp_%"] = row[f"{m}_loss_exp_%"]
                out[f"{m}_loss_sim_%"] = row[f"{m}_loss_sim_%"]
                out[f"e_{m}"] = row[f"{m}_loss_sim_%"] - row[f"{m}_loss_exp_%"]
            code = regime(row["Vmp_loss_sim_%"], row["Imp_loss_sim_%"])
            out["regime"] = code
            out["regime_name"] = regime_name(code)
            rows.append(out)
    return pd.DataFrame(rows)


# --- statistics -----------------------------------------------------------------------

def _t975(dof):
    try:
        from scipy.stats import t
        return float(t.ppf(0.975, dof))
    except ImportError:
        return np.nan


def stats(e, z=Z95):
    """Bias, spread and resulting interval of a set of residuals [% points].

    half_width     = z * sigma, the interval used in the energy model;
    half_width_t   = t(0.975, n-1) * sigma * sqrt(1 + 1/n), the same interval
                     corrected for the small sample (bias and sigma are themselves
                     estimated); given for information."""
    e = np.asarray(e, dtype=float)
    e = e[np.isfinite(e)]
    n = len(e)
    if n == 0:
        return pd.Series({"n": 0}, dtype=float)
    b = e.mean()
    sigma = e.std(ddof=1) if n > 1 else np.nan
    return pd.Series({
        "n": n, "MBE": b, "MAE": np.abs(e).mean(), "RMSE": np.sqrt(np.mean(e ** 2)),
        "sigma": sigma, "max_abs": np.abs(e).max(),
        "half_width": z * sigma,
        "half_width_t": _t975(n - 1) * sigma * np.sqrt(1 + 1 / n) if n > 1 else np.nan,
    })


def summary(res, metric="Pmp", by=None):
    """stats() for all residuals, or per group."""
    col = f"e_{metric}"
    if not by:
        return stats(res[col]).to_frame().T.assign(group="all")
    by = [by] if isinstance(by, str) else list(by)
    return res.groupby(by)[col].apply(stats).unstack().reset_index()


def _brown_forsythe(groups):
    """p-value of the Brown-Forsythe test (equal variances), NaN without scipy
    or with fewer than two groups of at least two points."""
    groups = [np.asarray(g, dtype=float) for g in groups if len(g) >= 2]
    if len(groups) < 2:
        return np.nan
    try:
        from scipy.stats import levene
    except ImportError:
        return np.nan
    return float(levene(*groups, center="median").pvalue)


def homogeneity(res, metric="Pmp", factors=("module", "shape", "orientation", "regime_name", "tau_class")):
    """Bias and spread of the residuals for each level of each factor,
    compared with the global values.

    sigma_ratio = sigma of the level / global sigma. p_equal_var is the
    Brown-Forsythe p-value of the factor: below 0.05, the spread is not the
    same for all levels and a single global sigma is questionable."""
    col = f"e_{metric}"
    glob = stats(res[col])
    frames = []
    for factor in factors:
        table = res.groupby(factor)[col].apply(stats).unstack().reset_index().rename(columns={factor: "level"})
        table.insert(0, "factor", factor)
        table["sigma_ratio"] = table["sigma"] / glob["sigma"]
        table["p_equal_var"] = _brown_forsythe([g[col] for _, g in res.groupby(factor)])
        frames.append(table)
    return pd.concat(frames, ignore_index=True)


def trends(res, metric="Pmp", variables=("A_cov", "tau_I", "tau_eff")):
    """Linear trend of the residuals against continuous variables. A clear
    slope means the bias is not constant over the range of the loss maps."""
    col = f"e_{metric}"
    rows = []
    for var in variables:
        slope, intercept, r2 = linear_fit(res[var], res[col])
        rows.append({"variable": var, "slope": slope, "intercept": intercept, "R2": r2,
                     "min": res[var].min(), "max": res[var].max()})
    return pd.DataFrame(rows)


# --- model and application to the loss maps ---------------------------------------------

def model(res, metric="Pmp", by=None, min_n=5, z=Z95):
    """Bias and sigma to apply to the simulated losses.

    by=None gives one global row. With a factor (e.g. "regime"), one row per
    level plus a last row "all" with the global values; a level with fewer
    than min_n residuals takes the global values (column 'source' says which)."""
    glob = stats(res[f"e_{metric}"], z)
    if not by:
        return pd.DataFrame([{"metric": metric, "group": "all", "n": glob["n"], "bias": glob["MBE"],
                              "sigma": glob["sigma"], "half_width": glob["half_width"], "source": "global"}])
    rows = []
    for level, sub in res.groupby(by):
        s = stats(sub[f"e_{metric}"], z)
        own = s["n"] >= min_n
        rows.append({"metric": metric, by: level, "n": s["n"],
                     "bias": s["MBE"] if own else glob["MBE"],
                     "sigma": s["sigma"] if own else glob["sigma"],
                     "half_width": (s["half_width"] if own else glob["half_width"]),
                     "source": by if own else "global (n < %d)" % min_n})
    rows.append({"metric": metric, by: "all", "n": glob["n"], "bias": glob["MBE"], "sigma": glob["sigma"],
                 "half_width": glob["half_width"], "source": "global"})
    return pd.DataFrame(rows)


def apply(sweep, unc, metric="Pmp", by=None, z=Z95):
    """Add the corrected loss and its 95 % interval to a loss-map table.

    New columns: <metric>_bias, <metric>_sigma, <metric>_loss_corr_%,
    <metric>_loss_lo_%, <metric>_loss_hi_%. Levels of `by` absent from the
    validation take the global values."""
    out = sweep.copy()
    if by == "regime" and "regime" not in out:
        out["regime"] = [regime(v, i) for v, i in zip(out["Vmp_loss_%"], out["Imp_loss_%"])]
    if by:
        table = unc.set_index(by)
        b = out[by].map(table["bias"]).fillna(table.loc["all", "bias"])
        s = out[by].map(table["sigma"]).fillna(table.loc["all", "sigma"])
    else:
        b = pd.Series(float(unc["bias"].iloc[0]), index=out.index)
        s = pd.Series(float(unc["sigma"].iloc[0]), index=out.index)
    loss = out[f"{metric}_loss_%"]
    out[f"{metric}_bias"] = b
    out[f"{metric}_sigma"] = s
    out[f"{metric}_loss_corr_%"] = loss - b
    out[f"{metric}_loss_lo_%"] = loss - b - z * s
    out[f"{metric}_loss_hi_%"] = loss - b + z * s
    return out


def coverage(res, unc, metric="Pmp", by=None, z=Z95):
    """Share of the measured losses that fall inside the interval built from
    their own simulation (in-sample check: should be close to 95 %)."""
    sim = res.rename(columns={f"{metric}_loss_sim_%": f"{metric}_loss_%"})
    est = apply(sim, unc, metric, by, z)
    inside = ((res[f"{metric}_loss_exp_%"] >= est[f"{metric}_loss_lo_%"])
              & (res[f"{metric}_loss_exp_%"] <= est[f"{metric}_loss_hi_%"]))
    return 100.0 * inside.mean()


def global_model(metric="Pmp", z=Z95):
    """Global bias and sigma from results/validation/residuals.csv, or None if
    the validation has not been run. Used to print the model uncertainty on
    the simulation figures (indicative)."""
    path = config.VALIDATION_DIR / "residuals.csv"
    if not path.exists():
        return None
    s = stats(pd.read_csv(path)[f"e_{metric}"], z)
    if not s["n"]:
        return None
    return {"metric": metric, "n": int(s["n"]), "bias": float(s["MBE"]), "sigma": float(s["sigma"]), "z": z}
