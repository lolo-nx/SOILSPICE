"""Roof field test: two twin SunPower modules in landscape, one clean (Lorraine 1) and one with the
80 mm rectangular band of one filter layer (Lorraine 2), compared with SPICE.

    soilspice.field.run(csv, out) (chapters/3_field_test.py)

The two modules are not identical, so the ratio of their powers only means something relative
to the same ratio before the band (calibration r0):

    r(t)    = (P_band / G_band) / (P_clean / G_clean)      pairs of records < 90 s apart
    r0      = median of r before the band
    dP(t)   = r(t) / r0 - 1                                the band's own effect

The clean module is the weather model (same roof, same instant), so SPICE needs no irradiance or
temperature model:

    P_model(t) = P_clean(t) x G_band/G_clean x r0 x (1 + dP_SPICE)

dP_SPICE is the STC loss of the equal height sweep (SunPower, rectangle, landscape, h = 80 mm),
interpolated at tau_I = 0.71. The records taken while the band was installed are dropped.

Outputs (out/):
    pairs.csv         every pair with its ratio and loss
    daily.csv         daily energy of the band module: expected clean, measured, SPICE
    loss_vs_G.csv     median loss per irradiance bin
    loss_by_hour.csv  median loss per half hour of the day
    summary.json      r0, SPICE loss, measured power and energy loss
    fig_field_test        (a) power of the band module, (b) daily energy, (c) daily energy lost
    fig_field_test_loss   (a) loss over time, (b) loss over the day, (c) cumulative energy
"""

import json

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from . import config, plots
from .style import LETTER_KW, STYLE

SWEEP_CSV = config.EQUAL_HEIGHT_CSV
UNC_CSV = config.UNCERTAINTY_DIR / "model_Pmp.csv"

CLEAN, BAND = "Lorraine 1", "Lorraine 2"
BAND_FROM, BAND_TO = pd.Timestamp("2026-09-28 15:30"), pd.Timestamp("2026-09-28 15:45")
H_MM, TAU_I = 80.0, 0.71               # band of the field test (one filter layer)
LOSS_MEASURED_STC = None              # optional: the indoor flash loss of this band, in %

MAX_PAIR_S = 90                        # the tracer sweeps the two modules about 30 s apart
G_MIN = 400.0                          # W/m2, power ratios (clear sky, stable)
G_ENERGY = 20.0                        # W/m2, energy integration (keeps mornings and evenings)
MAX_DG = 0.03                          # relative irradiance gap inside a pair (passing clouds)
MAX_DG_ENERGY = 0.10                   # looser for the energy, so that cloudy hours still count
SHADE_FROM_H = 15                      # the clean module is shaded from about 17:20 (local time)
SHADE_PR = 0.95                        # shaded if P/G of the clean module < 95 % of its midday value
MAX_GAP_MIN = 15                       # longer gaps between pairs are not integrated
G_BINS = np.arange(100, 1101, 100)

COL_CLEAN = plots.TAU_ANCHORS[1.00][0]          # same colours as the validation grid
COL_BAND, MK_BAND = plots.TAU_ANCHORS[0.70]


# --- data --------------------------------------------------------------------------

def load(path):
    d = pd.read_csv(path, parse_dates=["RecTime"])
    d = d.dropna(subset=["Pmpp", "Irr"])
    d = d[(d["Pmpp"] > 0) & (d["Irr"] > 0)]
    return {n: d[d["Name"] == n][["RecTime", "Pmpp", "Irr"]].sort_values("RecTime") for n in (CLEAN, BAND)}


def pair(series):
    """Each band record with the nearest clean record, within MAX_PAIR_S."""
    b = series[BAND].rename(columns={"RecTime": "t", "Pmpp": "P_b", "Irr": "G_b"})
    c = series[CLEAN].rename(columns={"RecTime": "t_c", "Pmpp": "P_c", "Irr": "G_c"})
    p = pd.merge_asof(b, c, left_on="t", right_on="t_c", direction="nearest",
                      tolerance=pd.Timedelta(seconds=MAX_PAIR_S)).dropna(subset=["P_c"])
    p = p[(p["t"] < BAND_FROM) | (p["t"] > BAND_TO)].copy()
    p["phase"] = np.where(p["t"] < BAND_FROM, "before", "after")
    p["dG"] = (p["G_b"] - p["G_c"]).abs() / (0.5 * (p["G_b"] + p["G_c"]))
    p["r"] = (p["P_b"] / p["G_b"]) / (p["P_c"] / p["G_c"])
    # late afternoon, an obstacle shades the clean module (its own P/G collapses, the band module's
    # does not): those pairs are dropped everywhere, before and after the band
    pr = p["P_c"] / p["G_c"]
    ref = pr[p["t"].dt.hour.between(11, 14)].groupby(p["t"].dt.date).median()
    p["clean_shaded"] = (p["t"].dt.hour >= SHADE_FROM_H) & (pr < SHADE_PR * p["t"].dt.date.map(ref))
    p["stable"] = (p[["G_b", "G_c"]].min(axis=1) >= G_MIN) & (p["dG"] < MAX_DG) & ~p["clean_shaded"]
    return p.reset_index(drop=True)


def spice_loss():
    """STC loss of the band [%] (negative) and the model uncertainty (95 % half width, bias)."""
    s = pd.read_csv(SWEEP_CSV)
    s = s[(s["module"] == "SUNPOWER") & (s["shape"] == "rectangle") & (s["orientation"] == "landscape")
          & (s["height_mm"] == H_MM)].sort_values("tau_I")
    dp = float(np.interp(TAU_I, s["tau_I"], s["Pmp_loss_%"]))
    u = pd.read_csv(UNC_CSV).iloc[0]
    return dp, float(u["half_width"]), float(u["bias"])


def integrate(t, p):
    """Energy [Wh] of a power series on irregular times, gaps longer than MAX_GAP_MIN skipped."""
    t = t.to_numpy().astype("datetime64[s]").astype(float) / 3600.0
    dt = np.diff(t)
    ok = dt <= MAX_GAP_MIN / 60
    mid = 0.5 * (p[1:] + p[:-1])
    return float(np.sum(mid[ok] * dt[ok]))


def daily_energy(p, r0, dp):
    """Daily energy of the band module after the band: expected clean, measured and SPICE."""
    a = p[(p["phase"] == "after") & (p[["G_b", "G_c"]].min(axis=1) >= G_ENERGY) & ~p["clean_shaded"]
          & (p["dG"] < MAX_DG_ENERGY)].copy()
    a["P_exp"] = a["P_c"] * a["G_b"] / a["G_c"] * r0       # band module if it were clean
    a["P_mod"] = a["P_exp"] * (1 + dp / 100)
    rows = []
    for day, g in a.groupby(a["t"].dt.date):
        if day == BAND_TO.date():                            # band installed at 15:30, too little left
            continue
        e = {k: integrate(g["t"], g[c].to_numpy()) for k, c in
             (("E_clean_Wh", "P_exp"), ("E_meas_Wh", "P_b"), ("E_spice_Wh", "P_mod"))}
        rows.append({"day": pd.Timestamp(day), **e, "n_pairs": len(g)})
    d = pd.DataFrame(rows)
    d["loss_meas_%"] = 100 * (d["E_meas_Wh"] / d["E_clean_Wh"] - 1)
    last = p.groupby(p["t"].dt.date)["t"].max()
    d["partial"] = [last[x.date()].hour < 17 for x in d["day"]]   # export stopped during the day
    return d[~d["partial"]].reset_index(drop=True)


def loss_by_hour(p, r0):
    """Median loss per hour of the day with the band (G >= 100 W/m2), with its quartiles."""
    a = p[(p["phase"] == "after") & ~p["clean_shaded"] & (p[["G_b", "G_c"]].min(axis=1) >= 100)].copy()
    a["dP"] = 100 * (a["r"] / r0 - 1)
    a["hour"] = a["t"].dt.hour + a["t"].dt.minute / 60
    a["bin"] = (a["hour"] // 0.5) * 0.5 + 0.25
    g = a.groupby("bin")["dP"]
    out = pd.DataFrame({"hour": g.median().index, "median": g.median().values, "q25": g.quantile(0.25).values,
                        "q75": g.quantile(0.75).values, "n": g.size().values})
    return out[out["n"] >= 10]


def bypass_plateau():
    """SPICE loss of the same band when it is opaque [%]: the plateau of regime II in landscape."""
    s = pd.read_csv(SWEEP_CSV)
    s = s[(s["module"] == "SUNPOWER") & (s["shape"] == "rectangle") & (s["orientation"] == "landscape")
          & (s["height_mm"] == H_MM) & (s["tau_I"] == 0.0)]
    return float(s["Pmp_loss_%"].iloc[0])


def loss_vs_irradiance(p, r0):
    a = p[(p["phase"] == "after") & (p["dG"] < MAX_DG) & ~p["clean_shaded"]].copy()
    a["dP"] = 100 * (a["r"] / r0 - 1)
    a["G"] = 0.5 * (a["G_b"] + a["G_c"])
    a["bin"] = pd.cut(a["G"], G_BINS)
    g = a.groupby("bin", observed=True)["dP"]
    out = pd.DataFrame({"G_mid": [b.mid for b in g.median().index], "median": g.median().values,
                        "q25": g.quantile(0.25).values, "q75": g.quantile(0.75).values,
                        "n": g.size().values})
    return out[out["n"] >= 10]


# --- figure ------------------------------------------------------------------------

def power_series(p, r0, dp):
    """Power of the band module on the pair grid: expected clean, measured, SPICE (NaN in gaps and shade)."""
    a = p[(p["phase"] == "after") & ~p["clean_shaded"] & (p["dG"] < MAX_DG_ENERGY)].copy()
    a["P_exp"] = a["P_c"] * a["G_b"] / a["G_c"] * r0
    a["P_mod"] = a["P_exp"] * (1 + dp / 100)
    gap = a["t"].diff() > pd.Timedelta(minutes=MAX_GAP_MIN)
    a.loc[gap, ["P_exp", "P_mod", "P_b"]] = np.nan          # breaks the lines at night and in gaps
    return a


def figure_energy(p, r0, dp, daily, path):
    """Main figure: (a) power of the band module, (b) daily energy, (c) daily energy lost."""
    a = power_series(p, r0, dp)
    a = a[a["t"] < daily["day"].max() + pd.Timedelta(days=1)]
    col_spice = plots.SIM_COLOR
    days = daily["day"] + pd.Timedelta(hours=12)
    alpha = np.where(daily["partial"], 0.35, 1.0)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(3, 1, figsize=(7.2, 5.6), sharex=True,
                                 gridspec_kw={"height_ratios": [1.3, 1, 1]})
        t0 = daily["day"].min()
        t1 = daily["day"].max() + pd.Timedelta(days=1)
        for ax in axes:
            for d in pd.date_range(t0, t1, freq="D")[:-1:2]:
                ax.axvspan(d, d + pd.Timedelta(days=1), color="0.94", lw=0, zorder=0)

        # (a) power
        ax = axes[0]
        ax.plot(a["t"], a["P_exp"], color=COL_CLEAN, lw=0.8, ls=(0, (3, 1.5)), label="SPR1, without band")
        ax.plot(a["t"], a["P_b"], color=COL_BAND, lw=1.6, alpha=0.55, label="SPR2, band, measured")
        ax.plot(a["t"], a["P_mod"], color=col_spice, lw=0.7,
                label=rf"SPR2, band, SPICE ($\Delta P_{{mp}}$ = {dp:.1f} %)")
        ax.set_ylabel("Power [W]")
        ax.set_ylim(0, 1.25 * np.nanmax(a["P_exp"]))
        ax.legend(loc="upper left", ncol=3, frameon=False, handlelength=1.8, columnspacing=1.2)

        # (b) daily energy
        ax = axes[1]
        w = pd.Timedelta(hours=6.5)
        for k, (col, c, _) in enumerate((("E_clean_Wh", COL_CLEAN, "SPR1, without band"),
                                           ("E_meas_Wh", COL_BAND, "SPR2, band, measured"),
                                           ("E_spice_Wh", col_spice, "SPR2, band, SPICE"))):
            for x, y, al in zip(days + (k - 1) * w, daily[col], alpha):
                ax.bar(x, y, width=w, color=c, alpha=al, lw=0)
        handles = [Patch(color=c, label=lab) for c, lab in
                   ((COL_CLEAN, "SPR1, without band"), (COL_BAND, "SPR2, band, measured"),
                    (col_spice, "SPR2, band, SPICE"))]
        ax.set_ylabel("Daily energy [Wh]")
        ax.set_ylim(0, 1.3 * daily["E_clean_Wh"].max())
        ax.legend(handles=handles, loc="upper left", ncol=4, frameon=False, handlelength=1.2, columnspacing=1.2)

        # (c) daily energy lost
        ax = axes[2]
        lost_m = daily["E_clean_Wh"] - daily["E_meas_Wh"]
        lost_s = daily["E_clean_Wh"] - daily["E_spice_Wh"]
        w = pd.Timedelta(hours=8)
        for x, y, ys, al in zip(days, lost_m, lost_s, alpha):
            ax.bar(x - w / 2, y, width=w, color=COL_BAND, alpha=al, lw=0)
            ax.bar(x + w / 2, ys, width=w, color=col_spice, alpha=al, lw=0)
        ax.set_ylabel("Energy lost [Wh]")
        ax.set_ylim(0, 1.25 * max(lost_m.max(), lost_s.max()))
        ax.legend(handles=[Patch(color=COL_BAND, label="measured"), Patch(color=col_spice, label="SPICE")],
                  loc="upper left", ncol=2, frameon=False, handlelength=1.2)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.set_xlim(t0, t1)
        ax.set_xlabel("Date (2026, local time)")

        for i, ax in enumerate(axes):
            ax.text(-0.075, 1.0, "abc"[i], transform=ax.transAxes, va="top", **LETTER_KW)
            ax.tick_params(pad=1.5)
        fig.align_ylabels(axes)
        fig.tight_layout(h_pad=0.6)
        plots.save(fig, path)


def figure(p, r0, dp, half, daily, hours, plateau, path):
    st = p[p["stable"]].copy()
    st["dP"] = 100 * (st["r"] / r0 - 1)
    g = st.groupby(st["t"].dt.floor("D"))["dP"]
    med = g.median()[g.size() >= 10]                      # days with a few clear pairs are noise
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={"width_ratios": [1.5, 1, 1]})

        # (a) loss over time
        ax = axes[0]
        ax.plot(st["t"], st["dP"], ls="none", marker=".", ms=1.5, color="0.75", rasterized=True)
        for phase, c in (("before", COL_CLEAN), ("after", COL_BAND)):
            m = med[(med.index < BAND_FROM.floor("D")) if phase == "before" else (med.index > BAND_FROM.floor("D"))]
            ax.plot(m.index + pd.Timedelta(hours=12), m.values, ls="none", marker="o" if phase == "before"
                    else MK_BAND, ms=3.5, color=c, mec="black", mew=0.3, zorder=4)
        t1 = st["t"].max()
        ax.fill_between([BAND_TO, t1], dp - half, dp + half, color=COL_BAND, alpha=0.15, lw=0)
        ax.plot([BAND_TO, t1], [dp, dp], ls="--", lw=1.1, color=COL_BAND)
        ax.plot([st["t"].min(), BAND_FROM], [0, 0], ls="--", lw=1.1, color=COL_CLEAN)
        ax.axvline(BAND_FROM, color="black", lw=0.6, ls=":")
        ax.text(BAND_FROM, 3.2, " band installed", fontsize=6.5, va="top")
        ax.set_ylim(-14, 4)
        ax.yaxis.set_major_locator(MultipleLocator(4))
        ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=(1, 8, 15, 22)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.set_ylabel(r"$\Delta P_{mp}$ [%]")
        ax.set_xlabel("Date (2026)")

        # (b) loss over the day
        ax = axes[1]
        ax.errorbar(hours["hour"], hours["median"], yerr=[hours["median"] - hours["q25"], hours["q75"] - hours["median"]],
                    ls="none", marker=MK_BAND, ms=3.5, color=COL_BAND, mec="black", mew=0.3, elinewidth=0.8,
                    capsize=1.5)
        ax.axhspan(dp - half, dp + half, color=COL_BAND, alpha=0.15, lw=0)
        ax.axhline(dp, ls="--", lw=1.1, color=COL_BAND)
        ax.axhline(plateau, ls=":", lw=1.0, color="0.3")
        ax.text(18.4, plateau + 0.8, "bypass plateau (SPICE)", fontsize=6.5, ha="right", color="0.3")
        ax.set_xlim(7.5, 18.5)
        ax.set_ylim(-30, 4)
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.set_xlabel("Hour of the day")

        # (c) cumulative energy after the band
        ax = axes[2]
        cum = daily.set_index("day")[["E_clean_Wh", "E_meas_Wh", "E_spice_Wh"]].cumsum() / 1000
        x = np.arange(1, len(cum) + 1)
        ax.plot(x, cum["E_clean_Wh"], ls="-", lw=1.1, color=COL_CLEAN, marker="o", ms=2.5, mew=0)
        ax.plot(x, cum["E_meas_Wh"], ls="none", marker=MK_BAND, ms=3.5, color=COL_BAND, mec="black", mew=0.3)
        ax.plot(x, cum["E_spice_Wh"], ls="--", lw=1.1, color=COL_BAND)
        ax.set_xlabel("Days with the band")
        ax.set_ylabel("Cumulative energy [kWh]")
        ax.set_xlim(0, len(cum) + 1)
        ax.set_ylim(bottom=0)
        tot = daily[["E_clean_Wh", "E_meas_Wh", "E_spice_Wh"]].sum()
        e_meas = 100 * (tot["E_meas_Wh"] / tot["E_clean_Wh"] - 1)
        ax.text(0.04, 0.96, f"measured {e_meas:+.1f} %\nSPICE {dp:+.1f} %", transform=ax.transAxes,
                fontsize=6.5, va="top")

        for i, ax in enumerate(axes):
            ax.text(0.97, 0.05, "abc"[i], transform=ax.transAxes, ha="right", **LETTER_KW)
            ax.tick_params(pad=1.5)
        handles = [Line2D([], [], color=COL_CLEAN, marker="o", ms=3.5, mec="black", mew=0.3, ls="none",
                          label="before the band (calibration)"),
                   Line2D([], [], color=COL_BAND, marker=MK_BAND, ms=3.5, mec="black", mew=0.3, ls="none",
                          label=r"with the band, measured"),
                   Line2D([], [], color=COL_BAND, ls="--", lw=1.1, label="SPICE"),
                   Patch(color=COL_BAND, alpha=0.15, lw=0, label="SPICE 95 % interval"),
                   Line2D([], [], color=COL_CLEAN, lw=1.1, label="expected without band")]
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, -0.06),
                   frameon=False, handletextpad=0.3, columnspacing=1.0)
        fig.tight_layout(rect=(0, 0.06, 1, 1), w_pad=0.8)
        plots.save(fig, path)


# --- main --------------------------------------------------------------------------

def run(csv=config.FIELD_CSV, out=config.FIELD_DIR, fig_dir=None):
    """Whole analysis; tables and summary.json in out, the two figures in fig_dir (default out)."""
    from types import SimpleNamespace
    a = SimpleNamespace(csv=csv, out=out)
    fig_dir = fig_dir or out
    a.out.mkdir(parents=True, exist_ok=True)

    p = pair(load(a.csv))
    cal = p[(p["phase"] == "before") & p["stable"]]
    aft = p[(p["phase"] == "after") & p["stable"]]
    r0 = float(cal["r"].median())
    dp, half, bias = spice_loss()
    p["dP_%"] = 100 * (p["r"] / r0 - 1)
    p.to_csv(a.out / "pairs.csv", index=False)

    daily = daily_energy(p, r0, dp)
    daily.to_csv(a.out / "daily.csv", index=False)
    bins = loss_vs_irradiance(p, r0)
    bins.to_csv(a.out / "loss_vs_G.csv", index=False)
    hours = loss_by_hour(p, r0)
    hours.to_csv(a.out / "loss_by_hour.csv", index=False)
    plateau = bypass_plateau()

    tot = daily[["E_clean_Wh", "E_meas_Wh", "E_spice_Wh"]].sum()
    power_loss = 100 * (aft["r"].median() / r0 - 1)
    energy_loss = 100 * (tot["E_meas_Wh"] / tot["E_clean_Wh"] - 1)
    summary = {
        "data": str(a.csv), "band": {"h_mm": H_MM, "tau_I": TAU_I, "from": str(BAND_FROM), "to": str(BAND_TO)},
        "calibration": {"r0": r0, "r0_iqr": [float(cal["r"].quantile(0.25)), float(cal["r"].quantile(0.75))],
                        "n_pairs": int(len(cal)), "days": [str(cal["t"].min()), str(cal["t"].max())]},
        "after": {"n_pairs": int(len(aft)), "days": int(len(daily))},
        "spice_loss_%": dp, "spice_half_width_pt": half, "spice_bias_pt": bias,
        "power_loss_%": float(power_loss),
        "power_loss_iqr_%": [float(100 * (aft["r"].quantile(q) / r0 - 1)) for q in (0.25, 0.75)],
        "energy_loss_%": float(energy_loss),
        "energy_lost_Wh": float(tot["E_clean_Wh"] - tot["E_meas_Wh"]),
        "energy_lost_Wh_per_day": float((tot["E_clean_Wh"] - tot["E_meas_Wh"]) / len(daily)),
        "energy_meas_minus_spice_Wh": float(tot["E_meas_Wh"] - tot["E_spice_Wh"]),
        "gap_energy_pt": float(energy_loss - dp),
        "n_pairs_clean_shaded": int(p["clean_shaded"].sum()),
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2))
    figure_energy(p, r0, dp, daily, fig_dir / "fig_field_test")
    figure(p, r0, dp, half, daily, hours, plateau, fig_dir / "fig_field_test_loss")

    print(f"\ncalibration  r0 = {r0:.4f} ({len(cal)} stable pairs, {cal['t'].min():%d %b} to {cal['t'].max():%d %b})")
    print(f"after band   {len(aft)} stable pairs, {len(daily)} full days")
    print(f"SPICE (STC)  dP = {dp:+.2f} %  (95 % interval +/- {half:.2f} pt)")
    print(f"measured     power {power_loss:+.2f} % (G >= {G_MIN:g} W/m2),  energy {energy_loss:+.2f} %")
    print(f"gap          energy - SPICE = {energy_loss - dp:+.2f} pt,  "
          f"{summary['energy_lost_Wh_per_day']:.0f} Wh lost per day")
    print(f"shaded pairs (clean module, evening): {int(p['clean_shaded'].sum())}")
    print(f"bypass plateau (SPICE, tau_I = 0): {plateau:+.2f} %")
    print(hours.round(2).to_string(index=False))
    return summary
