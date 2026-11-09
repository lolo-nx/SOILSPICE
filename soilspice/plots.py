"""Figures. Every function saves PNG (+ PDF) and closes the figure.

Colours are tied to an identity, never to a position in a list: a module or
a measured tau_I keeps the same colour and marker in every figure.
"""
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize, TwoSlopeNorm, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle

from . import geometry
from .analysis import REGIME_NAMES, linear_fit
from .modules import MODULE_ORDER, MODULES

STYLE = {
    "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold", "axes.labelsize": 12,
    "axes.edgecolor": "black", "axes.linewidth": 1.2, "axes.grid": False,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.major.size": 5, "ytick.major.size": 5,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 9, "legend.frameon": True, "legend.edgecolor": "black",
    "legend.framealpha": 1.0, "legend.fancybox": False,
}

MODULE_COLOR = {"YINGLI": "#009E73", "SUNPOWER": "#D55E00", "JASOLAR": "#E69F00", "LONGI": "#0072B2",
                "MATRIX": "#8e44ad"}
MODULE_MARKER = {"LONGI": "o", "SUNPOWER": "v", "YINGLI": "D", "JASOLAR": "^", "MATRIX": "s"}

# measured tau_I values: 0.707/0.71/0.713 -> 0.70, 0.539/0.545/0.549 -> 0.53, ...
# ordered like the band: light = clean, dark = opaque (ColorBrewer Blues, colour-blind safe and readable in grey)
TAU_ANCHORS = {1.00: ("#E3B34F", "o"), 0.70: ("#D08A2E", "v"), 0.53: ("#B3611B", "D"),
               0.45: ("#8A3F12", "^"), 0.00: ("#4E2209", "s")}
TAU_TOL = 0.025

LOSS_LABEL = {"Pmp": r"$(P_{mp}-P_{mp,clean})/P_{mp,clean}$ [%]",
              "Vmp": r"$(V_{mp}-V_{mp,clean})/V_{mp,clean}$ [%]",
              "Imp": r"$(I_{mp}-I_{mp,clean})/I_{mp,clean}$ [%]"}
PARAM_LATEX = {"Pmp": r"$P_{mp}$", "Vmp": r"$V_{mp}$", "Imp": r"$I_{mp}$",
               "Voc": r"$V_{oc}$", "Isc": r"$I_{sc}$", "FF": "FF"}

EXP_COLOR, SIM_COLOR = "#2E75B6", "#D9A521"
WIRE = "#c0392b"


# --- helpers ---------------------------------------------------------------------------

def save(fig, path, pdf=True, dpi=300, *, bbox_inches="tight"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches=bbox_inches)
    if pdf:
        fig.savefig(path.with_suffix(".pdf"), bbox_inches=bbox_inches)
    plt.close(fig)
    print(f"  saved {path.with_suffix('.png')}")


def tau_style(tau):
    """(colour, marker) of a measured tau_I; values outside the anchors get a viridis colour."""
    anchor = min(TAU_ANCHORS, key=lambda a: abs(a - tau))
    if abs(anchor - tau) <= TAU_TOL:
        return TAU_ANCHORS[anchor]
    return plt.get_cmap("viridis")(float(np.clip(tau, 0, 1))), "P"


def module_label(name, with_tech=True):
    m = MODULES.get(name)
    if m is None:
        return name
    return m.label if with_tech else m.pretty_name


def sort_modules(names):
    names = list(dict.fromkeys(names))
    return sorted(names, key=lambda n: MODULE_ORDER.index(n) if n in MODULE_ORDER else len(MODULE_ORDER))


def _symmetric_ylim(ax, values, pad=0.15):
    values = np.asarray([v for v in values if np.isfinite(v)])
    if values.size:
        lo, hi = min(0.0, values.min()), max(0.0, values.max())
        margin = pad * max(abs(lo), abs(hi), 1.0)
        ax.set_ylim(lo - margin, hi + margin)


# --- module drawing ------------------------------------------------------------------

def _diode(ax, x, y, size, horizontal=True):
    s = size / 2
    if horizontal:
        ax.add_patch(Polygon([(x - s, y - s), (x - s, y + s), (x + s, y)], color=WIRE, zorder=8))
        ax.plot([x + s, x + s], [y - s, y + s], color=WIRE, lw=1.6, zorder=8)
    else:
        ax.add_patch(Polygon([(x - s, y + s), (x + s, y + s), (x, y - s)], color=WIRE, zorder=8))
        ax.plot([x - s, x + s], [y - s, y - s], color=WIRE, lw=1.6, zorder=8)


def module_schematic(module, fcov, shape, orientation, height_mm, path, apex="right"):
    """Cells with their f_cov, strings, bypass diodes and the band outline."""
    lay = geometry.layout(module, orientation)
    W, H = lay.width, lay.height
    cx, cy = (lay.x0 + lay.x1) / 2, (lay.y0 + lay.y1) / 2
    cw, ch = lay.x1[0] - lay.x0[0], lay.y1[0] - lay.y0[0]
    # strings run along the portrait columns, or along the portrait rows when the
    # bypass groups are rows (matrix shingle); in landscape everything turns
    vertical = (orientation == "portrait") == (module.bypass_axis == "cols")
    n_strings = lay.n_cols if vertical else lay.n_rows

    with plt.rc_context(STYLE):
        fig_w = 7.5 if W / H < 2.5 else 16.0               # very elongated modules (matrix shingle)
        fig, ax = plt.subplots(figsize=(fig_w, max(fig_w * H / W, 3.0)))
        cell_pt = min(cw, ch) / (1.2 * W) * fig_w * 72     # smallest cell side in points
        label_fs = min(8.5, 0.55 * cell_pt)
        for i in range(lay.n_rows):
            for j in range(lay.n_cols):
                f = fcov[i, j]
                if np.isnan(f):                    # inactive cell (matrix contact row)
                    ax.add_patch(Rectangle((lay.x0[j], lay.y0[i]), cw, ch, facecolor="0.85", edgecolor="black",
                                           hatch="////", lw=1.1, zorder=3))
                    continue
                face = to_rgba("#d9a521", 0.28) if f > 1e-3 else "white"
                ax.add_patch(Rectangle((lay.x0[j], lay.y0[i]), cw, ch, facecolor=face, edgecolor="black",
                                       lw=1.1, zorder=3))
                if f > 1e-3 and label_fs >= 4.5:
                    ax.text(cx[j], cy[i], f"{f:.2g}", ha="center", va="center", fontsize=label_fs, zorder=5,
                            rotation=90 if cw < ch else 0,
                            path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])

        dot = 0.12 * min(cw, ch)
        for s in range(n_strings):
            xs, ys = (np.full(lay.n_rows, cx[s]), cy) if vertical else (cx, np.full(lay.n_cols, cy[s]))
            ax.plot(xs, ys, color=WIRE, lw=2.0, zorder=4)
            for x, y in zip(xs, ys):
                ax.add_patch(Rectangle((x - dot / 2, y - dot / 2), dot, dot, color=WIRE, zorder=6))

        # bypass diodes: one per group of strings, drawn at mid-length
        mid = (lay.n_rows if vertical else lay.n_cols) // 2
        for first, last in module.bypass_groups:
            if vertical:
                y = (lay.y1[mid - 1] + lay.y0[mid]) / 2
                ax.plot([cx[first], cx[last]], [y, y], color=WIRE, lw=2.0, zorder=4)
                _diode(ax, (cx[first] + cx[last]) / 2, y, 0.5 * min(cw, ch), horizontal=True)
            else:
                x = (lay.x1[mid - 1] + lay.x0[mid]) / 2
                ax.plot([x, x], [cy[first], cy[last]], color=WIRE, lw=2.0, zorder=4)
                _diode(ax, x, (cy[first] + cy[last]) / 2, 0.5 * min(cw, ch), horizontal=False)

        ax.add_patch(Rectangle((0, 0), W, H, facecolor="none", edgecolor="black", lw=1.6, zorder=2))
        if height_mm:
            h = min(height_mm, H)
            band = ([(0, 0), (W, 0), (W, h), (0, h)] if shape == "rectangle"
                    else [(0, 0), (W, 0), (W, h)] if apex == "right" else [(0, 0), (W, 0), (0, h)])
            ax.add_patch(Polygon(band, fill=False, edgecolor="0.2", lw=1.8, ls=(0, (5, 3)), zorder=7))

        ax.set_xlim(-0.1 * W, 1.1 * W)
        ax.set_ylim(-0.14 * H, 1.14 * H)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{module.label} - {shape}, {orientation} "
                     f"(A_cov = {geometry.acov(fcov):.1%}, h = {height_mm:.0f} mm)")
        save(fig, path, pdf=False, dpi=200)


# --- validation figures ------------------------------------------------------------------

def iv_comparison(result, path):
    """Measured points and simulated curves (I-V and P-V) of one configuration."""
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        for m in result.measured:
            color, marker = tau_style(m.tau_i)
            sim = result.simulated.get(m.tau_i)
            for ax, y_meas in zip(axes, (m.i, m.v * m.i)):
                ax.plot(m.v, y_meas, ls="none", marker=marker, ms=6, mew=0, markevery=0.05, color=color,
                        label=rf"$\tau_I$={m.tau_i:.2f} (EXP)")
            if sim is not None:
                for ax, y_sim in zip(axes, (sim.i, sim.v * sim.i)):
                    ax.plot(sim.v, y_sim, ls="--", lw=2.2, color=color, label=rf"$\tau_I$={m.tau_i:.2f} (SPICE)")
        for ax, ylabel, tag in zip(axes, ("Current [A]", "Power [W]"), ("(a)", "(b)")):
            ax.set_xlabel("Voltage [V]")
            ax.set_ylabel(ylabel)
            ax.set_ylim(bottom=0)
            ax.legend(ncol=2, edgecolor="none", framealpha=0.9)
            ax.annotate(tag, xy=(0.5, -0.16), xycoords="axes fraction", ha="center")
        fig.suptitle(f"{result.module.label} - {result.shape}, {result.orientation} "
                     f"(A_cov = {result.acov:.1%}) - points: measured, dashed: SPICE")
        fig.tight_layout()
        save(fig, path, dpi=150)


def loss_bars(losses, title, path):
    """Measured vs simulated Pmp/Vmp/Imp losses, one group of bars per tau_I."""
    losses = losses.sort_values("tau_I", ascending=False).reset_index(drop=True)
    x, w = np.arange(len(losses)), 0.34
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), sharex=True)
        for ax, metric, tag in zip(axes, ("Pmp", "Vmp", "Imp"), ("(a)", "(b)", "(c)")):
            exp = losses[f"{metric}_loss_exp_%"].to_numpy(dtype=float)
            sim = losses.get(f"{metric}_loss_sim_%", pd.Series(np.nan, index=losses.index)).to_numpy(dtype=float)
            bars = [ax.bar(x - w / 2, exp, w, color=EXP_COLOR, label="EXP", zorder=3),
                    ax.bar(x + w / 2, sim, w, color=SIM_COLOR, hatch="///", edgecolor="white", lw=0.6,
                           label="SPICE", zorder=3)]
            for group in bars:
                for b in group:
                    h = b.get_height()
                    if np.isfinite(h):
                        ax.annotate(f"{h:.1f}", (b.get_x() + b.get_width() / 2, h), xytext=(0, 3 if h >= 0 else -3),
                                    textcoords="offset points", ha="center", va="bottom" if h >= 0 else "top",
                                    fontsize=8.5)
            _symmetric_ylim(ax, np.concatenate([exp, sim]), pad=0.18)
            ax.axhline(0, color="black", lw=1.0)
            ax.set_ylabel(LOSS_LABEL[metric])
            ax.set_xticks(x)
            ax.set_xticklabels([rf"$\tau_I$={t:.2f}" for t in losses["tau_I"]])
            ax.legend(loc="upper left", edgecolor="none", framealpha=0.9)
            ax.annotate(tag, xy=(0.5, -0.20), xycoords="axes fraction", ha="center", fontweight="bold")
        fig.suptitle(title)
        fig.tight_layout(rect=[0, 0.03, 1, 0.92])
        save(fig, path)


def parity(losses, title, path, tol_pct=5.0):
    """Simulated vs measured loss, with the y = x line and a +/- tol band."""
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(15.6, 5))
        for ax, metric in zip(axes, ("Pmp", "Vmp", "Imp")):
            x = losses[f"{metric}_loss_exp_%"].to_numpy(dtype=float)
            y = losses[f"{metric}_loss_sim_%"].to_numpy(dtype=float)
            ok = np.isfinite(x) & np.isfinite(y)
            x, y = x[ok], y[ok]
            if not len(x):
                ax.set_axis_off()
                continue
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            pad = 0.1 * ((hi - lo) or 1.0)
            d = np.linspace(lo - pad, hi + pad, 50)
            ax.fill_between(d, d - tol_pct, d + tol_pct, color="gray", alpha=0.15, label=f"±{tol_pct:.0f} pts")
            ax.plot(d, d, "k--", lw=1.2, label="y = x")
            for module, sub in losses[ok].groupby("module"):
                ax.scatter(sub[f"{metric}_loss_exp_%"], sub[f"{metric}_loss_sim_%"], s=50,
                           color=MODULE_COLOR.get(module, "0.4"), marker=MODULE_MARKER.get(module, "o"),
                           edgecolor="black", zorder=3, label=module_label(module, False))
            ss_tot = np.sum((x - x.mean()) ** 2)
            r2 = 1 - np.sum((y - x) ** 2) / ss_tot if ss_tot > 0 else np.nan
            ax.text(0.05, 0.95, f"RMSE={np.sqrt(np.mean((y - x) ** 2)):.2f} pts\nMAE={np.mean(np.abs(y - x)):.2f} pts\n"
                                f"R²={r2:.3f}\nn={len(x)}", transform=ax.transAxes, va="top", fontsize=9,
                    bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"))
            ax.set(xlim=(d[0], d[-1]), ylim=(d[0], d[-1]), aspect="equal",
                   xlabel=f"{PARAM_LATEX[metric]} loss - measured [%]", ylabel=f"{PARAM_LATEX[metric]} loss - SPICE [%]")
            ax.legend(fontsize=8, loc="lower right")
        fig.suptitle(title)
        fig.tight_layout()
        save(fig, path)


def error_bars(param_errors, title, path, params=("Pmp", "Vmp", "Imp", "Voc", "Isc")):
    """Relative error SPICE vs measurement, one panel per configuration,
    one colour per tau_I."""
    sub = param_errors[param_errors["parameter"].isin(params)]
    configs = sorted(sub[["shape", "orientation"]].drop_duplicates().itertuples(index=False, name=None))
    taus = sorted(sub["tau_I"].unique(), reverse=True)
    n_rows = int(np.ceil(len(configs) / 2)) or 1
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(n_rows, 2, figsize=(13, 4.3 * n_rows), squeeze=False)
        w = 0.82 / max(len(taus), 1)
        for ax, (shape, orientation) in zip(axes.flat, configs):
            cfg = sub[(sub["shape"] == shape) & (sub["orientation"] == orientation)]
            for k, tau in enumerate(taus):
                vals = [cfg[(cfg["parameter"] == p) & (cfg["tau_I"] == tau)]["rel_error_%"].mean() for p in params]
                ax.bar(np.arange(len(params)) + (k - (len(taus) - 1) / 2) * w, vals, w,
                       color=tau_style(tau)[0], edgecolor="black", lw=0.4, label=rf"$\tau_I$={tau:.2f}")
            ax.axhline(0, color="black", lw=0.8)
            ax.set_xticks(range(len(params)))
            ax.set_xticklabels([PARAM_LATEX[p] for p in params])
            ax.set_ylabel(r"$X_{rel}$ [%]")
            ax.set_title(f"{shape} / {orientation}")
            ax.grid(axis="y", alpha=0.3)
        for ax in list(axes.flat)[len(configs):]:
            ax.set_axis_off()
        handles, labels = axes[0][0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=min(len(taus), 8), bbox_to_anchor=(0.5, 0.96),
                   edgecolor="none")
        fig.suptitle(title, y=1.0)
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        save(fig, path)


def power_ratio(df, title, path):
    styles = {("rectangle", "landscape"): ("#2E75B6", "-", "o", True),
              ("rectangle", "portrait"): ("#2E75B6", "--", "o", False),
              ("triangle", "landscape"): ("#E63946", "-", "s", True),
              ("triangle", "portrait"): ("#E63946", "--", "s", False)}
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(100, color="0.6", lw=1.0, ls=":")
        for (shape, orientation), (color, ls, marker, filled) in styles.items():
            sub = df[(df["shape"] == shape) & (df["orientation"] == orientation)].sort_values("tau_I")
            if len(sub):
                ax.plot(sub["tau_I"], sub["power_ratio_%"], color=color, ls=ls, marker=marker, ms=7, lw=2,
                        mfc=color if filled else "white", label=f"{shape} ({orientation})")
        ax.set_xlabel(r"$\tau_I$")
        ax.set_ylabel(r"$P_{max}/P_{clean}$ [%]")
        ax.legend(loc="lower right")
        fig.suptitle(title)
        fig.tight_layout()
        save(fig, path)


# --- comparing modules -------------------------------------------------------------------

def losses_vs_x(df, x_col, x_label, title, path, show_slope=True):
    """Pmp/Vmp/Imp losses of several modules against x (tau_I or tau_eff).
    Solid line + filled marker: measured (or SPICE when there is nothing else).
    Dashed line + open marker: SPICE, when *_loss_sim_% columns are filled.
    The legend box of each panel gives gamma = d(loss)/dx per 0.1 of x."""
    df = df.dropna(subset=[x_col])
    if df.empty:
        return
    sim_cols = [f"{m}_loss_sim_%" for m in ("Pmp", "Vmp", "Imp")]
    has_sim = all(c in df for c in sim_cols) and df[sim_cols].notna().any().any()
    names = sort_modules(df["module"])

    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
        for ax, metric in zip(axes, ("Pmp", "Vmp", "Imp")):
            values, slopes = [], []
            for name in names:
                sub = df[df["module"] == name].sort_values(x_col)
                color, marker = MODULE_COLOR.get(name, "0.4"), MODULE_MARKER.get(name, "o")
                x = sub[x_col].to_numpy(dtype=float)
                gammas = []
                for col, ls, face in ((f"{metric}_loss_%", "-", color), (f"{metric}_loss_sim_%", "--", "white")):
                    if col not in sub or (ls == "--" and not has_sim):
                        gammas.append(np.nan)
                        continue
                    y = pd.to_numeric(sub[col], errors="coerce").to_numpy(dtype=float)
                    ok = np.isfinite(x) & np.isfinite(y)
                    if not ok.any():
                        gammas.append(np.nan)
                        continue
                    values.extend(y[ok])
                    ax.plot(x[ok], y[ok], color=color, ls=ls, marker=marker, ms=7, mec=color, mfc=face, lw=1.6, zorder=3)
                    a, b, _r2 = linear_fit(x[ok], y[ok])
                    gammas.append(a)
                    if show_slope and ls == "-" and np.isfinite(a):
                        xf = np.array([x[ok].min(), x[ok].max()])
                        ax.plot(xf, a * xf + b, color=color, lw=1.0, ls=":", alpha=0.55, zorder=2)
                if np.isfinite(gammas).any():
                    slopes.append((name, color, *gammas))
            _symmetric_ylim(ax, values)
            ax.axhline(0, color="black", lw=1.0, zorder=1)
            ax.set_xlabel(x_label)
            ax.set_ylabel(LOSS_LABEL[metric])
            if show_slope and slopes:
                texts = []
                for name, color, g_exp, g_sim in slopes:
                    t = f"{module_label(name, False)}: γ={g_exp * 0.1:+.1f}"
                    if np.isfinite(g_sim):
                        t += f" (SPICE {g_sim * 0.1:+.1f})"
                    texts.append(Line2D([], [], ls="none", label=t))
                leg = ax.legend(handles=texts, loc="upper left", handlelength=0, handletextpad=0, fontsize=10)
                for txt, (_n, color, *_g) in zip(leg.get_texts(), slopes):
                    txt.set_color(color)
                    txt.set_fontweight("bold")

        handles = [Line2D([], [], color=MODULE_COLOR.get(n, "0.4"), marker=MODULE_MARKER.get(n, "o"), ms=7, lw=1.6,
                          label=module_label(n)) for n in names]
        if has_sim:
            handles += [Line2D([], [], color="0.35", ls="-", marker="o", label="EXP"),
                        Line2D([], [], color="0.35", ls="--", marker="o", mfc="white", label="SPICE")]
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 0.92))
        fig.suptitle(title, y=0.98)
        fig.tight_layout(rect=[0, 0.03, 1, 0.85])
        save(fig, path)


def _grid(sub, value):
    """2-D array, rows = A_cov (increasing), columns = tau_I (increasing).
    Rows use the target A_cov when present: the achieved one differs by ~1e-4
    between runs of the same band (bisection)."""
    index = "target_A_cov" if "target_A_cov" in sub else "A_cov"
    table = sub.pivot_table(index=index, columns="tau_I", values=value, aggfunc="first")
    return table.to_numpy(dtype=float), list(table.index), list(table.columns)


HEATMAPS = (("Pmp_loss_%", "Pmp", "YlOrBr_r", False),
            ("Vmp_loss_%", "Vmp", "RdBu_r", True),     # Vmp can also increase -> diverging
            ("Imp_loss_%", "Imp", "PuBuGn_r", False))


def loss_map_grid(sweep, module_name, shape, path, scale=None, unc=None):
    """Loss maps of one module and shape: rows = portrait / landscape,
    columns = Pmp, Vmp, Imp. Colour limits come from `scale` (default: the
    whole sweep), so that the same colour means the same loss for every module."""
    scale = sweep if scale is None else scale
    vmax = {col: max(np.nanmax(np.abs(scale[col].to_numpy(dtype=float))), 1e-6) for col, *_ in HEATMAPS}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.0), squeeze=False, constrained_layout=True)
        images = {}
        for r, orientation in enumerate(("portrait", "landscape")):
            sub = sweep[(sweep["module"] == module_name) & (sweep["shape"] == shape)
                        & (sweep["orientation"] == orientation)]
            for c, (col, short, cmap_name, diverging) in enumerate(HEATMAPS):
                ax = axes[r][c]
                if sub.empty:
                    ax.set_axis_off()
                    continue
                grid, acovs, taus = _grid(sub, col)
                cmap = plt.get_cmap(cmap_name)
                norm = (TwoSlopeNorm(vmin=-vmax[col], vcenter=0, vmax=vmax[col]) if diverging
                        else Normalize(vmin=-vmax[col], vmax=0))
                images[c] = ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto", origin="lower")
                for (i, j), v in np.ndenumerate(grid):
                    if np.isfinite(v):
                        rgb = cmap(norm(v))[:3]
                        dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 0.55
                        ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7.5,
                                color="white" if dark else "black")
                ax.set_xticks(range(len(taus)))
                ax.set_yticks(range(len(acovs)))
                ax.set_xticklabels([f"{t:.1f}" for t in taus] if r == 1 else [])
                ax.set_yticklabels([f"{a * 100:.0f}" for a in acovs] if c == 0 else [])
                if r == 1:
                    ax.set_xlabel(r"$\tau_I$")
                if c == 0:
                    ax.set_ylabel(r"$A_{cov}$ [%]")
                    ax.annotate(orientation.capitalize(), xy=(-0.38, 0.5), xycoords="axes fraction", rotation=90,
                                ha="center", va="center", fontsize=12, fontweight="bold")
                if r == 0:
                    ax.set_title(LOSS_LABEL[short], fontsize=11)
        for c, im in images.items():
            fig.colorbar(im, ax=axes[:, c].tolist(), shrink=0.9, pad=0.02, aspect=30)
        fig.suptitle(f"{module_label(module_name)} - {shape}", fontweight="bold")
        add_uncertainty_note(fig, unc)
        save(fig, path)


REGIME_FILL = {0: "#e5effa", 1: "#EEBBDD", 2: "#ded9f0"}
REGIME_EDGE = {0: "#2a78d6", 1: "#606060", 2: "#4a3aa7"}


PAPER_ORDER = ("YINGLI", "SUNPOWER", "JASOLAR", "LONGI", "MATRIX")
ATLAS_ROWS = (("rectangle", "portrait"), ("triangle", "portrait"),
              ("rectangle", "landscape"), ("triangle", "landscape"))
ISO_LOSS_LEVELS = (-80, -60, -45, -35, -25, -15, -5)


def uncertainty_text(unc, on_m=False):
    """One line describing the SPICE model uncertainty (from uncertainty.global_model)."""
    if not unc:
        return ""
    b, s, z = unc["bias"], unc["sigma"], unc.get("z", 1.96)
    text = (f"SPICE model uncertainty on $\\Delta P_{{mp}}$ (validation, n = {unc['n']} flash tests): "
            f"bias $b$ = {b:+.2f} pt, $\\sigma$ = {s:.2f} pt; "
            f"real loss within [{-b - z * s:+.1f}, {-b + z * s:+.1f}] pts of the SPICE value (95 %)")
    if on_m:
        text += f";  on $M$: $\\pm{z * s:.1f}/S$"
    return text


def add_uncertainty_note(fig, unc, on_m=False, y=-0.01, fontsize=9):
    text = uncertainty_text(unc, on_m)
    if text:
        fig.text(0.01, y, text, ha="left", va="top", fontsize=fontsize, style="italic", color="0.25")


def _edges(values):
    v = np.asarray(values, dtype=float)
    mid = (v[1:] + v[:-1]) / 2
    return np.r_[v[0] - (mid[0] - v[0]), mid, v[-1] + (v[-1] - mid[-1])]


def regime_atlas(sweep, path, rows=ATLAS_ROWS, module_order=PAPER_ORDER, contours=True,
                 levels=ISO_LOSS_LEVELS, unc=None):
    """Electrical regime of every state with the iso-loss contours of Pmp.
    Columns = modules (Yingli, SunPower, JA Solar, LONGi), rows = rectangle
    portrait, triangle portrait, rectangle landscape, triangle landscape.
    Cell colour = regime (I current-limited, II clamped bypass, III mixed),
    black lines = dPmp/Pmp,clean [%] interpolated on the (tau_I, A_cov) grid."""
    names = [n for n in module_order if n in set(sweep["module"])]
    names += [n for n in sort_modules(sweep["module"]) if n not in names]
    rows = [r for r in rows if not sweep[(sweep["shape"] == r[0]) & (sweep["orientation"] == r[1])].empty]
    cmap = ListedColormap([REGIME_FILL[k] for k in (0, 1, 2)])
    cmap.set_bad("white")
    norm = Normalize(vmin=-0.5, vmax=2.5)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(len(rows), len(names), figsize=(3.5 * len(names) + 0.8, 2.75 * len(rows) + 0.6),
                                 squeeze=False, sharex=True, sharey=True)
        for r, (shape, orientation) in enumerate(rows):
            for c, name in enumerate(names):
                ax = axes[r][c]
                sub = sweep[(sweep["module"] == name) & (sweep["shape"] == shape) & (sweep["orientation"] == orientation)]
                if sub.empty:
                    ax.set_axis_off()
                    continue
                grid, acovs, taus = _grid(sub, "regime")
                a_pct, t = 100 * np.asarray(acovs), np.asarray(taus)
                ax.pcolormesh(_edges(t), _edges(a_pct), np.ma.masked_invalid(grid), cmap=cmap, norm=norm,
                              edgecolors="white", linewidth=0.8)
                if contours:
                    loss, _, _ = _grid(sub, "Pmp_loss_%")
                    lv = [l for l in levels if np.nanmin(loss) < l < np.nanmax(loss)]
                    if lv:
                        cs = ax.contour(t, a_pct, loss, levels=sorted(lv), colors="#111111", linewidths=1.0,
                                        linestyles="solid")
                        for txt in ax.clabel(cs, inline=True, fontsize=7.5, fmt="%d"):
                            txt.set_path_effects([pe.withStroke(linewidth=2.2, foreground="white")])
                ax.set_xlim(_edges(t)[0], _edges(t)[-1])
                ax.set_ylim(_edges(a_pct)[0], _edges(a_pct)[-1])
                ax.set_yticks(a_pct)
                ax.set_yticklabels([f"{a:.0f}" for a in a_pct], fontsize=8)
                ax.set_xticks(t[::2])
                ax.set_xticklabels([f"{v:.1f}" for v in t[::2]], fontsize=8)
                if r == len(rows) - 1:
                    ax.set_xlabel(r"$\tau_I$")
                if c == 0:
                    ax.set_ylabel(f"{shape}\n{orientation}\n" + r"$A_{cov}$ [%]", fontsize=10)
                if r == 0:
                    ax.set_title(module_label(name), fontsize=11)
        handles = [Rectangle((0, 0), 1, 1, facecolor=REGIME_FILL[k], edgecolor=REGIME_EDGE[k], lw=1.5,
                             label=REGIME_NAMES[k]) for k in (0, 1, 2)]
        if contours:
            handles.append(Line2D([], [], color="#111111", lw=1.2, label=r"$\Delta P_{mp}$ iso-loss contour [%]"))
        fig.tight_layout(rect=(0, 0.05 if unc else 0.03, 1, 1))
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=10,
                   bbox_to_anchor=(0.5, 0.012 if unc else -0.01))
        add_uncertainty_note(fig, unc, y=0.01)
        save(fig, path)


# --- resilience metrics (module comparison on the loss maps) -----------------------------

RES_MARKER = {"rectangle": "o", "triangle": "^"}
RES_LINESTYLE = {"rectangle": "-", "triangle": "--"}


def _modules_in(df, module_order=PAPER_ORDER):
    names = [n for n in module_order if n in set(df["module"])]
    return names + [n for n in sort_modules(df["module"]) if n not in names]


def _module_handles(names):
    return [Rectangle((0, 0), 1, 1, facecolor=MODULE_COLOR.get(n, "0.5"), edgecolor="black", lw=0.6,
                      label=module_label(n)) for n in names]


def _group_axis(ax, groups, shapes):
    """x axis of the grouped bars: orientation under each group, shape under each pair."""
    for g in range(1, len(groups)):
        new_shape = groups[g][0] != groups[g - 1][0]
        ax.axvline(g - 0.5, color="black" if new_shape else "0.75", lw=1.0 if new_shape else 0.7,
                   ls="-" if new_shape else "--")
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([o.capitalize() for _, o in groups])
    ax.tick_params(axis="x", length=0)
    for sh in shapes:
        idx = [g for g, (s_, _) in enumerate(groups) if s_ == sh]
        ax.text(np.mean(idx), -0.12, sh.capitalize(), transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=12, fontweight="bold")


def loss_vs_optical_loss(df, path, unc=None):
    """dPmp/Pmp,clean of every simulated state against the optical loss S,
    one panel per orientation; dashed line = uniform soiling (M = 1)."""
    names = _modules_in(df)
    smax = np.nanmax(df["S_%"])
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
        for n, (ax, orientation) in enumerate(zip(axes, ("portrait", "landscape"))):
            s = np.linspace(0, smax * 1.05, 50)
            ax.plot(s, -s, "k--", lw=0.9, zorder=1)
            ax.text(smax, -smax + 2.5, "$M = 1$", ha="right", va="bottom")
            for name in names:
                for shape in ("rectangle", "triangle"):
                    sub = df[(df["module"] == name) & (df["shape"] == shape) & (df["orientation"] == orientation)]
                    ax.scatter(sub["S_%"], sub["Pmp_loss_%"], s=40, marker=RES_MARKER[shape], lw=0.5,
                               color=MODULE_COLOR.get(name, "0.5"), edgecolor="white", zorder=3)
            ax.set_xlim(0, smax * 1.05)
            ax.set_xlabel(r"Optical loss $S = A_{cov}(1-\tau_I)$ [%]")
            ax.text(0.03, 0.04, f"({chr(97 + n)}) {orientation.capitalize()}", transform=ax.transAxes,
                    fontsize=13, fontweight="bold")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel(r"$\Delta P_{mp}/P_{mp,clean}$ [%]")
        handles = _module_handles(names)
        handles += [Line2D([], [], ls="", marker=RES_MARKER[s], ms=8, color="0.35", label=s.capitalize())
                    for s in ("rectangle", "triangle")]
        handles.append(Line2D([], [], color="k", ls="--", lw=0.9, label=r"$\Delta P = -S$"))
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.tight_layout(rect=(0, 0.02, 1, 0.94))
        add_uncertainty_note(fig, unc)
        save(fig, path)


def amplification_maps(df, shape, path, unc=None, vmax=20.0):
    """Loss amplification M on the (tau_I, A_cov) grid, rows = portrait /
    landscape, columns = modules. Hatched: bypass diode conducting
    (/// regime II, xxx regime III)."""
    names = _modules_in(df)
    cmap = plt.get_cmap("YlOrRd")
    norm = Normalize(vmin=0, vmax=vmax)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, len(names), figsize=(3.6 * len(names) + 1.2, 7.4), squeeze=False,
                                 sharex=True, sharey=True)
        im = None
        for r, orientation in enumerate(("portrait", "landscape")):
            for c, name in enumerate(names):
                ax = axes[r][c]
                sub = df[(df["module"] == name) & (df["shape"] == shape) & (df["orientation"] == orientation)]
                if sub.empty:
                    ax.set_axis_off()
                    continue
                M, acovs, taus = _grid(sub, "M")
                reg, _, _ = _grid(sub, "regime")
                im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto", origin="lower")
                for (i, j), v in np.ndenumerate(M):
                    if not np.isfinite(v):
                        continue
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8,
                            color="white" if norm(v) > 0.6 else "black")
                    if reg[i, j] > 0:
                        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=0, ec="0.2", alpha=0.5,
                                               hatch="///" if reg[i, j] == 1 else "xxx"))
                ax.set_xticks(range(0, len(taus), 2))
                ax.set_xticklabels([f"{t:.1f}" for t in taus[::2]], fontsize=8)
                ax.set_yticks(range(len(acovs)))
                ax.set_yticklabels([f"{a * 100:.0f}" for a in acovs], fontsize=8)
                if r == 1:
                    ax.set_xlabel(r"$\tau_I$")
                if c == 0:
                    ax.set_ylabel(f"{orientation.capitalize()}\n" + r"$A_{cov}$ [%]")
                if r == 0:
                    ax.set_title(module_label(name), fontsize=11)
        fig.suptitle(f"Loss amplification $M$ - {shape}", fontweight="bold", y=0.995)
        fig.legend(handles=[Rectangle((0, 0), 1, 1, fc="white", ec="0.3", hatch="////", label="Regime II (clamped bypass)"),
                            Rectangle((0, 0), 1, 1, fc="white", ec="0.3", hatch="xxxx", label="Regime III (mixed)")],
                   loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.955))
        fig.tight_layout(rect=(0, 0.02, 0.92, 0.91))
        if im is not None:
            cax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
            fig.colorbar(im, cax=cax).set_label(r"Loss amplification $M$ [–]")
        add_uncertainty_note(fig, unc, on_m=True)
        save(fig, path)


def resilience_scorecard(summary, path, unc=None):
    """Grouped bars of M0, S* and the maximum power loss, one bar per module."""
    names = _modules_in(summary)
    shapes = [s for s in ("rectangle", "triangle") if s in set(summary["shape"])]
    groups = [(s, o) for s in shapes for o in ("portrait", "landscape")]
    look = summary.set_index(["module", "shape", "orientation"])
    panels = (("M0", r"(a) Initial amplification $M_0$ [–]"),
              ("S_star_%", r"(b) Critical optical loss $S^*$ [%]"),
              ("max_loss_%", r"(c) Maximum power loss [%]"))
    width = 0.8 / len(names)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(19, 5.8))
        for ax, (col, title) in zip(axes, panels):
            values = summary[col].to_numpy(dtype=float)
            vmax = np.nanmax(values) if np.isfinite(values).any() else 1.0
            for g, key in enumerate(groups):
                for k, name in enumerate(names):
                    if (name, *key) not in look.index:
                        continue
                    v = float(look.loc[(name, *key), col])
                    x = g + (k - (len(names) - 1) / 2) * width
                    ax.bar(x, 0 if np.isnan(v) else v, width * 0.92, color=MODULE_COLOR.get(name, "0.5"),
                           edgecolor="black", lw=0.6, zorder=3)
                    ax.text(x, (0 if np.isnan(v) else v) + vmax * 0.012, "–" if np.isnan(v) else f"{v:.1f}",
                            ha="center", va="bottom", fontsize=8)
            _group_axis(ax, groups, shapes)
            ax.set_ylim(0, vmax * 1.18)
            ax.grid(axis="y", alpha=0.3, zorder=0)
            ax.set_title(title, loc="left", fontweight="bold")
        fig.legend(handles=_module_handles(names), loc="upper center", ncol=len(names), frameon=False,
                   bbox_to_anchor=(0.5, 1.02))
        fig.text(0.99, 0.01, "– : no bypass activation in the simulated range", ha="right", fontsize=9, style="italic")
        fig.tight_layout(rect=(0, 0.03, 1, 0.93), w_pad=2.0)
        add_uncertainty_note(fig, unc)
        save(fig, path)


def average_power(df, summary, s_full, path, unc=None):
    """Top: mean normalized power against S (1 % bins), one panel per
    orientation, with the ideal line 1 - S. Bottom: P_an,S on [0, s_full]."""
    from .resilience import power_curve
    names = _modules_in(df)
    shapes = [s for s in ("rectangle", "triangle") if s in set(df["shape"])]
    groups = [(s, o) for s in shapes for o in ("portrait", "landscape")]
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(15, 10.5))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.12)
        top = [fig.add_subplot(gs[0, i]) for i in range(2)]
        top[1].sharey(top[0])
        top[1].tick_params(labelleft=False)
        for n, (ax, orientation) in enumerate(zip(top, ("portrait", "landscape"))):
            s = np.linspace(0, s_full * 1.05, 50)
            ax.fill_between(s, 1 - s / 100, 0, color="0.93", zorder=0)
            ax.plot(s, 1 - s / 100, color="k", lw=1.0)
            ax.text(s_full * 0.97, 1 - s_full * 0.97 / 100 + 0.01, "Ideal", ha="right", va="bottom")
            for name in names:
                for shape in shapes:
                    sub = df[(df["module"] == name) & (df["shape"] == shape) & (df["orientation"] == orientation)]
                    if sub.empty:
                        continue
                    x, y = power_curve(sub["S_%"], sub["Pmp_loss_%"], s_max=s_full)
                    ax.plot(x, y, color=MODULE_COLOR.get(name, "0.5"), ls=RES_LINESTYLE[shape], lw=2,
                            marker=RES_MARKER[shape], ms=5, mec="white", mew=0.5)
            ax.set_xlim(0, s_full * 1.05)
            ax.set_ylim(-0.02, 1.03)
            ax.set_xlabel(r"Optical loss $S = A_{cov}(1-\tau_I)$ [%]")
            ax.text(0.03, 0.04, f"({chr(97 + n)}) {orientation.capitalize()}", transform=ax.transAxes,
                    fontsize=13, fontweight="bold")
        top[0].set_ylabel(r"Mean normalized power $\overline{P}/P_{clean}$ [–]")

        ax = fig.add_subplot(gs[1, :])
        look = {(r.module, r.shape, r.orientation): r for r in summary.itertuples()}
        width = 0.8 / len(names)
        for g, key in enumerate(groups):
            for k, name in enumerate(names):
                row = look.get((name, *key))
                if row is None:
                    continue
                x = g + (k - (len(names) - 1) / 2) * width
                ax.bar(x, row.P_an_S, width * 0.92, color=MODULE_COLOR.get(name, "0.5"), edgecolor="black", lw=0.6,
                       zorder=3)
                ax.text(x, row.P_an_S + 0.01, f"{row.P_an_S:.2f}", ha="center", va="bottom", fontsize=9)
        _group_axis(ax, groups, shapes)
        ax.axhline(1, color="k", lw=1)
        ax.text(len(groups) - 0.52, 1.01, "Ideal", ha="right", va="bottom")
        ax.set_ylim(0, 1.12)
        ax.set_ylabel(r"$P_{an,S}$ [–]")
        ax.set_title(r"(c) Average normalized power under soiling  $P_{an,S}=\int P\,dS\,/\int P_{ideal}\,dS$"
                     + f"  —  $S \\in [0,\\ {s_full:.0f}$ %]", loc="left", fontweight="bold")
        handles = _module_handles(names)
        handles += [Line2D([], [], color="0.35", ls=RES_LINESTYLE[s], marker=RES_MARKER[s], lw=2, label=s.capitalize())
                    for s in shapes]
        handles.append(Line2D([], [], color="k", lw=1, label="Ideal (uniform soiling)"))
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 0.99))
        fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.1)
        add_uncertainty_note(fig, unc, y=0.035)
        save(fig, path)


def average_power_table(summary, windows, path, unc=None, vmin=0.5):
    """P_an,S on the windows [0, S_ref]: one row per module, grouped by band
    configuration, one column per window."""
    names = _modules_in(summary)
    configs = [(s, o) for s in ("rectangle", "triangle") for o in ("portrait", "landscape")
               if not summary[(summary["shape"] == s) & (summary["orientation"] == o)].empty]
    cols = [f"P_an_S_0to{w:g}" for w in windows]
    look = summary.set_index(["module", "shape", "orientation"])
    values, labels = [], []
    for key in configs:
        for name in names:
            values.append([look.loc[(name, *key), c] if (name, *key) in look.index else np.nan for c in cols])
            labels.append(module_label(name))
    values = np.array(values, dtype=float)
    cmap = plt.get_cmap("RdYlGn")
    norm = Normalize(vmin=vmin, vmax=1.0)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(0.75 * len(windows) + 5.5, 0.42 * len(values) + 1.6))
        im = ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")
        for (i, j), v in np.ndenumerate(values):
            if np.isfinite(v):
                rgb = cmap(norm(v))[:3]
                dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 0.45
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.5, color="white" if dark else "black")
        n = len(names)
        for g in range(1, len(configs)):
            ax.axhline(g * n - 0.5, color="black", lw=1.5)
        ax.set_yticks(range(len(values)))
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xticks(range(len(windows)))
        ax.set_xticklabels([f"0–{w:g}" for w in windows])
        ax.set_xlabel(r"Integration window $[0, S_{ref}]$ [%]")
        ax.tick_params(length=0)
        for g, (shape, orientation) in enumerate(configs):
            ax.text(-0.42, g * n + (n - 1) / 2, f"{shape.capitalize()}\n{orientation.capitalize()}",
                    transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=10, fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02).set_label(r"$P_{an,S}$ [–]")
        fig.tight_layout()
        add_uncertainty_note(fig, unc, fontsize=8)
        save(fig, path)


def curve_error_summary(curves, path, module_order=("YINGLI", "SUNPOWER", "JASOLAR", "LONGI")):
    """Whole-curve agreement of all modules, one group of bars per configuration:
    left, accuracy on the P-V curve (nRMSE); right, bias on the I-V curve (nMBE).
    Each bar is the mean over the tau_I conditions of that configuration."""
    configs = [("rectangle", "landscape"), ("rectangle", "portrait"),
               ("triangle", "landscape"), ("triangle", "portrait")]
    names = [n for n in module_order if n in set(curves["module"])]
    mean = curves.groupby(["module", "shape", "orientation"])[["nRMSE_PV_%", "nMBE_IV_%"]].mean()
    x, w = np.arange(len(configs)), 0.8 / max(len(names), 1)

    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        for ax, col, title, ylabel in ((axes[0], "nRMSE_PV_%", "Accuracy — $P$–$V$ curve", "nRMSE [%]"),
                                       (axes[1], "nMBE_IV_%", "Bias — $I$–$V$ curve", "nMBE [%]")):
            for k, name in enumerate(names):
                vals = [mean[col].get((name, s, o), np.nan) for s, o in configs]
                ax.bar(x + (k - (len(names) - 1) / 2) * w, vals, w, color=MODULE_COLOR.get(name, "0.5"),
                       edgecolor="black", lw=0.6, label=module_label(name), zorder=3)
            ax.axhline(0, color="black", lw=1.0)
            ax.grid(axis="y", alpha=0.3, zorder=0)
            ax.set_xticks(x)
            ax.set_xticklabels([f"{s}\n{o}" for s, o in configs])
            ax.set_ylabel(ylabel)
            ax.set_title(title)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=len(names), frameon=False, fontsize=13,
                   bbox_to_anchor=(0.5, 1.06))
        fig.tight_layout()
        save(fig, path)


# --- model uncertainty -------------------------------------------------------------------

SHAPE_MARKER = {"rectangle": "s", "triangle": "^"}
REGIME_ORDER = [REGIME_NAMES[k] for k in (0, 1, 2)]


def _residual_lines(ax, b, sigma, z=1.96):
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(b, color="0.25", lw=1.4, label=f"bias b = {b:+.2f} pts")
    ax.axhspan(b - z * sigma, b + z * sigma, color="0.5", alpha=0.15, lw=0,
               label=f"b ± {z:.2f}σ = ±{z * sigma:.2f} pts")


def _residual_points(ax, df, x, col, jitter=0.0, rng=None):
    """Colour = module, marker = shape, filled = portrait, open = landscape."""
    for (module, shape, orientation), sub in df.groupby(["module", "shape", "orientation"]):
        xs = sub[x].to_numpy(dtype=float)
        if jitter:
            xs = xs + rng.uniform(-jitter, jitter, len(xs))
        color = MODULE_COLOR.get(module, "0.4")
        ax.scatter(xs, sub[col], s=46, marker=SHAPE_MARKER.get(shape, "o"), linewidths=1.2, zorder=3,
                   facecolor=color if orientation == "portrait" else "white", edgecolor=color)


def _residual_legend(fig, modules):
    handles = [Line2D([], [], ls="", marker="o", ms=8, color=MODULE_COLOR.get(m, "0.4"), label=module_label(m, False))
               for m in modules]
    handles += [Line2D([], [], ls="", marker=SHAPE_MARKER[s], ms=8, mfc="0.3", mec="0.3", label=s)
                for s in ("rectangle", "triangle")]
    handles += [Line2D([], [], ls="", marker="o", ms=8, mfc="0.3", mec="0.3", label="portrait (filled)"),
                Line2D([], [], ls="", marker="o", ms=8, mfc="white", mec="0.3", label="landscape (open)")]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.03))


def residual_diagnostics(res, metric, path, z=1.96):
    """Residuals e = loss_SPICE - loss_EXP against the continuous variables
    (top row) and per group (bottom row), with the global bias and the
    b +/- z sigma band. Used to decide whether one global sigma is enough."""
    col = f"e_{metric}"
    e = res[col].to_numpy(dtype=float)
    b, sigma = np.nanmean(e), np.nanstd(e, ddof=1)
    rng = np.random.default_rng(0)
    modules = sort_modules(res["module"])
    ylabel = f"$e$ = {PARAM_LATEX[metric]} loss SPICE − EXP [pts]"

    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), sharey=True)
        for ax, (x, xlabel) in zip(axes[0], (("A_cov", r"$A_{cov}$ [%]"), ("tau_I", r"$\tau_I$"),
                                             ("tau_eff", r"$\tau_{eff}$ (most shaded cell / sub-string)"))):
            df = res.assign(A_cov=res["A_cov"] * 100) if x == "A_cov" else res
            _residual_lines(ax, b, sigma, z)
            _residual_points(ax, df, x, col, jitter=0.012 if x == "tau_I" else 0.0, rng=rng)
            ax.set_xlabel(xlabel)
            ax.grid(alpha=0.3, zorder=0)

        groups = (("module", modules, lambda m: module_label(m, False)),
                  ("config", [f"{s}\n{o}" for s in ("rectangle", "triangle") for o in ("portrait", "landscape")], str),
                  ("regime_name", [r for r in REGIME_ORDER if r in set(res["regime_name"])], str))
        res = res.assign(config=res["shape"] + "\n" + res["orientation"])
        for ax, (factor, levels, fmt) in zip(axes[1], groups):
            levels = [lv for lv in levels if lv in set(res[factor])]
            pos = {lv: k for k, lv in enumerate(levels)}
            _residual_lines(ax, b, sigma, z)
            _residual_points(ax, res.assign(_x=res[factor].map(pos)), "_x", col, jitter=0.18, rng=rng)
            for lv, k in pos.items():
                sub = res.loc[res[factor] == lv, col].to_numpy(dtype=float)
                if len(sub) > 1:
                    m, s = sub.mean(), sub.std(ddof=1)
                    ax.errorbar(k + 0.33, m, yerr=z * s, fmt="_", ms=14, color="black", capsize=4, lw=1.4, zorder=4)
                ax.text(k, 0.02, f"n={len(sub)}", transform=ax.get_xaxis_transform(), ha="center", fontsize=8,
                        color="0.3")
            ax.set_xticks(range(len(levels)))
            ax.set_xticklabels([fmt(lv) for lv in levels], fontsize=9)
            ax.set_xlim(-0.6, len(levels) - 0.4)
            ax.grid(axis="y", alpha=0.3, zorder=0)
        for ax in axes[:, 0]:
            ax.set_ylabel(ylabel)
        axes[1, 2].text(0.98, 0.98, "black: group mean ± 1.96σ", transform=axes[1, 2].transAxes,
                        ha="right", va="top", fontsize=8)
        axes[0, 0].legend(loc="lower left", fontsize=8)
        _residual_legend(fig, modules)
        fig.tight_layout()
        save(fig, path)


def loss_with_uncertainty(sweep, metric, module_name, shape, path, acov_targets=None):
    """Corrected loss and 95 % interval against tau_I, one panel per
    orientation, one line per A_cov target (output of uncertainty.apply)."""
    sub = sweep[(sweep["module"] == module_name) & (sweep["shape"] == shape)]
    if sub.empty:
        return
    targets = sorted(acov_targets or sub["target_A_cov"].unique())
    cmap = plt.get_cmap("viridis")
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
        for ax, orientation in zip(axes, ("portrait", "landscape")):
            for k, target in enumerate(targets):
                d = sub[(sub["orientation"] == orientation) & (sub["target_A_cov"] == target)].sort_values("tau_I")
                if d.empty:
                    continue
                color = cmap(k / max(len(targets) - 1, 1) * 0.9)
                ax.fill_between(d["tau_I"], d[f"{metric}_loss_lo_%"], d[f"{metric}_loss_hi_%"], color=color,
                                alpha=0.18, lw=0)
                ax.plot(d["tau_I"], d[f"{metric}_loss_corr_%"], "-o", ms=4, lw=2, color=color,
                        label=fr"$A_{{cov}}$ = {target:.0%}")
                ax.plot(d["tau_I"], d[f"{metric}_loss_%"], ":", lw=1, color=color)
            ax.axhline(0, color="black", lw=0.8)
            ax.set(xlabel=r"$\tau_I$", title=f"{shape}, {orientation}")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel(f"{PARAM_LATEX[metric]} loss [%]")
        axes[1].legend(fontsize=8, loc="lower right", title="solid: SPICE − b, band: 95 %\ndotted: raw SPICE",
                       title_fontsize=8)
        fig.suptitle(module_label(module_name))
        fig.tight_layout()
        save(fig, path)


def residual_overview(res, metric, path, z=1.96, module_order=("YINGLI", "SUNPOWER", "JASOLAR", "LONGI")):
    """Global agreement of the losses, all flashes together.
    Left: SPICE vs measured loss with y = x, the bias and the b +/- z sigma band.
    Right: distribution of the residuals e = SPICE - EXP with the normal law
    N(b, sigma) and the statistics MBE, MAE, RMSE, sigma."""
    col = f"e_{metric}"
    x = res[f"{metric}_loss_exp_%"].to_numpy(dtype=float)
    y = res[f"{metric}_loss_sim_%"].to_numpy(dtype=float)
    e = res[col].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(e)
    x, y, e = x[ok], y[ok], e[ok]
    b, sigma = e.mean(), e.std(ddof=1)
    mae, rmse = np.abs(e).mean(), np.sqrt(np.mean(e ** 2))
    names = [n for n in module_order if n in set(res["module"])]
    names += [n for n in sort_modules(res["module"]) if n not in names]

    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [1, 1.15]})

        # parity
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max(), 0.0)
        pad = 0.06 * (hi - lo)
        d = np.array([lo - pad, hi + pad])
        ax1.fill_between(d, d + b - z * sigma, d + b + z * sigma, color="0.5", alpha=0.15, lw=0,
                         label=f"$b \\pm {z:.2f}\\sigma$")
        ax1.plot(d, d, "k--", lw=1.1, label="$y = x$")
        ax1.plot(d, d + b, color="0.25", lw=1.4, label=f"$y = x + b$")
        for name in names:
            sub = res[res["module"] == name]
            ax1.scatter(sub[f"{metric}_loss_exp_%"], sub[f"{metric}_loss_sim_%"], s=42,
                        color=MODULE_COLOR.get(name, "0.5"), edgecolor="black", lw=0.6, zorder=3)
        ax1.set(xlim=d, ylim=d, aspect="equal",
                xlabel=f"{PARAM_LATEX[metric]} loss, measured [%]", ylabel=f"{PARAM_LATEX[metric]} loss, SPICE [%]")
        ax1.grid(alpha=0.3, zorder=0)
        ax1.legend(loc="upper left", fontsize=9)
        ax1.set_title(f"Losses — {PARAM_LATEX[metric]}")

        # distribution of the residuals
        width = 0.5
        edges = np.arange(np.floor(e.min() / width) * width, np.ceil(e.max() / width) * width + width, width)
        bottom = np.zeros(len(edges) - 1)
        for name in names:
            counts, _ = np.histogram(e[(res["module"].to_numpy()[ok]) == name], bins=edges)
            ax2.bar(edges[:-1], counts, width=width, bottom=bottom, align="edge", color=MODULE_COLOR.get(name, "0.5"),
                    edgecolor="white", lw=1.0, zorder=3)
            bottom += counts
        ax2.axvspan(b - z * sigma, b + z * sigma, color="0.5", alpha=0.15, lw=0, zorder=0)
        ax2.axvline(0, color="black", lw=0.8)
        ax2.axvline(b, color="0.25", lw=1.6, zorder=4, label="bias $b$")
        for k, s in enumerate((-1, 1)):
            ax2.axvline(b + s * z * sigma, color="0.25", lw=1.0, ls="--", zorder=4,
                        label=f"$b \\pm {z:.2f}\\sigma$" if k == 0 else None)
        g = np.linspace(edges[0] - width, edges[-1] + width, 300)
        ax2.plot(g, len(e) * width * np.exp(-0.5 * ((g - b) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi)),
                 color="black", lw=1.4, label=f"$\\mathcal{{N}}(b, \\sigma)$")
        ax2.text(0.98, 0.97,
                 f"n = {len(e)}\nMBE (b) = {b:+.2f} pts\nMAE = {mae:.2f} pts\nRMSE = {rmse:.2f} pts\n"
                 f"σ = {sigma:.2f} pts\n{z:.2f}σ = ±{z * sigma:.2f} pts\nmax |e| = {np.abs(e).max():.2f} pts",
                 transform=ax2.transAxes, ha="right", va="top", fontsize=10,
                 bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"))
        ax2.set(xlabel=f"$e$ = {PARAM_LATEX[metric]} loss SPICE − EXP [% points]", ylabel="number of flashes")
        ax2.grid(axis="y", alpha=0.3, zorder=0)
        ax2.legend(loc="upper left", fontsize=9)
        ax2.set_title("Residuals")

        handles = [Rectangle((0, 0), 1, 1, facecolor=MODULE_COLOR.get(n, "0.5"), edgecolor="black", lw=0.6,
                             label=module_label(n)) for n in names]
        fig.legend(handles=handles, loc="upper center", ncol=len(names), frameon=False, fontsize=13,
                   bbox_to_anchor=(0.5, 1.07))
        fig.tight_layout()
        save(fig, path)


# --- equal band height -------------------------------------------------------------------

HEIGHT_CONFIGS = (("rectangle", "portrait"), ("rectangle", "landscape"),
                  ("triangle", "portrait"), ("triangle", "landscape"))


def _first_cell_mm(name, orientation):
    """Distance from the lower module edge to the first cell (margin + frame)."""
    left, _right, bottom, _top = MODULES[name].offsets
    return bottom if orientation == "portrait" else left


def loss_vs_height(df, path, value="Pmp_loss_%", taus=(0.5, 0.0), unc=None, ylabel=None):
    """Loss (or A_cov) against the band height, same h for every module.
    Rows = tau_I, columns = band configuration, one line per module."""
    names = _modules_in(df)
    taus = [t for t in taus if np.isclose(df["tau_I"], t).any()] or [None]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(len(taus), 4, figsize=(17, 3.9 * len(taus) + 1.0), sharex=True, sharey="row",
                                 squeeze=False)
        for r, tau in enumerate(taus):
            for c, (shape, orientation) in enumerate(HEIGHT_CONFIGS):
                ax = axes[r][c]
                sub = df[(df["shape"] == shape) & (df["orientation"] == orientation)]
                if tau is not None:
                    sub = sub[np.isclose(sub["tau_I"], tau)]
                ax.axhline(0, color="black", lw=0.8)
                for name in names:
                    d = sub[sub["module"] == name].sort_values("height_mm")
                    if d.empty:
                        continue
                    y = d[value] * (100 if value == "A_cov" else 1)
                    color = MODULE_COLOR.get(name, "0.5")
                    ax.plot(d["height_mm"], y, "-", color=color, lw=2, marker=MODULE_MARKER.get(name, "o"),
                            ms=5, mec="white", mew=0.5, zorder=3)
                ax.grid(alpha=0.3)
                if r == 0:
                    ax.set_title(f"{shape}, {orientation}", fontsize=11)
                if r == len(taus) - 1:
                    ax.set_xlabel("band height $h$ [mm]")
                if c == 0:
                    lab = ylabel or (r"$A_{cov}$ [%]" if value == "A_cov" else r"$\Delta P_{mp}/P_{mp,clean}$ [%]")
                    ax.set_ylabel(lab + ("" if tau is None else f"\n$\\tau_I$ = {tau:g}"))
        handles = _module_handles(names)
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 1.0))
        fig.tight_layout(rect=(0, 0.03 if unc else 0, 1, 0.94))
        if value != "A_cov":
            add_uncertainty_note(fig, unc, y=0.02)
        save(fig, path)


def height_scorecard(summary, path, tau=0.5, unc=None):
    """Grouped bars of the equal-height metrics at one tau_I: band height for
    1 % and 5 % loss, and P_an,h. A hatched bar reaching the top means the
    threshold is not reached within the simulated heights."""
    h_max = summary.attrs.get("h_max", 100.0)
    sub = summary[summary["tau_I"].astype(str) == str(tau)]
    names = _modules_in(sub)
    shapes = [s for s in ("rectangle", "triangle") if s in set(sub["shape"])]
    groups = [(s, o) for s in shapes for o in ("portrait", "landscape")]
    look = sub.set_index(["module", "shape", "orientation"])
    h_cols = [c for c in sub.columns if c.startswith("h") and c.endswith("_mm")]
    panels = [(c, f"({chr(97 + k)}) Band height for {c[1:-3]} % loss " + r"$h_{" + c[1:-3] + r"}$ [mm]")
              for k, c in enumerate(h_cols)]
    panels.append(("P_an_h_tau", f"({chr(97 + len(h_cols))}) Mean normalized power " + r"$P_{an,h}$"
                   + f", $h$ = 0–{h_max:g} mm [–]"))
    width = 0.8 / max(len(names), 1)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(panels), figsize=(6.4 * len(panels), 5.8))
        for ax, (col, title) in zip(np.atleast_1d(axes), panels):
            is_h = not col.startswith("P_an")
            top = h_max if is_h else 1.0
            for g, key in enumerate(groups):
                for k, name in enumerate(names):
                    if (name, *key) not in look.index:
                        continue
                    v = float(look.loc[(name, *key), col])
                    x = g + (k - (len(names) - 1) / 2) * width
                    beyond = is_h and not np.isfinite(v)
                    ax.bar(x, top if beyond else v, width * 0.92, color=MODULE_COLOR.get(name, "0.5"),
                           edgecolor="black", lw=0.6, zorder=3, hatch="///" if beyond else None,
                           alpha=0.45 if beyond else 1.0)
                    label = f">{h_max:g}" if beyond else (f"{v:.0f}" if is_h else f"{v:.2f}")
                    ax.text(x, (top if beyond else v) + top * 0.012, label, ha="center", va="bottom",
                            fontsize=8 if is_h else 6.5)
            _group_axis(ax, groups, shapes)
            ax.set_ylim(0, top * 1.15)
            if not is_h:
                ax.axhline(1, color="black", lw=1)
            ax.grid(axis="y", alpha=0.3, zorder=0)
            ax.set_title(title, loc="left", fontweight="bold", fontsize=11)
        fig.legend(handles=_module_handles(names), loc="upper center", ncol=len(names), frameon=False,
                   bbox_to_anchor=(0.5, 1.02))
        fig.text(0.99, 0.01, f"same band height on every module (margins included), $\\tau_I$ = {tau:g}; "
                             f"hatched: threshold not reached for $h \\leq$ {h_max:g} mm",
                 ha="right", fontsize=9, style="italic")
        fig.tight_layout(rect=(0, 0.04, 1, 0.93), w_pad=2.0)
        add_uncertainty_note(fig, unc, y=0.035)
        save(fig, path)


def average_power_height(df, summary, path, h_max=100.0, unc=None):
    """Equal-height counterpart of average_power. Top: mean normalized power
    against the optical band height h_opt = h (1 - tau_I), all h and tau_I
    together, one panel per orientation. Bottom: P_an,h on [0, h_max]."""
    from .resilience import OPT_BIN_MM, optical_height, power_curve
    names = _modules_in(df)
    shapes = [s for s in ("rectangle", "triangle") if s in set(df["shape"])]
    groups = [(s, o) for s in shapes for o in ("portrait", "landscape")]
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(15, 10.5))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.12)
        top = [fig.add_subplot(gs[0, i]) for i in range(2)]
        top[1].sharey(top[0])
        top[1].tick_params(labelleft=False)
        for n, (ax, orientation) in enumerate(zip(top, ("portrait", "landscape"))):
            ax.axhline(1, color="k", lw=1.0)
            for name in names:
                for shape in shapes:
                    sub = df[(df["module"] == name) & (df["shape"] == shape) & (df["orientation"] == orientation)]
                    if sub.empty:
                        continue
                    x, y = power_curve(optical_height(sub["height_mm"], sub["tau_I"]), sub["Pmp_loss_%"],
                                       s_max=h_max, bin_width=OPT_BIN_MM)
                    ax.plot(x, y, color=MODULE_COLOR.get(name, "0.5"), ls=RES_LINESTYLE[shape], lw=2,
                            marker=RES_MARKER[shape], ms=4, mec="white", mew=0.5)
            ax.set_xlim(0, h_max * 1.03)
            ax.set_ylim(-0.02, 1.05)
            ax.set_xlabel(r"Optical band height $h_{opt} = h\,(1-\tau_I)$ [mm]")
            ax.text(0.03, 0.04, f"({chr(97 + n)}) {orientation.capitalize()}", transform=ax.transAxes,
                    fontsize=13, fontweight="bold")
            ax.grid(alpha=0.3)
        top[0].set_ylabel(r"Mean normalized power $\overline{P}/P_{clean}$ [–]")

        ax = fig.add_subplot(gs[1, :])
        sub = summary.drop_duplicates(["module", "shape", "orientation"]).set_index(["module", "shape", "orientation"])
        width = 0.8 / len(names)
        for g, key in enumerate(groups):
            for k, name in enumerate(names):
                if (name, *key) not in sub.index:
                    continue
                v = float(sub.loc[(name, *key), "P_an_h"])
                x = g + (k - (len(names) - 1) / 2) * width
                ax.bar(x, v, width * 0.92, color=MODULE_COLOR.get(name, "0.5"), edgecolor="black", lw=0.6, zorder=3)
                ax.text(x, v + 0.01, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
        _group_axis(ax, groups, shapes)
        ax.axhline(1, color="k", lw=1)
        ax.text(len(groups) - 0.52, 1.01, "no loss", ha="right", va="bottom")
        ax.set_ylim(0, 1.12)
        ax.set_ylabel(r"$P_{an,h}$ [–]")
        ax.set_title(r"(c) Mean normalized power at equal band height  $P_{an,h}=\frac{1}{h_{max}}\int P\,dh_{opt}$"
                     + f"  —  $h_{{opt}} \\in [0,\\ {h_max:g}$ mm]", loc="left", fontweight="bold")
        handles = _module_handles(names)
        handles += [Line2D([], [], color="0.35", ls=RES_LINESTYLE[s], marker=RES_MARKER[s], lw=2, label=s.capitalize())
                    for s in shapes]
        fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 0.99))
        fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.1)
        add_uncertainty_note(fig, unc, y=0.035)
        save(fig, path)
