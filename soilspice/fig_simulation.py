"""Figures of chapter 2 (electrical resilience): regime map and atlas, initial amplification, mean power on
windows of the optical loss, power loss against the optical loss for the four modules. Inputs: the loss maps of
run_simulation in chapter 2 (resilience_states.csv, resilience_metrics.csv, same_acov_tauI_sweep.csv)."""
from fractions import Fraction

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.path import Path

from . import modules, plots, resilience
from .analysis import REGIME_DI_PCT, REGIME_DV_PCT, REGIME_NAMES
from .style import LETTER_KW, ORDER, ORIENTATIONS, STYLE

# optical_loss_grid: (orientation, shape, colour, marker)
RECTANGLE_MARKER = Path([(-1, -0.65), (1, -0.65), (1, 0.65), (-1, 0.65), (-1, -0.65)],
                        [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
OPTICAL_CONFIGS = [
    ("portrait", "rectangle", "#08519C", RECTANGLE_MARKER),
    ("portrait", "triangle", "#6BAED6", "^"),
    ("landscape", "rectangle", "#A4161A", RECTANGLE_MARKER),
    ("landscape", "triangle", "#F4A259", "^"),
]

SHAPES = ("rectangle", "triangle")
CONFIGS = [(s, o) for s in SHAPES for o in ORIENTATIONS]
SHAPE_MARKER = {"rectangle": RECTANGLE_MARKER, "triangle": "^"}
SHAPE_LS = {"rectangle": "-", "triangle": "--"}
LETTER = 7.5                                        # panel letters, as in the validation figures
EDGE = dict(edgecolor="black", lw=0.4)
REGIME_LABEL = {0: "I, current limited", 1: "II, bypass clamped", 2: "III, mixed"}
# Original pale blue and lavender; replace the salmon/orange regime with pale purple.
REGIME_FILL = {0: "#e5effa", 1: "#EEBBDD", 2: "#ded9f0"}
REGIME_EDGE = {0: "#2a78d6", 1: "#606060", 2: "#4a3aa7"}
S_LABEL = r"Optical loss $S = A_{\mathrm{cov}}(1-\tau_I)$ [%]"
LOSS_P = r"Loss of $P_{\mathrm{mp}}$ [%]"


def letter(ax, k, x=0.03, y=0.95, va="top", **kw):
    ax.text(x, y, "abcdefghijklmnop"[k], transform=ax.transAxes, va=va, **{**LETTER_KW, **kw})


def bold(text):
    return r"$\bf{" + text.replace(" ", r"\ ") + r"}$"


def module_handles():
    return [Rectangle((0, 0), 1, 1, facecolor=plots.MODULE_COLOR[n], label=plots.module_label(n), **EDGE)
            for n in ORDER]


def top_legend(fig, handles, ncol=None, y=1.0):
    fig.legend(handles=handles, loc="upper center", ncol=ncol or len(handles), bbox_to_anchor=(0.5, y),
               frameon=False, handlelength=1.4, handletextpad=0.4, columnspacing=1.2)


def two_row_legend(fig, styles, y=1.0, dy=0.055):
    """Modules on the first row, line and marker styles on the second."""
    top_legend(fig, module_handles(), y=y)
    top_legend(fig, styles, y=y - dy)


def config_ylabel(shape, orientation, unit=r"$A_{\mathrm{cov}}$ [%]"):
    return bold(shape.capitalize()) + "\n" + bold(orientation) + "\n" + unit


def sel(df, name, shape, orientation):
    return df[(df["module"] == name) & (df["shape"] == shape) & (df["orientation"] == orientation)]


def group_axis(ax, y=-0.14, fontsize=None):
    """x axis of grouped bars: orientation under each group, shape in bold under each pair."""
    ax.axvline(1.5, color="0.6", lw=0.6, ls="--")
    ax.set_xlim(-0.5, len(CONFIGS) - 0.5)
    ax.set_xticks(range(len(CONFIGS)))
    ax.set_xticklabels([o.capitalize() for _, o in CONFIGS], fontsize=fontsize)
    ax.tick_params(axis="x", length=0)
    for i, s in enumerate(SHAPES):
        ax.text(2 * i + 0.5, y, s.capitalize(), transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontweight="bold")


def bars(ax, look, col, w=0.8 / len(ORDER), missing="–"):
    vmax = 0
    for g, key in enumerate(CONFIGS):
        for k, name in enumerate(ORDER):
            v = look.get((name, *key), {}).get(col, np.nan)
            x = g + (k - (len(ORDER) - 1) / 2) * w
            if np.isfinite(v):
                ax.bar(x, v, w, color=plots.MODULE_COLOR[name], zorder=3, **EDGE)
                vmax = max(vmax, v)
            else:
                ax.text(x, 0, missing, ha="center", va="bottom", fontsize=7)
    return vmax


# --- figures -----------------------------------------------------------------------------

def regime_atlas(df, path, levels=(-80, -60, -45, -35, -25, -15, -5), *,
                 show_regime_legend=True, figsize=(7.2, 3.6), compact=False):
    cmap = ListedColormap([REGIME_FILL[k] for k in (0, 1, 2)])
    norm = Normalize(vmin=-0.5, vmax=2.5)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(len(CONFIGS), len(ORDER), figsize=figsize, sharex=True, sharey=True)
        for r, (shape, orientation) in enumerate(CONFIGS):
            for c, name in enumerate(ORDER):
                ax, sub = axes[r, c], sel(df, name, shape, orientation)
                grid, acovs, taus = plots._grid(sub, "regime")
                a, t = 100 * np.asarray(acovs), np.asarray(taus)
                ax.pcolormesh(plots._edges(t), plots._edges(a), np.ma.masked_invalid(grid), cmap=cmap, norm=norm,
                              edgecolors="white", linewidth=0.5)
                loss, _, _ = plots._grid(sub, "Pmp_loss_%")
                lv = sorted(l for l in levels if np.nanmin(loss) < l < np.nanmax(loss))
                if lv:
                    cs = ax.contour(t, a, loss, levels=lv, colors="black", linewidths=0.6, linestyles="solid")
                    for txt in ax.clabel(cs, inline=True, fontsize=5, fmt="%d"):
                        txt.set_path_effects([pe.withStroke(linewidth=1.5, foreground="white")])
                ax.set_xlim(plots._edges(t)[[0, -1]])
                ax.set_ylim(plots._edges(a)[[0, -1]])
                ax.set_yticks(a[::2])                      # 2, 6, 10: less crowded on short panels
                ax.set_xticks([0, 0.4, 0.8] if compact else [0, 0.2, 0.4, 0.6, 0.8])
                ax.tick_params(pad=1.5)                    # no panel letters on this figure
                if r == 0:
                    title = plots.module_label(name)
                    if compact:
                        ax.set_title(title.replace(" (", "\n("), fontsize=6.5)
                    else:
                        ax.set_title(title)
                if r == len(CONFIGS) - 1:
                    ax.set_xlabel(r"$\tau_I$")
                if c == 0:
                    ax.set_ylabel(config_ylabel(shape, orientation), labelpad=2)
        handles = [Rectangle((0, 0), 1, 1, facecolor=REGIME_FILL[k], edgecolor=REGIME_EDGE[k], lw=0.8,
                             label=REGIME_LABEL[k]) for k in (0, 1, 2)]
        if not show_regime_legend:
            handles = []
        handles.append(Line2D([], [], color="black", lw=0.6, label=r"Loss of $P_{\mathrm{mp}}$ [%]"))
        top_legend(fig, handles)
        fig.tight_layout(rect=(0, 0, 1, 0.955), h_pad=0.3, w_pad=0.3)
        plots.save(fig, path, bbox_inches=None if compact else "tight")


def loss_vs_S(df, path):
    smax = np.nanmax(df["S_%"])
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
        for k, (ax, orientation) in enumerate(zip(axes, ORIENTATIONS)):
            s = np.array([0, smax * 1.05])
            ax.plot(s, -s, color="black", lw=0.7, ls="--", zorder=1)
            for name in ORDER:
                for shape in SHAPES:
                    sub = sel(df, name, shape, orientation)
                    ax.plot(sub["S_%"], sub["Pmp_loss_%"], ls="none", marker=SHAPE_MARKER[shape], ms=3.0, mew=0,
                            color=plots.MODULE_COLOR[name], alpha=0.9, zorder=3)
            ax.set_xlim(0, smax * 1.05)
            ax.set_ylim(-102, 3)
            ax.set_xlabel(S_LABEL)
            ax.set_title(orientation.capitalize(), fontweight="bold")
            ax.tick_params(pad=1.5)
            letter(ax, k, y=0.05, va="bottom")
        axes[0].set_ylabel(LOSS_P, labelpad=2)
        handles = [Line2D([], [], ls="none", marker=SHAPE_MARKER[s], ms=3.5, mew=0, color="0.4",
                          label=s.capitalize()) for s in SHAPES]
        handles.append(Line2D([], [], color="black", ls="--", lw=0.7, label=r"Uniform soiling, $M = 1$"))
        two_row_legend(fig, handles, y=1.01, dy=0.07)
        fig.tight_layout(rect=(0, 0, 1, 0.87), w_pad=1.0)
        plots.save(fig, path)


def resilience_scorecard(metrics, path):
    look = {(r["module"], r["shape"], r["orientation"]): r for r in metrics.to_dict("records")}
    panels = (("M0", r"Initial amplification $M_0$ [–]"),
              ("S_star_%", r"Critical optical loss $S^*$ [%]"),
              ("max_loss_%", r"Maximum loss of $P_{\mathrm{mp}}$ [%]"))
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
        for k, (ax, (col, ylabel)) in enumerate(zip(axes, panels)):
            vmax = bars(ax, look, col)
            group_axis(ax, fontsize=6)
            ax.set_ylim(0, vmax * 1.15)
            ax.set_ylabel(ylabel, labelpad=2)
            ax.tick_params(axis="y", pad=1.5)
            letter(ax, k)
        top_legend(fig, module_handles())
        fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=1.2)
        plots.save(fig, path)


def initial_amplification(metrics, path):
    """M0 only (panel (a) of the scorecard), one column wide."""
    look = {(r["module"], r["shape"], r["orientation"]): r for r in metrics.to_dict("records")}
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(3.5, 2.6))
        vmax = bars(ax, look, "M0")
        group_axis(ax, fontsize=6.5)
        ax.set_ylim(0, vmax * 1.12)
        ax.set_ylabel(r"Initial amplification $M_0$ [–]", labelpad=2)
        ax.tick_params(axis="y", pad=1.5)
        fig.legend(handles=module_handles(), loc="upper center", ncol=2, bbox_to_anchor=(0.55, 1.0), frameon=False,
                   handlelength=1.2, handletextpad=0.4, columnspacing=1.2)
        fig.tight_layout(rect=(0, 0, 1, 0.86))
        plots.save(fig, path)


def average_power(df, metrics, path, s_full=10.0):
    look = {(r["module"], r["shape"], r["orientation"]): r for r in metrics.to_dict("records")}
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(7.2, 4.6))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1])
        top = [fig.add_subplot(gs[0, 0])]
        top.append(fig.add_subplot(gs[0, 1], sharey=top[0]))
        top[1].tick_params(labelleft=False)
        for k, (ax, orientation) in enumerate(zip(top, ORIENTATIONS)):
            s = np.array([0, s_full * 1.05])
            ax.plot(s, 1 - s / 100, color="black", lw=0.7, ls=":")
            for name in ORDER:
                for shape in SHAPES:
                    sub = sel(df, name, shape, orientation)
                    x, y = resilience.power_curve(sub["S_%"], sub["Pmp_loss_%"], s_max=s_full)
                    ax.plot(x, y, color=plots.MODULE_COLOR[name], ls=SHAPE_LS[shape], lw=1.0,
                            marker=SHAPE_MARKER[shape], ms=2.6, mew=0.3, mec="black")
            ax.set_xlim(0, s_full * 1.05)
            ax.set_ylim(0, 1.03)
            ax.set_xlabel(S_LABEL)
            ax.set_title(orientation.capitalize(), fontweight="bold")
            ax.tick_params(pad=1.5)
            letter(ax, k, y=0.05, va="bottom")
        top[0].set_ylabel(r"Mean normalized power $\overline{P}/P_{\mathrm{clean}}$ [–]", labelpad=2)

        ax = fig.add_subplot(gs[1, :])
        bars(ax, look, "P_an_S")
        group_axis(ax, y=-0.16)
        ax.axhline(1, color="black", lw=0.7, ls=":")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(r"$P_{\mathrm{an},S}$ on $S \in [0, %g]$%% [–]" % s_full, labelpad=2)
        ax.tick_params(axis="y", pad=1.5)
        letter(ax, 2, x=0.01)

        handles = [Line2D([], [], color="0.3", ls=SHAPE_LS[s], lw=1.0, marker=SHAPE_MARKER[s], ms=3, mew=0.3,
                          mec="black", label=s.capitalize()) for s in SHAPES]
        handles.append(Line2D([], [], color="black", ls=":", lw=0.7, label="Uniform soiling"))
        two_row_legend(fig, handles, y=1.01, dy=0.045)
        fig.tight_layout(rect=(0, 0, 1, 0.92), h_pad=1.2, w_pad=0.8)
        plots.save(fig, path)


def _average_power_colormap():
    """Shared original grey--gold--green scale for the average-power figures."""
    return LinearSegmentedColormap.from_list("klasen", [(0.0, "#7d7d7d"), (0.3, "#b9a35a"), (0.55, "#f8c421"),
                                                        (0.78, "#8fb12b"), (1.0, "#1d6e1d")])


def average_power_windows(metrics, path, windows=resilience.PAN_WINDOWS, vmin=0.5, order=ORDER):
    """order: modules drawn in each group (ORDER = the four tested modules; add
    "MATRIX" for the version with the matrix shingle, triangle = apex right)."""
    look = metrics.set_index(["module", "shape", "orientation"])
    cols = [f"P_an_S_0to{w:g}" for w in windows]
    values = np.array([[look.loc[(n, *key), c] if (n, *key) in look.index else np.nan for c in cols]
                       for key in CONFIGS for n in order], dtype=float)
    # grey -> gold -> green, as the power gain maps of Klasen et al. (Prog. Photovolt. 2022, Fig. 14)
    cmap = _average_power_colormap()
    norm = Normalize(vmin=vmin, vmax=1.0)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(3.4, 5.0 * len(order) / len(ORDER)))
        im = ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")
        for (i, j), v in np.ndenumerate(values):
            if np.isfinite(v):
                rgb = cmap(norm(v))[:3]
                dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 0.45
                txt = f"{v:.2f}" if v >= 0.995 else f"{v:.2f}"[1:]      # .87 instead of 0.87: fits the cell
                ax.text(j, i, txt, ha="center", va="center", fontsize=5.5, color="black")
        n = len(order)
        for g in range(1, len(CONFIGS)):
            ax.axhline(g * n - 0.5, color="black", lw=0.8)
        ax.set_yticks(range(len(values)))
        ax.set_yticklabels([plots.module_label(m, with_tech=False) for _ in CONFIGS for m in order])
        ax.set_xticks(range(len(windows)))
        ax.set_xticklabels([f"{w:g}" for w in windows])
        ax.set_xlabel(r"Window $[0, S_{\mathrm{ref}}]$, $S_{\mathrm{ref}}$ [%]")
        ax.tick_params(length=0, pad=1.5)
        for g, (shape, orientation) in enumerate(CONFIGS):
            ax.text(-0.5, g * n + (n - 1) / 2, f"{shape.capitalize()}\n{orientation}", transform=ax.get_yaxis_transform(),
                    ha="center", va="center", fontweight="bold", rotation=90)
        fig.colorbar(im, ax=ax, orientation="horizontal", location="top", fraction=0.03, pad=0.02,
                     aspect=30).set_label(r"$P_{\mathrm{an},S}$ [–]", labelpad=2)
        fig.tight_layout()
        plots.save(fig, path)


def average_power_windows_horizontal(metrics, path, windows=resilience.PAN_WINDOWS):
    """Module-grouped power table with configuration symbols and black numbers."""
    keys = [(name, shape, orientation) for name in ORDER for shape, orientation in CONFIGS]
    selected = metrics[metrics.module.isin(ORDER)]
    if len(selected) != 16 or selected.duplicated(["module", "shape", "orientation"]).any():
        raise ValueError("Expected exactly 16 module/shape/orientation configurations")
    cols = [f"P_an_S_0to{w:g}" for w in windows]
    ordered = selected.set_index(["module", "shape", "orientation"]).loc[keys, cols]
    values = ordered.to_numpy(float).T
    if not np.isfinite(values).all():
        raise ValueError("All saved average-power values must be finite")
    cmap = _average_power_colormap()
    norm = Normalize(vmin=.5, vmax=1.)
    with plt.rc_context({**STYLE, "xtick.labelsize": 6, "ytick.labelsize": 7}):
        fig = plt.figure(figsize=(183 / 25.4, 58 / 25.4))
        bottom, height = 13 / 58, 39 / 58
        ax = fig.add_axes([.09, bottom, .78, height])
        im = ax.pcolormesh(np.arange(17) - .5, np.arange(len(windows) + 1) - .5, values,
                          cmap=cmap, norm=norm,
                          shading="flat", edgecolors="none", linewidth=0, rasterized=False)
        for (i, j), value in np.ndenumerate(values):
            label = f"{value:.2f}" if value >= .995 else f"{value:.2f}"[1:]
            ax.text(j, i, label, ha="center", va="center", fontsize=7, color="black")
        ax.set(xlim=(-.5, 15.5), ylim=(len(windows) - .5, -.5))
        ax.set_xticks([])
        ax.set_yticks(range(len(windows)))
        ax.set_yticklabels([f"{w:g}" for w in windows])
        ax.set_ylabel(r"Window $[0,S_{\mathrm{ref}}]$, $S_{\mathrm{ref}}$ [%]", labelpad=6)
        ax.tick_params(length=0, pad=4)
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        for j, (_, shape, orientation) in enumerate(keys):
            ax.plot(j, -.065, transform=ax.get_xaxis_transform(), clip_on=False,
                    ls="none", marker=SHAPE_MARKER[shape], ms=4, mew=.65,
                    mec="0.25", mfc="none" if orientation == "portrait" else "0.25")
        for group, name in enumerate(ORDER):
            if group:
                ax.axvline(group * 4 - .5, color="black", lw=.8)
            ax.text(group * 4 + 1.5, 1.035, plots.module_label(name, with_tech=False),
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, fontweight="bold")
        handles = [Line2D([], [], ls="none", marker=SHAPE_MARKER[shape], ms=4, mew=.65,
                          mec="0.25", mfc="none" if orientation == "portrait" else "0.25",
                          label=f"{shape.capitalize()}, {orientation}") for shape, orientation in CONFIGS]
        legend = fig.legend(handles=handles, loc="center", ncol=4,
                            bbox_to_anchor=(.48, 4.5 / 58), bbox_transform=fig.transFigure,
                            fontsize=6.2, frameon=False, borderpad=0,
                            handletextpad=.5, columnspacing=1.7, handlelength=1.4)
        cax = fig.add_axes([.90, bottom, .016, height])
        cb = fig.colorbar(im, cax=cax, orientation="vertical", ticks=np.linspace(.5, 1, 6))
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor("face")
        cb.ax.yaxis.set_ticks_position("right")
        cb.ax.yaxis.set_label_position("right")
        cb.set_label(r"$P_{\mathrm{an},S}$ [–]", labelpad=4)
        cb.ax.tick_params(labelsize=6.5, length=2, pad=2)
        cb.outline.set_linewidth(.5)
        np.testing.assert_array_equal(np.asarray(im.get_array()).reshape(values.shape), values)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for text in [*ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels(), ax.yaxis.label,
                     cb.ax.yaxis.label, *cb.ax.get_yticklabels(), *legend.get_texts()]:
            box = text.get_window_extent(renderer)
            if box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height:
                raise ValueError(f"Text clipped in horizontal power figure: {text.get_text()}")
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path.with_suffix(".pdf"))
        fig.savefig(path.with_suffix(".png"), dpi=600)
        plt.close(fig)
    return ordered.reset_index()


def resilience_overview(metrics, path, windows=resilience.PAN_WINDOWS):
    """Saved M0 values beside the existing average-power heatmap; no metric calculation."""
    keys = [(name, shape, orientation) for shape, orientation in CONFIGS for name in ORDER]
    selected = metrics[metrics.module.isin(ORDER)].copy()
    if len(selected) != 16 or selected.duplicated(["module", "shape", "orientation"]).any():
        raise ValueError("Expected exactly four modules × two shapes × two orientations")
    lookup = selected.set_index(["module", "shape", "orientation"])
    cols = [f"P_an_S_0to{w:g}" for w in windows]
    ordered = lookup.loc[keys, ["M0", *cols]].copy()
    m0, values = ordered.M0.to_numpy(float), ordered[cols].to_numpy(float)
    if not np.isfinite(m0).all() or not np.isfinite(values).all():
        raise ValueError("All 16 saved M0 values and average-power windows must be finite")
    cmap, norm = _average_power_colormap(), Normalize(vmin=.5, vmax=1.)
    with plt.rc_context({**STYLE, "font.size": 7, "ytick.labelsize": 7, "xtick.labelsize": 7}):
        fig = plt.figure(figsize=(183 / 25.4, 145 / 25.4))
        ax_m0 = fig.add_axes([.215, .105, .16, .715])
        ax_pan = fig.add_axes([.405, .105, .57, .715], sharey=ax_m0)
        # pcolormesh retains the original colour values as vector PDF cells.
        im = ax_pan.pcolormesh(np.arange(len(windows) + 1) - .5, np.arange(17) - .5,
                              values, cmap=cmap, norm=norm, shading="flat", rasterized=False,
                              edgecolors="none", linewidth=0)
        m0_text, pan_text = [], []
        for i, value in enumerate(m0):
            m0_text.append(ax_m0.text(.5, i, f"{value:.1f}", ha="center", va="center", fontsize=7.5))
        for (i, j), value in np.ndenumerate(values):
            label = f"{value:.2f}" if value >= .995 else f"{value:.2f}"[1:]
            pan_text.append(ax_pan.text(j, i, label, ha="center", va="center", fontsize=7, color="black"))
        ax_m0.set_xlim(0, 1)
        ax_m0.set_ylim(15.5, -.5)
        ax_m0.set_xticks([])
        ax_m0.set_yticks(range(16))
        ax_m0.set_yticklabels([plots.module_label(name, with_tech=False) for name, _, _ in keys])
        ax_m0.tick_params(length=0, pad=6)
        ax_pan.tick_params(length=0, pad=3, labelleft=False)
        ax_pan.set_xlim(-.5, len(windows) - .5)
        ax_pan.set_xticks(range(len(windows)))
        ax_pan.set_xticklabels([f"{w:g}" for w in windows])
        ax_pan.set_xlabel(r"Window $[0,S_{\mathrm{ref}}]$, $S_{\mathrm{ref}}$ [%]", labelpad=7)
        for ax in (ax_m0, ax_pan):
            ax.spines["top"].set_visible(True)
            ax.spines["right"].set_visible(True)
            for group in range(1, 4):
                ax.axhline(group * 4 - .5, color="black", lw=.8)
        for spine in ax_m0.spines.values():
            spine.set_color("0.75")
        # A white gap keeps the numerical M0 column separate from the colour scale.
        for group, (shape, orientation) in enumerate(CONFIGS):
            centre = group * 4 + 1.5
            ax_m0.text(-.80, centre, f"{shape.capitalize()}\n{orientation}",
                       transform=ax_m0.get_yaxis_transform(), ha="center", va="center",
                       fontweight="bold", fontsize=7, rotation=90, clip_on=False)
        fig.text(.295, .965, "(a) Initial amplification", ha="center", va="top", fontsize=8, fontweight="bold")
        fig.text(.295, .915, r"$M_0$ [–]", ha="center", va="center", fontsize=8)
        fig.text(.69, .965, "(b) Average normalised power", ha="center", va="top", fontsize=8, fontweight="bold")
        cax = fig.add_axes([.405, .865, .57, .022])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=np.linspace(.5, 1, 6))
        cb.solids.set_rasterized(False)
        cb.set_label(r"$P_{\mathrm{an},S}$ [–]", labelpad=4)
        cb.ax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=6.5, pad=2, length=2)
        cb.outline.set_linewidth(.5)
        np.testing.assert_array_equal(np.asarray(im.get_array()).reshape(values.shape), values)
        assert [text.get_text() for text in m0_text] == [f"{value:.1f}" for value in m0]
        assert [text.get_text() for text in pan_text] == [f"{v:.2f}" if v >= .995 else f"{v:.2f}"[1:]
                                                        for v in values.flat]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for text in [*fig.texts, *ax_m0.texts, *ax_m0.get_yticklabels(), *ax_pan.get_xticklabels(),
                     ax_pan.xaxis.label, cb.ax.xaxis.label, *cb.ax.get_xticklabels()]:
            box = text.get_window_extent(renderer)
            if box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height:
                raise ValueError(f"Text clipped at final article size: {text.get_text()}")
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path.with_suffix(".pdf"))
        fig.savefig(path.with_suffix(".png"), dpi=600)
        plt.close(fig)
    ordered["M0_display"] = [f"{value:.1f}" for value in m0]
    return ordered.reset_index()


def amplification(df, shape, path, vmax=18.0):
    cmap, norm = plt.get_cmap("YlOrRd"), Normalize(vmin=0, vmax=vmax)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, len(ORDER), figsize=(7.2, 3.7), sharex=True, sharey=True)
        im = None
        for r, orientation in enumerate(ORIENTATIONS):
            for c, name in enumerate(ORDER):
                ax, sub = axes[r, c], sel(df, name, shape, orientation)
                M, acovs, taus = plots._grid(sub, "M")
                reg, _, _ = plots._grid(sub, "regime")
                a, t = 100 * np.asarray(acovs), np.asarray(taus)
                im = ax.pcolormesh(plots._edges(t), plots._edges(a), M, cmap=cmap, norm=norm,
                                   edgecolors="white", linewidth=0.5)
                dt, da = t[1] - t[0], a[1] - a[0]
                for (i, j), v in np.ndenumerate(reg):
                    if v > 0:
                        ax.add_patch(Rectangle((t[j] - dt / 2, a[i] - da / 2), dt, da, fill=False, lw=0, ec="0.25",
                                               hatch="////" if v == 1 else "xxxx"))
                ax.set_yticks(a)
                ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
                ax.tick_params(pad=1.5)
                letter(ax, r * len(ORDER) + c, x=0.04, y=0.96,
                       bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.8))
                if r == 0:
                    ax.set_title(plots.module_label(name))
                else:
                    ax.set_xlabel(r"$\tau_I$")
                if c == 0:
                    ax.set_ylabel(bold(orientation.capitalize()) + "\n" + r"$A_{\mathrm{cov}}$ [%]", labelpad=2)
        top_legend(fig, [Rectangle((0, 0), 1, 1, fc="white", ec="0.25", lw=0.4, hatch="////", label=REGIME_LABEL[1]),
                         Rectangle((0, 0), 1, 1, fc="white", ec="0.25", lw=0.4, hatch="xxxx", label=REGIME_LABEL[2])])
        fig.tight_layout(rect=(0, 0, 0.91, 0.93), h_pad=0.5, w_pad=0.3)
        cax = fig.add_axes([0.925, 0.14, 0.015, 0.72])
        fig.colorbar(im, cax=cax).set_label(r"Loss amplification $M$ [–]")
        plots.save(fig, path)


LOSS_MAPS = (("Pmp_loss_%", r"Loss of $P_{\mathrm{mp}}$ [%]", "YlOrBr_r", False),
             ("Vmp_loss_%", r"Change of $V_{\mathrm{mp}}$ [%]", "RdBu_r", True),
             ("Imp_loss_%", r"Loss of $I_{\mathrm{mp}}$ [%]", "PuBuGn_r", False))


def loss_maps(df, name, shape, path):
    vmax = {col: max(np.nanmax(np.abs(df[col].to_numpy(dtype=float))), 1e-6) for col, *_ in LOSS_MAPS}
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 3.9), sharex=True, sharey=True)
        for c, (col, label, cmap_name, diverging) in enumerate(LOSS_MAPS):
            cmap = plt.get_cmap(cmap_name)
            norm = (TwoSlopeNorm(vmin=-vmax[col], vcenter=0, vmax=vmax[col]) if diverging
                    else Normalize(vmin=-vmax[col], vmax=0))
            for r, orientation in enumerate(ORIENTATIONS):
                ax = axes[r, c]
                grid, acovs, taus = plots._grid(sel(df, name, shape, orientation), col)
                a, t = 100 * np.asarray(acovs), np.asarray(taus)
                im = ax.pcolormesh(plots._edges(t), plots._edges(a), grid, cmap=cmap, norm=norm,
                                   edgecolors="white", linewidth=0.5)
                for (i, j), v in np.ndenumerate(grid):
                    if np.isfinite(v):
                        rgb = cmap(norm(v))[:3]
                        dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 0.55
                        ax.text(t[j], a[i], f"{v:.0f}", ha="center", va="center", fontsize=5,
                                color="white" if dark else "black")
                ax.set_yticks(a)
                ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
                ax.tick_params(pad=1.5)
                ax.tick_params(length=0)
                letter(ax, r * 3 + c, x=-0.02, y=1.02, va="bottom")
                if r == 1:
                    ax.set_xlabel(r"$\tau_I$")
                if c == 0:
                    ax.set_ylabel(bold(orientation.capitalize()) + "\n" + r"$A_{\mathrm{cov}}$ [%]", labelpad=2)
            fig.colorbar(im, ax=axes[:, c].tolist(), orientation="horizontal", location="top", fraction=0.05,
                         pad=0.06, aspect=30).set_label(label, labelpad=2)
        fig.suptitle(f"{plots.module_label(name)}, {shape} band", fontweight="bold", fontsize=7)
        plots.save(fig, path)


def regime_map(df, path, dv=REGIME_DV_PCT, di=REGIME_DI_PCT, *,
               legend_right=False, legend_below=False):
    """How the regimes are defined: every simulated state in the (dVmp, dImp) plane, on the
    regions of analysis.regime (I: dVmp > dv, II: dVmp <= dv and |dImp| < di, III: the rest)."""
    with plt.rc_context(STYLE):
        figsize = (7.2, 2.6) if legend_right else (2.7, 6.5) if legend_below else (3.5, 3.1)
        fig, ax = plt.subplots(figsize=figsize)
        x0, x1, y0, y1 = -75, 22, -103, 13
        box = dict(lw=0, zorder=0)
        ax.add_patch(Rectangle((dv, y0), x1 - dv, y1 - y0, fc=REGIME_FILL[0], **box))
        ax.add_patch(Rectangle((x0, y0), dv - x0, y1 - y0, fc=REGIME_FILL[2], **box))
        ax.add_patch(Rectangle((x0, -di), dv - x0, 2 * di, fc=REGIME_FILL[1], ec=REGIME_EDGE[1],
                               lw=0.6, zorder=1))
        ax.axvline(dv, color="0.3", lw=0.6, zorder=1)
        ax.axhline(0, color="black", lw=0.5, zorder=1)
        for name in ORDER:
            for shape in SHAPES:
                for o in ORIENTATIONS:
                    sub = sel(df, name, shape, o)
                    filled = o == "landscape"
                    ax.plot(sub["Vmp_loss_%"], sub["Imp_loss_%"], ls="none", marker=SHAPE_MARKER[shape], ms=2.8,
                            mew=0 if filled else 0.5, mec=plots.MODULE_COLOR[name],
                            mfc=plots.MODULE_COLOR[name] if filled else "none", alpha=0.9, zorder=3)
        kw = dict(fontsize=7, fontweight="bold", color="0.25")
        ax.text(x1 - 2, y0 + 3, "I", ha="right", va="bottom", **kw)
        ax.text(x0 + 2, di + 1.0, "II", ha="left", va="bottom", **kw)
        ax.text(x0 + 2, y0 + 3, "III", ha="left", va="bottom", **kw)
        ax.text(dv + 1, y0 + 3,
                rf"$\Delta V_{{\mathrm{{mp}}}}/V_{{\mathrm{{mp,clean}}}} = {dv:g}$%",
                rotation=90, ha="left", va="bottom",
                fontsize=6.5, color="0.3")
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_xlabel(r"$\Delta V_{\mathrm{mp}}/V_{\mathrm{mp,clean}}$ [%]")
        ax.set_ylabel(r"$\Delta I_{\mathrm{mp}}/I_{\mathrm{mp,clean}}$ [%]", labelpad=2)
        ax.tick_params(pad=1.5)
        handles = [Line2D([], [], ls="none", marker="o", ms=3.5, mfc=plots.MODULE_COLOR[n], mew=0,
                          label=plots.module_label(n, with_tech=False)) for n in ORDER]
        # one entry per (shape, orientation): marker = band shape, open = portrait, filled = landscape
        handles += [Line2D([], [], ls="none", marker=SHAPE_MARKER[sh], ms=3.5, mec="0.4",
                           mfc="0.4" if o == "landscape" else "none", mew=0 if o == "landscape" else 0.5,
                           label=f"{sh.capitalize()}, {o}") for sh in SHAPES for o in ORIENTATIONS]
        handles += [Rectangle((0, 0), 1, 1, facecolor=REGIME_FILL[k], edgecolor=REGIME_EDGE[k], lw=0.8,
                              label=REGIME_LABEL[k]) for k in (0, 1, 2)]
        if legend_right or legend_below:
            if legend_below:
                fig.subplots_adjust(left=0.24, right=0.97, bottom=0.61, top=0.90)
                groups = ((handles[:4], "Module", 0.50, 2),
                          (handles[4:8], "Band configuration", 0.37, 1),
                          (handles[8:], "Regime", 0.19, 1))
            else:
                fig.subplots_adjust(left=0.08, right=0.49, bottom=0.18, top=0.96)
                groups = ((handles[:4], "Module", 0.98, 2),
                          (handles[4:8], "Band configuration", 0.70, 1),
                          (handles[8:], "Regime", 0.38, 1))
            for group, title, y, ncol in groups:
                fig.legend(handles=group, title=title,
                           loc="upper center" if legend_below else "upper left", ncol=ncol,
                           bbox_to_anchor=(0.5 if legend_below else 0.55, y),
                           frameon=False, borderaxespad=0,
                           fontsize=7 if legend_below else 6.5, title_fontsize=7, handlelength=1.2,
                           handletextpad=0.4, columnspacing=1.2)
        else:
            ax.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 1.02), frameon=False,
                      borderaxespad=0, handlelength=1.2, handletextpad=0.4, columnspacing=0.9)
            fig.tight_layout(pad=0.6)
        plots.save(fig, path, bbox_inches=None if legend_right or legend_below else "tight")


BIN = 1.0   # [%] bin width of S


def binned(S, loss, width=BIN):
    """Mean, min and max of the loss in bins of S, starting at (0, 0)."""
    edges = np.arange(0.0, np.nanmax(S) + width + 1e-9, width)
    idx = np.clip(np.digitize(S, edges) - 1, 0, len(edges) - 2)
    rows = [(0.0, 0.0, 0.0, 0.0)]
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.any():
            rows.append((S[sel].mean(), loss[sel].mean(), loss[sel].min(), loss[sel].max()))
    return np.array(rows).T


def optical_loss_grid(sweep, path, names=ORDER, show_pan=False):
    """2 x 2 figure with positive power losses and a shared legend on top. With show_pan,
    P_an,S of each configuration written in its panel and its definition once below
    (off by default: P_an,S is in the average_power_windows table of the paper).
    Horizontal guides mark the idealized loss of the smallest bypass substring:
    one third for the balanced modules, one quarter for SunPower. These are
    architecture reference levels, not exact diode-activation thresholds.
    The Y axis runs from 0 % at the bottom to 100 % at the top.
    """
    names = [n for n in names if (sweep["module"] == n).any()]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.6), sharex=True, sharey=True,
                                 gridspec_kw={"wspace": 0.08, "hspace": 0.28})
        fig.subplots_adjust(left=0.09, right=0.99, bottom=0.1, top=0.83)
        s_ref = None
        for k, (ax, name) in enumerate(zip(axes.flat, names)):
            mx = resilience.add_optical_loss(sweep[sweep["module"] == name])
            summ = resilience.summary(mx)
            pan = summ.set_index(["orientation", "shape"])["P_an_S"]
            s_ref = summ.attrs["s_full"]
            s_max = 10.0
            s = np.linspace(0, s_max * 1.04, 20)
            ax.plot(s, s, color="0.35", lw=0.9, ls="--", zorder=2)
            # uniform soiling line: named in the legend only (a label in the panel hits the zero line)
            lines = []
            for ori, shape, color, marker in OPTICAL_CONFIGS:
                g = mx[(mx["orientation"] == ori) & (mx["shape"] == shape)]
                if g.empty:
                    continue
                x, mean, lo, hi = binned(g["S_%"].to_numpy(float), -g["Pmp_loss_%"].to_numpy(float))
                ax.fill_between(x, lo, hi, color=color, alpha=0.18, lw=0, zorder=1)
                ax.plot(x, mean, color=color, lw=1.0, ls="-.", marker=marker, ms=4,
                        mfc="white" if ori == "portrait" else color, mec=color,
                        mew=0.9, zorder=4, label=f"{ori.capitalize()}, {shape}")
                lines.append((ori, shape, color))
            module = modules.get(name)
            if module.bypass_groups:
                total = module.n_rows if module.bypass_axis == "rows" else module.n_cols
                fraction = Fraction(min(last - first + 1 for first, last in module.bypass_groups), total)
                level = 100.0 * float(fraction)
                ax.axhline(level, color="0.55", lw=0.8, ls=":", zorder=0.8)
                ax.text(0.025, level + 2.0,
                        rf"${fraction.numerator}/{fraction.denominator}$", transform=ax.get_yaxis_transform(),
                        fontsize=6.5, color="0.35", ha="left", va="bottom", zorder=6)
            if show_pan:
                # P_an,S of the four configurations, bottom-left corner
                y0 = 0.05
                ax.text(0.04, y0 + 0.07 * len(lines), r"$P_{an,S}$", transform=ax.transAxes, fontsize=7,
                        ha="left", va="bottom")
                for j, (ori, shape, color) in enumerate(reversed(lines)):
                    ax.text(0.04, y0 + 0.07 * j, f"{ori.capitalize()}, {shape}: {pan[(ori, shape)]:.2f}",
                            transform=ax.transAxes, fontsize=7, color=color, ha="left", va="bottom",
                            fontweight="bold", zorder=6,
                            bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.6))
            ax.set_title(f"{chr(97 + k)}   {modules.get(name).label}", loc="left", fontsize=7,
                         fontweight="bold", pad=4)
            ax.set(xlim=(-0.3, s_max * 1.04), ylim=(-5, 105))
            ax.set_xticks(np.arange(0, s_max + 0.1, 2))
            ax.set_yticks(np.arange(0, 101, 20))
        fig.supxlabel(r"Optical loss $S = A_\mathrm{cov}(1-\tau_I)$ [%]",
                      x=0.54, y=0.01, fontsize=7)
        # Shared axis labels for the four panels.
        fig.supylabel("Power loss (%)", x=0.015, y=0.465, fontsize=7)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        handles.append(plt.Line2D([], [], color="0.35", lw=0.9, ls="--"))
        labels.append("Uniform soiling approximation")
        handles.append(plt.Line2D([], [], color="0.55", lw=0.8, ls=":"))
        labels.append("Bypass loss reference")
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.895), ncol=3, fontsize=7,
                   handlelength=2.4, columnspacing=1.2, borderpad=0.5)
        if show_pan:
            fig.text(0.5, 0.035,
                     r"$P_{an,S} = \dfrac{\int_0^{S_\mathrm{ref}} P/P_\mathrm{clean}\,dS}{\int_0^{S_\mathrm{ref}} (1-S)\,dS}$"
                     rf"      $P/P_\mathrm{{clean}}$: mean power in 1 % bins of $S$;  $S_\mathrm{{ref}}$ = {s_ref:.0f} %;"
                     r"  $P_{an,S}$ = 1: loss equal to the blocked light",
                     ha="center", va="top", fontsize=7)
        plots.save(fig, path, dpi=600)


