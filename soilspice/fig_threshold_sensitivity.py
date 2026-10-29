"""Additional 3 × 3 washing-threshold figures; existing paper figures are untouched."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from . import energy
from .energy_sensitivity import THRESHOLDS
from .style import MODULE_COLOR, ORDER, ORIENTATIONS, PRETTY, STYLE

REFERENCE_THRESHOLD = 15.
METRICS = ("energy_none_band_extra_loss_pct", "margin10_npv25_eur_kwp", "manual_gain_pct_clean_revenue")
Y_LABELS = ("Additional band energy loss\n[% of clean yield]",
            "Net present value of +10 mm margin\n[€ kWp$^{-1}$, 25 years]",
            "Manual-cleaning gain, total soiling\n[% of clean revenue]")
WIDTH_MM, HEIGHT_MM = 183., 170.


def limits(values, *, nonnegative=False):
    """Limits include zero, negative cleaning gains and a small readable padding."""
    lo, hi = min(0., np.min(values)), max(0., np.max(values))
    span = hi - lo
    if span == 0:
        return (0, 1) if nonnegative else (-1, 1)
    return (0 if nonnegative else lo - .08 * span, hi + .10 * span)


def draw(data, market, path):
    """Every line is a configuration. No uncertainty bars or long caption in the image."""
    q = data[data.market == market]
    if len(q) != 120 or q.duplicated(["site", "module", "orientation", "R_H_mm_event"]).any():
        raise ValueError("expected 3 sites × 8 configurations × 5 thresholds")
    for _, g in q.groupby(["site", "module", "orientation"]):
        if tuple(sorted(g.R_H_mm_event)) != THRESHOLDS:
            raise ValueError("each curve requires all five simulated thresholds")
    # Energy and cleaning: a common Neuchâtel/site scale hides its small effects.
    # Ground-market margin values remain readable on one common scale. Rooftop
    # values span a wider range, so Neuchâtel also gets its own margin scale there.
    separate_rows = (0, 2) if market == "ground" else (0, 1, 2)
    axes_limits = []
    for row, metric in enumerate(METRICS):
        common_data = q[q.site != "Neuchâtel"] if row in separate_rows else q
        common = limits(common_data[metric].to_numpy(), nonnegative=row == 0)
        local = limits(q[q.site == "Neuchâtel"][metric].to_numpy(), nonnegative=row == 0) if row in separate_rows else common
        axes_limits.append([local, common, common])
    style = {**STYLE, "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
             "ytick.labelsize": 6.5, "legend.fontsize": 7}
    with plt.rc_context(style):
        fig, axes = plt.subplots(3, 3, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), sharex=True)
        fig.subplots_adjust(left=.115, right=.985, bottom=.095, top=.855, hspace=.28, wspace=.34)
        for col, site in enumerate(energy.SITES):
            for row, metric in enumerate(METRICS):
                ax = axes[row, col]
                ax.axvline(REFERENCE_THRESHOLD, color="0.65", lw=.65, ls=":", zorder=0)
                if row > 0:
                    ax.axhline(0, color="0.3", lw=.6, zorder=1)
                for module in ORDER:
                    for orientation, ls, marker in (("portrait", "-", "o"), ("landscape", "--", "^")):
                        g = q[(q.site == site) & (q.module == module) & (q.orientation == orientation)].sort_values("R_H_mm_event")
                        ax.plot(g.R_H_mm_event, g[metric], color=MODULE_COLOR[module], ls=ls, marker=marker,
                                ms=2.7, lw=.9, mew=.35, zorder=3)
                ax.set_xlim(4, 26)
                ax.set_ylim(*axes_limits[row][col])
                ax.set_xticks(THRESHOLDS)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
                ax.tick_params(pad=1.5)
                ax.text(.025, 1.025, "abcdefghi"[row * 3 + col], transform=ax.transAxes,
                        va="bottom", fontweight="bold", fontsize=8)
                if row == 0:
                    ax.set_title(site, pad=5)
                if col == 0:
                    ax.set_ylabel(Y_LABELS[row], labelpad=3)
                if row == 2:
                    ax.set_xlabel(r"$R_H$ [mm/event]", labelpad=3)
        hd = [Line2D([], [], color=MODULE_COLOR[m], lw=1.2, label=PRETTY[m]) for m in ORDER]
        fig.legend(handles=hd, loc="upper center", bbox_to_anchor=(.55, .962), ncol=4,
                   handlelength=1.8, columnspacing=1.7)
        hd = [Line2D([], [], color="0.3", ls=ls, marker=mk, ms=3, label=o.capitalize())
              for o, ls, mk in (("portrait", "-", "o"), ("landscape", "--", "^"))]
        hd.append(Line2D([], [], color="0.6", lw=.8, ls=":", label=r"Illustrative reference: $R_H=15$"))
        fig.legend(handles=hd, loc="upper center", bbox_to_anchor=(.55, .924), ncol=3,
                   columnspacing=1.3)
        fig.suptitle(f"{'Ground' if market == 'ground' else 'Rooftop'} market · h = 60 mm · f = 0.1 · rectangular band",
                     fontsize=8, y=.987)
        scale_note = ("Neuchâtel: separate vertical scale in rows 1 and 3."
                      if market == "ground" else "Neuchâtel: separate vertical scale in all three rows.")
        fig.text(.55, .025, scale_note, ha="center", fontsize=7)
        # Exact final article dimensions; all curves/text remain vector in PDF.
        fig.savefig(path.with_suffix(".pdf"))
        fig.savefig(path.with_suffix(".png"), dpi=600)
        layout = dict(width_mm=WIDTH_MM, height_mm=HEIGHT_MM, market=market, panel_count=9,
                      data_curves_per_panel=8, reference_mm_event=REFERENCE_THRESHOLD,
                      separate_neuchatel_scale_rows_1based=[row + 1 for row in separate_rows],
                      almeria_doha_scales_shared=True,
                      y_limits={metric: {site: list(axes_limits[row][col]) for col, site in enumerate(energy.SITES)}
                                for row, metric in enumerate(METRICS)})
        # Check plotted arrays, not only the exported input frame.
        for row, metric in enumerate(METRICS):
            for col, site in enumerate(energy.SITES):
                curves = axes[row, col].lines[-8:]
                for line, (module, orientation) in zip(curves, [(m, o) for m in ORDER for o in ORIENTATIONS]):
                    g = q[(q.site == site) & (q.module == module) & (q.orientation == orientation)].sort_values("R_H_mm_event")
                    np.testing.assert_array_equal(line.get_xdata(), g.R_H_mm_event)
                    np.testing.assert_array_equal(line.get_ydata(), g[metric])
        plt.close(fig)
    return layout
