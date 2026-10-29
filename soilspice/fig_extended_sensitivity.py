"""Extended Data sensitivity figures from saved, corrected chapter-4 results only."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from .energy import SITES, TOTAL_THRESHOLDS, REFERENCE_TOTAL
from .style import SITE_COLOR, STYLE

WIDTH_MM, HEIGHT_MM = 183., 125.
LOW_COLOR, HIGH_COLOR = "#0072B2", "#E69F00"
PARAMETERS = (
    ("f", "Trapped fraction $f$\n0.03 / 0.3"),
    ("rate", "Soiling rate $k$\n×0.5 / ×2"),
    ("alpha", "Optical attenuation $\\alpha$\n×0.75 / ×1.25"),
    ("light", "Transfer threshold $R_L$\n0.5 / 2 mm/event"),
    ("cost", "Cleaning cost $C$\n×0.5 / ×2"),
)
METRICS = ("energy_none_band_extra_loss_pct", "margin10_npv25_eur_kwp")
PLOT_STYLE = {**STYLE, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 7}


def finish(fig, axes, path):
    """Keep exact article dimensions and check text is not clipped by the page."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text = list(fig.texts)
    for ax in axes.flat:
        text.extend([ax.title, ax.xaxis.label, ax.yaxis.label, *ax.texts])
        # Locators may create ticks beyond the view limits; Matplotlib does not
        # draw those labels, so only inspect the ticks actually on the page.
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks():
                if lo - 1e-12 <= tick.get_loc() <= hi + 1e-12:
                    text.extend([tick.label1, tick.label2])
    for legend in fig.legends:
        text.extend(legend.get_texts())
    for item in text:
        if item.get_visible() and item.get_text():
            box = item.get_window_extent(renderer)
            if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                raise ValueError(f"Text outside final figure: {item.get_text()!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=600)
    result = dict(width_mm=WIDTH_MM, height_mm=HEIGHT_MM, panels=6, png_dpi=600,
                  y_limits={"abcdef"[i]: list(ax.get_ylim()) for i, ax in enumerate(axes.flat)})
    plt.close(fig)
    return result


def parameter_sensitivity(cases, path):
    """Paired horizontal bars extend from the reference; endpoints are saved medians."""
    with plt.rc_context(PLOT_STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), sharey=True)
        fig.subplots_adjust(left=.205, right=.985, top=.86, bottom=.13, hspace=.57, wspace=.22)
        for row, metric in enumerate(("band_loss", "band_cost")):
            for col, site in enumerate(SITES):
                ax = axes[row, col]
                values = cases[cases.site == site].set_index("case")[metric]
                reference = values["reference"]
                for i, (parameter, _) in enumerate(PARAMETERS):
                    for level, offset, color in (("low", -.16, LOW_COLOR), ("high", .16, HIGH_COLOR)):
                        value = values[f"{parameter} {level}"]
                        bars = ax.barh(i + offset, value - reference, left=reference,
                                       height=.28, facecolor=color, edgecolor="black", lw=.3, zorder=2)
                        np.testing.assert_allclose(bars[0].get_x() + bars[0].get_width(), value, atol=1e-12)
                ax.axvline(reference, color="black", lw=.9, zorder=4)
                ax.set_xlim(0, max(values) * 1.10)
                ax.set_ylim(len(PARAMETERS) - .48, -.52)
                ax.set_yticks(range(len(PARAMETERS)))
                ax.set_yticklabels([label for _, label in PARAMETERS], fontsize=7)
                ax.tick_params(axis="y", length=0, pad=5)
                ax.tick_params(axis="x", pad=2)
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
                ax.text(.01, 1.035, "abcdef"[row * 3 + col], transform=ax.transAxes,
                        fontweight="bold", fontsize=8, va="bottom")
                if row == 0:
                    ax.set_title(site, pad=7)
        fig.legend(handles=[Patch(facecolor=LOW_COLOR, edgecolor="black", lw=.3, label="Lower value"),
                            Patch(facecolor=HIGH_COLOR, edgecolor="black", lw=.3, label="Higher value"),
                            Line2D([], [], color="black", lw=.9, label="Reference")],
                   loc="upper center", bbox_to_anchor=(.59, .985), ncol=3, columnspacing=2)
        fig.text(.59, .505, "Additional annual band energy loss [% of clean yield]", ha="center", fontsize=7)
        fig.text(.59, .058, "Additional annual soiling cost, trigger rule [% of clean revenue]",
                 ha="center", fontsize=7)
        return finish(fig, axes, path)


def washing_threshold_sensitivity(aggregates, path):
    """Only the five computed thresholds: median and configuration min--max envelope."""
    with plt.rc_context(PLOT_STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), sharex=True)
        fig.subplots_adjust(left=.12, right=.985, top=.86, bottom=.14, hspace=.35, wspace=.32)
        for row, metric in enumerate(METRICS):
            # A shared economic scale is readable. The small Neuchâtel energy loss
            # needs its own scale; Almería and Doha retain the same energy scale.
            other = aggregates[aggregates.site != "Neuchâtel"] if row == 0 else aggregates
            shared_max = other[f"{metric}_max"].max()
            shared_min = min(0., other[f"{metric}_min"].min())
            for col, site in enumerate(SITES):
                ax = axes[row, col]
                q = aggregates[aggregates.site == site].sort_values("R_H_mm_event")
                x = q.R_H_mm_event.to_numpy()
                np.testing.assert_array_equal(x, TOTAL_THRESHOLDS)
                lo, med, hi = [q[f"{metric}_{stat}"].to_numpy() for stat in ("min", "median", "max")]
                ax.axvline(REFERENCE_TOTAL, color="0.5", ls=":", lw=.8, zorder=1)
                if row == 1:
                    ax.axhline(0, color="0.3", lw=.6, zorder=1)
                ax.fill_between(x, lo, hi, facecolor=SITE_COLOR[site], alpha=.20, lw=0, zorder=2)
                line, = ax.plot(x, med, color=SITE_COLOR[site], marker="o", ms=3.2, lw=1.15, zorder=3)
                np.testing.assert_array_equal(line.get_xdata(), x)
                np.testing.assert_array_equal(line.get_ydata(), med)
                upper = hi.max() if row == 0 and col == 0 else shared_max
                lower = shared_min
                span = upper - lower
                ax.set_ylim(0 if row == 0 else lower - .07 * span, upper + .10 * span)
                ax.set_xlim(4, 26)
                ax.set_xticks(TOTAL_THRESHOLDS)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
                ax.tick_params(pad=2)
                ax.text(.015, 1.035, "abcdef"[row * 3 + col], transform=ax.transAxes,
                        fontweight="bold", fontsize=8, va="bottom")
                if row == 0:
                    ax.set_title(site, color=SITE_COLOR[site], pad=7)
                else:
                    ax.set_xlabel(r"$R_H$ [mm/event]", labelpad=4)
                if col == 0:
                    ax.set_ylabel(("Additional annual band energy loss\n[% of clean yield]" if row == 0 else
                                   "Net present value of +10 mm margin\n[€ kWp$^{-1}$, 25 years]"), labelpad=5)
        fig.legend(handles=[Line2D([], [], color="0.25", marker="o", ms=3, lw=1.15, label="Median"),
                            Patch(facecolor="0.5", alpha=.20, lw=0, label="Min–max across configurations"),
                            Line2D([], [], color="0.5", ls=":", lw=.8, label=r"Reference: $R_H=15$ mm/event")],
                   loc="upper center", bbox_to_anchor=(.55, .985), ncol=3, columnspacing=1.4)
        return finish(fig, axes, path)
