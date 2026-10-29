"""Figure style of the paper, following the Nature Energy guidelines: Arial (Helvetica), 5-7 pt text, bold lowercase
panel letters (a, b, c) at 8 pt, thin lines (0.5 pt axes), no frame around legends, outward ticks, no top and right
axis lines, widths of 89 mm (one column, 3.5 in) or 183 mm (two columns, 7.2 in), colour-blind-safe colours
(Okabe-Ito), editable text in the PDF (Type 42 fonts)."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle

STYLE = {
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "font.size": 7, "axes.titlesize": 7, "axes.titleweight": "bold", "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6, "legend.title_fontsize": 6,
    "axes.linewidth": 0.5, "axes.edgecolor": "black", "axes.grid": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "lines.linewidth": 1.0, "lines.markersize": 3, "patch.linewidth": 0.5,
    "legend.frameon": False, "legend.handlelength": 1.6, "legend.borderaxespad": 0.3,
}
ORDER = ("YINGLI", "SUNPOWER", "JASOLAR", "LONGI")
ORIENTATIONS = ("portrait", "landscape")
MODULE_COLOR = {"YINGLI": "#009E73", "SUNPOWER": "#D55E00", "JASOLAR": "#E69F00", "LONGI": "#0072B2"}   # Okabe-Ito
PRETTY = {"YINGLI": "Yingli", "SUNPOWER": "SunPower", "JASOLAR": "JA Solar", "LONGI": "LONGi"}
TECH = {"YINGLI": "n-Si", "SUNPOWER": "IBC", "JASOLAR": "TOPCon", "LONGI": "HJT BC"}
EDGE = dict(edgecolor="black", lw=0.4)
LETTER = 8                                   # panel letters: 8 pt bold lowercase
LETTER_KW = dict(fontsize=LETTER, fontweight="bold")
# the three sites: line style on the maps, fill on the bars (as in the first version of the bars)
SITES = [("Neuchâtel", "Neuchâtel (ERA5)", "-", "#9db4c0"),
         ("Almería", "Almería (ERA5)", "--", "#e3a857"),
         ("Doha", "Doha (ERA5)", ":", "#b5452f")]
# binned cycle loss, sampled from YlOrBr (the P_mp loss maps of the paper), white for < 0.5 %
BOUNDS = [0, 0.5, 1, 2, 4, 8, 16, 100]
_ylorbr = plt.get_cmap("YlOrBr")
COLORS = ["#ffffff"] + [tuple(_ylorbr(v)[:3]) for v in np.linspace(0.12, 0.95, len(BOUNDS) - 2)]
LOSS_CMAP, LOSS_NORM = ListedColormap(COLORS), BoundaryNorm(BOUNDS, len(COLORS))


def module_label(name, with_tech=True):
    return f"{PRETTY[name]} ({TECH[name]})" if with_tech else PRETTY[name]


def letter(ax, k, x=0.03, y=0.95, va="top", **kw):
    ax.text(x, y, "abcdefghijklmnop"[k], transform=ax.transAxes, va=va, **{**LETTER_KW, **kw})


def bold(text):
    return r"$\bf{" + text.replace(" ", r"\ ") + r"}$"


def edges(values):
    v = np.asarray(values, dtype=float)
    mid = (v[1:] + v[:-1]) / 2
    return np.r_[v[0] - (mid[0] - v[0]), mid, v[-1] + (v[-1] - mid[-1])]


def module_handles():
    return [Rectangle((0, 0), 1, 1, facecolor=MODULE_COLOR[n], label=module_label(n), **EDGE) for n in ORDER]


def top_legend(fig, handles, ncol=None, y=1.0):
    fig.legend(handles=handles, loc="upper center", ncol=ncol or len(handles), bbox_to_anchor=(0.5, y),
               frameon=False, handlelength=1.4, handletextpad=0.4, columnspacing=1.2)


def group_axis(ax, x, y=-0.14):
    """x axis of grouped bars: orientation under each bar group, module in bold under each pair."""
    ax.set_xticks(x)
    ax.set_xticklabels([o.capitalize() for _ in ORDER for o in ORIENTATIONS])
    ax.tick_params(axis="x", length=0)
    for k, m in enumerate(ORDER):
        ax.text(x[2 * k:2 * k + 2].mean(), y, PRETTY[m], transform=ax.get_xaxis_transform(), ha="center",
                va="top", fontweight="bold")
        if k:
            ax.axvline((x[2 * k - 1] + x[2 * k]) / 2, color="0.6", lw=0.6, ls="--")


def loss_colorbar(fig, axes, label):
    sm = plt.cm.ScalarMappable(cmap=LOSS_CMAP, norm=LOSS_NORM)
    cb = fig.colorbar(sm, ax=axes, orientation="horizontal", location="top", fraction=0.04, pad=0.07, aspect=60,
                      spacing="uniform", ticks=BOUNDS[:-1])
    cb.ax.set_xticklabels([f"{b:g}" for b in BOUNDS[:-1]])
    cb.set_label(label, labelpad=2)
    cb.outline.set_linewidth(0.5)
    return cb


def save(fig, path, dpi=600):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path.with_suffix('.png')}")


def letter_above(ax, k, x=-0.12, y=1.02):
    """Panel letter above the upper left corner of the axes (Nature style)."""
    ax.text(x, y, "abcdefghijklmnop"[k], transform=ax.transAxes, va="bottom", ha="left", **LETTER_KW)


def letter_inside(ax, k, x=0.03, y=0.97):
    """Panel letter in the upper left corner of the axes."""
    ax.text(x, y, "abcdefghijklmnop"[k], transform=ax.transAxes, va="top", **LETTER_KW)

# site colours of the curves and bars of chapter 4
SITE_COLOR = {"Neuchâtel": "#3A6EA5", "Almería": "#E0A458", "Doha": "#A44A3F"}   # wet -> dry: blue, ochre, terracotta
