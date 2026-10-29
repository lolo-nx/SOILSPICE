"""Photos of the roof field test: (a) SPR1 without band, (b) SPR2 with the 80 mm band (tau_I = 0.71).

    soilspice.photos.figure(out_path)   (chapters/3_field_test.py)

The two photos are in data/field_test/photos/. The SPR2 photo is mirrored (left-right), so the band runs along the
left edge of the module and the label sits in the lower left corner. If a photo is changed, move the arrow with
ARROW_TAIL / ARROW_HEAD (axes fraction, 0 = left or bottom, 1 = right or top).
"""
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageOps

from . import config
from .style import LETTER_KW, STYLE, save

PHOTO_DIR = config.INPUT_DIR / "field_test" / "photos"
PHOTO_CLEAN = PHOTO_DIR / "SPR1_without_band.jpg"
PHOTO_BAND = PHOTO_DIR / "SPR2_with_band_mirror.jpg"

ASPECT = 296 / 338                  # height / width of each panel (same crop for both photos)
LABEL_XY = (0.04, 0.10)             # lower left corner of the band label (axes fraction)
ARROW_TAIL = (0.22, 0.20)           # top of the label
ARROW_HEAD = (0.27, 0.42)           # on the band
BOX = dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.9)


def load(path, aspect=ASPECT):
    """Photo with its EXIF rotation applied, centre-cropped to height / width = aspect."""
    img = ImageOps.exif_transpose(Image.open(path)).convert("RGB")
    w, h = img.size
    if h / w > aspect:
        nh = round(w * aspect)
        img = img.crop((0, (h - nh) // 2, w, (h - nh) // 2 + nh))
    else:
        nw = round(h / aspect)
        img = img.crop(((w - nw) // 2, 0, (w - nw) // 2 + nw, h))
    return np.asarray(img)


def figure(path, clean=PHOTO_CLEAN, band=PHOTO_BAND):
    with plt.rc_context(STYLE):
        width = 183 / 25.4
        gap = 0.02
        fig = plt.figure(figsize=(width, width * (1 - gap) / 2 * ASPECT))
        axes = [fig.add_axes((0, 0, (1 - gap) / 2, 1)), fig.add_axes(((1 + gap) / 2, 0, (1 - gap) / 2, 1))]
        for k, (ax, photo, title) in enumerate(zip(axes, (clean, band), ("SPR1, without band", "SPR2, with band"))):
            ax.imshow(load(photo), extent=(0, 1, 0, 1), aspect="auto", interpolation="lanczos")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(True)
                s.set_linewidth(0.5)
            ax.text(0.03, 0.965, "ab"[k], transform=ax.transAxes, va="top", ha="left", bbox=BOX, **LETTER_KW)
            ax.text(0.97, 0.965, title, transform=ax.transAxes, va="top", ha="right", fontsize=7,
                    fontweight="bold", bbox=BOX)
        ax = axes[1]
        ax.annotate("", xy=ARROW_HEAD, xytext=ARROW_TAIL, xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.0, mutation_scale=7))
        ax.text(*LABEL_XY, "band, $h$ = 80 mm\n$\\tau_I$ = 0.71", transform=ax.transAxes, va="bottom",
                ha="left", fontsize=6, linespacing=1.3, multialignment="center", bbox=BOX)
        save(fig, path)
