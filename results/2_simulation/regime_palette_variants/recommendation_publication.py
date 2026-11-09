"""Render the proposed publication palette and diagnostic viewing proofs.

Run from SOILSPICE with MPLBACKEND=Agg and a writable MPLCONFIGDIR.
The canonical manuscript figure is preserved.
"""
import json
from unittest.mock import patch

import numpy as np
import pandas as pd
import pymupdf
from matplotlib.colors import to_rgb

import build as comparison

OUT = comparison.OUT
PALETTE = {0: "#E6F1F5", 1: "#D1A6C0", 2: "#BFD2DE"}
# Machado et al. (2009), severity 100, from the colorspacious upstream data:
# https://github.com/njsmith/colorspacious/blob/master/colorspacious/cvd.py
MATRICES = {
    "Protanopie (simulation)": np.array([[.152286, 1.052583, -.204868],
                                        [.114503, .786281, .099216],
                                        [-.003882, -.048116, 1.051998]]),
    "Deuteranopie (simulation)": np.array([[.367322, .860646, -.227968],
                                          [.280085, .672501, .047413],
                                          [-.011820, .042940, .968881]]),
}


def linear(rgb):
    return np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)


def encoded(rgb):
    return np.where(rgb <= .0031308, 12.92 * rgb, 1.055 * np.maximum(rgb, 0) ** (1 / 2.4) - .055)


def main():
    protected = {p: comparison.sha(p) for p in comparison.config.CH2_DIR.rglob("*")
                 if p.is_file() and OUT not in p.parents}
    source = comparison.config.SWEEP_DIR / "resilience_states.csv"
    states = pd.read_csv(source)
    states = states[states.module.isin(comparison.ORDER)]
    baseline = json.loads((OUT / "I_III_proches" / "verification.json").read_text())
    assert comparison.sha(source) == baseline["source_sha256"]

    def save_checked(fig, *args, **kwargs):
        assert comparison.data_signature(fig) == baseline["data_signature"]
        comparison.SAVE(fig, *args, **kwargs)

    target = OUT / "PROPOSITION_PUBLICATION"
    with patch.dict(comparison.fig_regimes.REGIME_FILL, PALETTE), patch.object(comparison.plots, "save", save_checked):
        report = comparison.fig_regimes.draw(states, target)
    assert report["saved_states_plotted"] == report["scatter_points_plotted"] == 800
    assert report["atlas_panels"] == 16

    colours = np.array([to_rgb(c) for c in PALETTE.values()])
    luminance = linear(colours) @ np.array([.2126, .7152, .0722])
    report["CIELAB_lightness_D65"] = (116 * np.cbrt(luminance) - 16).tolist()
    report["contrast_black_text_on_regime_fill"] = ((luminance + .05) / .05).tolist()
    report["diagnostics_are_not_accessibility_certification"] = True
    report["cvd_source"] = "https://github.com/njsmith/colorspacious/blob/master/colorspacious/cvd.py"
    report["nature_guidance"] = [
        "https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/",
        "https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/",
    ]
    with pymupdf.open(target.with_suffix(".pdf")) as doc, pymupdf.open() as proof:
        p = doc[0]
        assert not p.get_images()
        p.get_pixmap(dpi=130).save(OUT / "PROPOSITION_PUBLICATION_preview.png")
        pix = p.get_pixmap(dpi=110, alpha=False)
        original = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3) / 255
        versions = [("Couleurs", pix), ("Niveaux de gris", p.get_pixmap(dpi=110, colorspace=pymupdf.csGRAY))]
        for name, matrix in MATRICES.items():
            data = np.ascontiguousarray(np.round(np.clip(encoded(linear(original) @ matrix.T), 0, 1) * 255).astype("uint8"))
            versions.append((name, pymupdf.Pixmap(pymupdf.csRGB, pix.width, pix.height, data.tobytes(), False)))
        page = proof.new_page(width=420 * comparison.MM, height=297 * comparison.MM)
        for i, (name, image) in enumerate(versions):
            x, y = i % 2 * 210, i // 2 * 148.5
            page.insert_text(((x + 4) * comparison.MM, (y + 7) * comparison.MM), name, fontsize=9)
            page.insert_image(pymupdf.Rect((x + 4) * comparison.MM, (y + 11) * comparison.MM,
                                           (x + 206) * comparison.MM, (y + 144) * comparison.MM), pixmap=image)
        proof.save(OUT / "PROPOSITION_PUBLICATION_verification_visuelle.pdf", garbage=4, deflate=True)
        page.get_pixmap(dpi=110).save(OUT / "PROPOSITION_PUBLICATION_verification_visuelle.png")
    for p, expected in protected.items():
        assert p.is_file() and comparison.sha(p) == expected, p
    report.update(source_sha256=comparison.sha(source), main_exports_preserved=True,
                  plotted_data_and_axis_limits_identical=True)
    (OUT / "PROPOSITION_PUBLICATION_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Publication proposal and diagnostic proofs generated; main exports preserved.", flush=True)


if __name__ == "__main__":
    main()
