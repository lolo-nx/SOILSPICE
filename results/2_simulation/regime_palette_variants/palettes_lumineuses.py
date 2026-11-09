"""Compare three fresher colour proposals against the previous publication proposal."""
import json
from unittest.mock import patch

import numpy as np
import pandas as pd
import pymupdf
from matplotlib.colors import to_rgb

import build as comparison
from recommendation_publication import MATRICES, linear, encoded

OUT = comparison.OUT / "Plus_lumineux"
PALETTES = [
    ("A", "Lagon et framboise", ("#E1F5F4", "#EF9FBD", "#A3DADB")),
    ("B", "Ciel et rose", ("#E1F1FF", "#E995B5", "#AFD1F1")),
    ("C", "Lavande et menthe", ("#F0ECFA", "#B4E2C0", "#C9B9EB")),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    protected = {p: comparison.sha(p) for p in comparison.config.CH2_DIR.rglob("*")
                 if p.is_file() and comparison.OUT not in p.parents}
    source = comparison.config.SWEEP_DIR / "resilience_states.csv"
    states = pd.read_csv(source)
    states = states[states.module.isin(comparison.ORDER)]
    baseline = json.loads((comparison.OUT / "I_III_proches/verification.json").read_text())
    assert comparison.sha(source) == baseline["source_sha256"]
    report = {"source_sha256": comparison.sha(source), "regime_order": ["I", "II", "III"], "variants": []}

    def save_checked(fig, *args, **kwargs):
        assert comparison.data_signature(fig) == baseline["data_signature"]
        comparison.SAVE(fig, *args, **kwargs)

    original_cell = comparison.fig_regimes._atlas_cell

    def quieter_grid(ax, *args, **kwargs):
        result = original_cell(ax, *args, **kwargs)
        # Cover antialiasing seams with an almost invisible stroke in the cell's
        # own colour; cell coordinates and categorical values are unchanged.
        ax.collections[0].set_edgecolor("face")
        ax.collections[0].set_linewidth(.05)
        return result

    for letter, title, colours in PALETTES:
        palette = dict(enumerate(colours))
        target = OUT / f"{letter}_{title.lower().replace(' ', '_')}"
        with patch.dict(comparison.fig_regimes.REGIME_FILL, palette), \
                patch.object(comparison.fig_regimes, "_atlas_cell", quieter_grid), \
                patch.object(comparison.plots, "save", save_checked):
            checked = comparison.fig_regimes.draw(states, target)
        assert checked["saved_states_plotted"] == checked["scatter_points_plotted"] == 800
        assert checked["atlas_panels"] == 16
        y = linear(np.array([to_rgb(c) for c in colours])) @ np.array([.2126, .7152, .0722])
        report["variants"].append(dict(letter=letter, title=title, palette=palette, pdf=target.name + ".pdf",
                                       lightness=(116 * np.cbrt(y) - 16).tolist(),
                                       black_text_contrast=((y + .05) / .05).tolist()))
        print(f"Rendered {letter}: {title}", flush=True)

    reference = comparison.OUT / "PROPOSITION_PUBLICATION.pdf"
    reference_report = json.loads((comparison.OUT / "PROPOSITION_PUBLICATION_verification.json").read_text())
    entries = [dict(letter="R", title="Proposition precedente", palette=reference_report["palette"], pdf=str(reference))]
    entries += [{**v, "pdf": str(OUT / v["pdf"])} for v in report["variants"]]
    mm = comparison.MM
    with pymupdf.open() as board, pymupdf.open() as pages:
        p = board.new_page(width=400 * mm, height=400 * mm)
        for i, entry in enumerate(entries):
            x, y = i % 2 * 200, i // 2 * 200
            p.insert_text(((x + 4) * mm, (y + 7) * mm), f"{entry['letter']} | {entry['title']}",
                          fontsize=10, fontname="hebo")
            for j, colour in enumerate(entry["palette"].values()):
                rect = pymupdf.Rect((x + 4 + j * 24) * mm, (y + 10) * mm,
                                    (x + 25 + j * 24) * mm, (y + 14) * mm)
                p.draw_rect(rect, fill=to_rgb(colour), color=(.4, .4, .4), width=.4)
                p.insert_text((rect.x0 + 2 * mm, rect.y1 - mm), ["I", "II", "III"][j], fontsize=7)
            with pymupdf.open(entry["pdf"]) as doc:
                p.show_pdf_page(pymupdf.Rect((x + 3) * mm, (y + 18) * mm, (x + 197) * mm, (y + 196) * mm), doc, 0)
                pages.insert_pdf(doc)
        board.save(OUT / "COMPARAISON.pdf", garbage=4, deflate=True)
        p.get_pixmap(dpi=125).save(OUT / "COMPARAISON.png")
        pages.set_toc([[1, f"{e['letter']} - {e['title']}", i + 1] for i, e in enumerate(entries)])
        pages.save(OUT / "VERSIONS_COMPLETES.pdf", garbage=4, deflate=True)

    with pymupdf.open() as proof:
        for entry in report["variants"]:
            with pymupdf.open(OUT / entry["pdf"]) as doc:
                page = doc[0]
                assert not page.get_images()
                page.get_pixmap(dpi=125).save(OUT / (entry["letter"] + "_preview.png"))
                pix = page.get_pixmap(dpi=105, alpha=False)
                original = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3) / 255
                views = [("Couleurs", pix), ("Gris", page.get_pixmap(dpi=105, colorspace=pymupdf.csGRAY))]
                for name, matrix in MATRICES.items():
                    data = np.ascontiguousarray(np.round(np.clip(encoded(linear(original) @ matrix.T), 0, 1) * 255).astype("uint8"))
                    views.append((name, pymupdf.Pixmap(pymupdf.csRGB, pix.width, pix.height, data.tobytes(), False)))
                p = proof.new_page(width=420 * mm, height=297 * mm)
                for j, (name, image) in enumerate(views):
                    x, y = j % 2 * 210, j // 2 * 148.5
                    p.insert_text(((x + 4) * mm, (y + 7) * mm), f"{entry['letter']} | {entry['title']} | {name}", fontsize=9)
                    p.insert_image(pymupdf.Rect((x + 4) * mm, (y + 11) * mm, (x + 206) * mm, (y + 144) * mm), pixmap=image)
        proof.save(OUT / "VERIFICATION_VISUELLE.pdf", garbage=4, deflate=True)

    for p, expected in protected.items():
        assert p.is_file() and comparison.sha(p) == expected, p
    report.update(main_exports_preserved=True, plotted_data_and_axes_identical=True,
                  atlas_white_grid_removed=True, diagnostic_cvd_model="Machado 2009, severity 100")
    (OUT / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    (OUT / "README.md").write_text(
        "# Propositions plus lumineuses\n\n"
        "Trois variantes conservent I et III dans une meme famille et distinguent II.\n"
        "Les cartes utilisent des aplats sans quadrillage blanc ; les contours et les donnees sont conserves.\n\n"
        "[Comparaison avec la proposition precedente](COMPARAISON.pdf) · "
        "[Versions completes](VERSIONS_COMPLETES.pdf) · [Epreuves de diagnostic](VERIFICATION_VISUELLE.pdf)\n\n"
        "Les figures principales ne sont pas modifiees. Chaque variante conserve les 800 etats et 16 cartes.\n\n"
        "Reproduction depuis SOILSPICE :\n\n```sh\n"
        "MPLCONFIGDIR=/tmp/soilspice-matplotlib MPLBACKEND=Agg python3 results/2_simulation/regime_palette_variants/palettes_lumineuses.py\n```\n")
    print("Three brighter proposals complete; main figure preserved.", flush=True)


if __name__ == "__main__":
    main()
