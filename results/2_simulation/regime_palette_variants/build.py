"""Compare supplied pale regime palettes without changing the main exports.

From SOILSPICE: MPLBACKEND=Agg python3 results/2_simulation/regime_palette_variants/build.py
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pymupdf

OUT = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from soilspice import config, fig_regimes, plots
from soilspice.analysis import regime
from soilspice.style import ORDER

COLOURS = {"bleu": "#AACCEE", "cyan": "#CCEEFF", "vert": "#BBDDBB",
           "jaune": "#EEEEBB", "rose": "#FFBBCC", "violet": "#EEBBDD", "gris": "#DDDDDD"}
TRIOS = [
    ("A", "cyan_jaune_violet", ("cyan", "jaune", "violet")),
    ("B", "cyan_vert_rose", ("cyan", "vert", "rose")),
    ("C", "gris_cyan_rose", ("gris", "cyan", "rose")),
    ("D", "cyan_jaune_vert", ("cyan", "jaune", "vert")),
    ("E", "bleu_jaune_rose", ("bleu", "jaune", "rose")),
    ("F", "gris_vert_violet", ("gris", "vert", "violet")),
]
RELATED_TRIOS = [
    ("A", "cyan_jaune_bleu", ("cyan", "jaune", "bleu")),
    ("B", "cyan_rose_bleu", ("cyan", "rose", "bleu")),
    ("C", "cyan_vert_bleu", ("cyan", "vert", "bleu")),
    ("D", "rose_cyan_violet", ("rose", "cyan", "violet")),
    ("E", "rose_vert_violet", ("rose", "vert", "violet")),
    ("F", "rose_jaune_violet", ("rose", "jaune", "violet")),
]
MM = 72 / 25.4
SAVE = plots.save


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_signature(fig):
    """Ignore colours; retain axis limits, response points and categorical grids."""
    digest = hashlib.sha256()
    for ax in fig.axes:
        digest.update(np.asarray([*ax.get_xlim(), *ax.get_ylim()], dtype=float).tobytes())
        for line in ax.lines:
            digest.update(np.asarray(line.get_xydata(), dtype=float).tobytes())
        for collection in ax.collections:
            values = collection.get_array()
            if values is not None:
                digest.update(np.asarray(np.ma.filled(values, np.nan)).tobytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--related", action="store_true", help="Similar hues for I and III; distinct hue for II")
    args = parser.parse_args()
    out = OUT / "I_III_proches" if args.related else OUT
    trios = RELATED_TRIOS if args.related else TRIOS
    out.mkdir(parents=True, exist_ok=True)
    protected = {p: sha(p) for p in config.CH2_DIR.rglob("*")
                 if p.is_file() and OUT not in p.parents}
    source = config.SWEEP_DIR / "resilience_states.csv"
    states = pd.read_csv(source)
    states = states[states.module.isin(ORDER)]
    np.testing.assert_array_equal(states.regime,
                                  [regime(v, i) for v, i in zip(states["Vmp_loss_%"], states["Imp_loss_%"])])
    signatures = []

    def save_checked(fig, *args, **kwargs):
        signatures.append(data_signature(fig))
        if len(signatures) > 1:
            assert signatures[-1] == signatures[0], "Variant changed plotted data or axis limits"
        SAVE(fig, *args, **kwargs)

    report = {"source": str(source.relative_to(ROOT)), "source_sha256": sha(source),
              "regime_order": ["I: current limited", "II: bypass clamped", "III: mixed"],
              "supplied_colours": COLOURS, "I_III_related_hues": args.related, "variants": []}
    for letter, slug, names in trios:
        palette = {i: COLOURS[name] for i, name in enumerate(names)}
        target = out / f"{letter}_{slug}"
        with patch.dict(fig_regimes.REGIME_FILL, palette), patch.object(plots, "save", save_checked):
            checked = fig_regimes.draw(states, target)
        assert checked["scatter_points_plotted"] == checked["saved_states_plotted"] == len(states)
        assert checked["atlas_panels"] == 16
        with pymupdf.open(target.with_suffix(".pdf")) as doc:
            assert not doc[0].get_images()
        report["variants"].append({"letter": letter, "names_I_II_III": names,
                                   "palette": palette, "pdf": target.name + ".pdf",
                                   "verified_states": len(states), "atlas_panels": 16})
        print(f"Variant {letter}: {' / '.join(names)}", flush=True)

    with pymupdf.open() as board, pymupdf.open() as all_pages:
        page = board.new_page(width=420 * MM, height=297 * MM)
        for index, variant in enumerate(report["variants"]):
            x, y = index % 3 * 140, index // 3 * 148.5
            title = f"{variant['letter']} | " + " / ".join(variant["names_I_II_III"])
            if index == 0 and not args.related:
                title += " (actuel)"
            page.insert_text(((x + 4) * MM, (y + 7) * MM), title, fontsize=9, fontname="hebo")
            for number, colour in variant["palette"].items():
                rgb = tuple(int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
                rect = pymupdf.Rect((x + 4 + number * 23) * MM, (y + 10) * MM,
                                   (x + 24 + number * 23) * MM, (y + 14) * MM)
                page.draw_rect(rect, fill=rgb, color=(.4, .4, .4), width=.4)
                page.insert_text((rect.x0 + 2 * MM, rect.y1 - 1 * MM),
                                 ["I", "II", "III"][number], fontsize=7)
            with pymupdf.open(out / variant["pdf"]) as figure:
                page.show_pdf_page(pymupdf.Rect((x + 3) * MM, (y + 18) * MM,
                                               (x + 137) * MM, (y + 144.5) * MM), figure, 0)
                all_pages.insert_pdf(figure)
        all_pages.set_toc([[1, f"{v['letter']} - {' / '.join(v['names_I_II_III'])}", i + 1]
                           for i, v in enumerate(report["variants"])])
        board.save(out / "COMPARAISON_TRIOS.pdf", garbage=4, deflate=True)
        page.get_pixmap(dpi=140).save(out / "COMPARAISON_TRIOS.png")
        all_pages.save(out / "TOUS_LES_TRIOS.pdf", garbage=4, deflate=True)

    for path, expected in protected.items():
        assert path.is_file() and sha(path) == expected, path
    report.update(identical_plotted_data_and_axes=True, data_signature=signatures[0],
                  preexisting_files_unchanged=len(protected), original_exports_preserved=True)
    (out / "verification.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    rows = ["# Variantes des couleurs de régime", "",
            "Les couleurs sont exactement celles de la palette fournie. L’ordre est I / II / III.",
            ("I et III utilisent des teintes proches ; II utilise une autre couleur. Les figures principales sont conservées."
             if args.related else "A reprend la proposition actuelle ; les figures principales ne sont pas modifiées."), "",
            "[Planche comparative](COMPARAISON_TRIOS.pdf) · [Six pages à taille complète](TOUS_LES_TRIOS.pdf)", "",
            "| Variante | I : current limited | II : bypass clamped | III : mixed | Figure |",
            "|---|---|---|---|---|"]
    for v in report["variants"]:
        cells = [f"{name} `{v['palette'][i]}`" for i, name in enumerate(v["names_I_II_III"])]
        rows.append(f"| {v['letter']} | {' | '.join(cells)} | [PDF]({v['pdf']}) |")
    rows += ["", "Chaque variante conserve les 800 états, les 16 cartes, les limites d’axes et la mise en page.",
             "Les PDF sont vectoriels et les PNG individuels sont exportés à 600 dpi.", "",
             "Reproduction depuis SOILSPICE :", "",
             "```sh", "MPLCONFIGDIR=/tmp/soilspice-matplotlib MPLBACKEND=Agg python3 results/2_simulation/regime_palette_variants/build.py" + (" --related" if args.related else ""), "```", ""]
    (out / "README.md").write_text("\n".join(rows))
    print("Six variants and comparison sheet complete; original exports preserved.", flush=True)


if __name__ == "__main__":
    main()
