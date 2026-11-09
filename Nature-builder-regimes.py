"""Draw the regime map and atlas as an integrated 183 × 170 mm vector figure.

Run: python Nature-builder-regimes.py
The full response scatter and its colour and symbol legend sit above the atlas.
The atlas uses module columns and vertical configuration rows, with its contour key below.
No metric is recalculated.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import numpy as np
import pymupdf

from soilspice import config, fig_regimes
from soilspice.style import ORDER
from soilspice.analysis import regime


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=config.CH2_DIR / "fig_regimes_panel")
    args = parser.parse_args()
    source = config.SWEEP_DIR / "resilience_states.csv"
    if not source.is_file():
        parser.error(f"Source data not found: {source}")
    outputs = {args.output.with_suffix(ext) for ext in (".pdf", ".png")}
    protected = {p: digest(p) for p in config.CH2_DIR.rglob("*") if p.is_file()
                 and p.suffix in (".csv", ".tex", ".pdf", ".png") and p not in outputs}
    states = pd.read_csv(source)
    states = states[states["module"].isin(ORDER)]
    assert not states.duplicated(["module", "shape", "orientation", "A_cov", "tau_I"]).any()
    np.testing.assert_array_equal(states.regime,
                                  [regime(v, i) for v, i in zip(states["Vmp_loss_%"], states["Imp_loss_%"])])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = fig_regimes.draw(states, args.output)
    assert report["scatter_points_plotted"] == report["saved_states_plotted"] == len(states)
    with pymupdf.open(args.output.with_suffix(".pdf")) as doc:
        assert len(doc) == 1 and not doc[0].get_images()
        np.testing.assert_allclose([doc[0].rect.width, doc[0].rect.height],
                                  [183 / 25.4 * 72, 170 / 25.4 * 72], atol=1e-3)
    for p, expected in protected.items():
        assert digest(p) == expected, f"A pre-existing data file or standalone figure changed: {p}"
    report.update(source=str(source.relative_to(config.REPO_DIR)), source_sha256=digest(source),
                  preexisting_exports_unchanged=len(protected),
                  standalone_regime_map_and_atlas_preserved=True,
                  pdf_sha256=digest(args.output.with_suffix(".pdf")),
                  png_sha256=digest(args.output.with_suffix(".png")))
    (config.SWEEP_DIR / "fig_regimes_panel_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    caption = (
        r"\caption{\textbf{Electrical regimes under bottom-edge soiling.} "
        r"\textbf{a}, Simulated changes in maximum-power voltage and current relative to the clean module; "
        rf"all {len(states)} saved states are shown. Background colours distinguish current-limited, "
        r"bypass-clamped and mixed regimes. The legend identifies all four band configurations: "
        r"rectangles and triangles denote band shape, with open symbols for portrait and filled symbols "
        r"for landscape. \textbf{b}, Corresponding regime maps versus band transmittance $\tau_I$ "
        r"and active-cell coverage $A_{\mathrm{cov}}$. Coloured frames around the single-line module "
        r"names and technologies key the module colours; "
        r"vertical row labels use the same four configuration symbols as a. Black contours show "
        r"relative maximum-power change in percent, keyed below b; negative values indicate power loss. "
        r"All classifications and power values are unchanged from the saved simulations.\label{fig:regimes}}"
    )
    (config.SWEEP_DIR / "fig_regimes_panel_caption.tex").write_text(caption + "\n")
    print(f"Verified {len(states)} saved states, 16 atlas panels, unchanged classifications and source figures.")


if __name__ == "__main__":
    main()
