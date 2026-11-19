"""Separate normalized version of validation panel a; uses saved results only.

Run from the repository root: python chapters/normalized_validation_iv.py
Both EXP and SPICE use the clean measured flash-test Voc and Isc, preserving
soiling losses and model discrepancies rather than scaling each curve itself.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from soilspice import config, fig_validation
from soilspice.style import ORDER, ORIENTATIONS


def main():
    curves = pd.read_csv(config.CH1_DIR / "data" / "iv_grid_rectangle.csv")
    factors = []
    for name in ORDER:
        detail = pd.read_csv(config.VALIDATION_DIR / name / f"validation_detail_{name}.csv")
        for orientation in ORIENTATIONS:
            clean = detail[(detail["shape"] == "rectangle")
                           & (detail.orientation == orientation)
                           & np.isclose(detail.tau_I, 1.0)]
            reference = clean.set_index("parameter")["measured"]
            voc, isc = float(reference["Voc"]), float(reference["Isc"])
            if not np.isfinite([voc, isc]).all() or min(voc, isc) <= 0:
                raise ValueError(f"Invalid clean reference: {name}, {orientation}")
            mask = (curves.module == name) & (curves.orientation == orientation)
            curves.loc[mask, "V"] /= voc
            curves.loc[mask, "I"] /= isc
            factors.append(dict(module=name, orientation=orientation,
                                Voc_clean_V=voc, Isc_clean_A=isc))
    output = config.CH1_DIR / "iv_grid_rectangle_normalized"
    curves.to_csv(config.CH1_DIR / "data" / "iv_grid_rectangle_normalized.csv", index=False)
    pd.DataFrame(factors).to_csv(config.CH1_DIR / "data" / "iv_normalization_references.csv", index=False)
    fig_validation.iv_grid(curves, output, normalized=True)


if __name__ == "__main__":
    main()
