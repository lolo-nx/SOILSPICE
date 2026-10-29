"""Draw the resilience comparison using saved metrics, without rerunning SPICE.

    MPLBACKEND=Agg python chapters/fig_resilience.py
    MPLBACKEND=Agg python chapters/fig_resilience.py --horizontal-power-only
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import pymupdf

from soilspice import config, fig_simulation as F, resilience
from soilspice.style import ORDER


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def horizontal_power():
    source = config.SWEEP_DIR / "resilience_metrics.csv"
    output = config.CH2_DIR / "average_power_windows_horizontal"
    targets = (output.with_suffix(".pdf"), output.with_suffix(".png"))
    protected = {p: digest(p) for p in config.CH2_DIR.rglob("*") if p.is_file()
                 and p.suffix in (".csv", ".tex", ".pdf", ".png") and p not in targets}
    metrics = pd.read_csv(source)
    plotted = F.average_power_windows_horizontal(metrics, output)
    columns = [f"P_an_S_0to{w:g}" for w in resilience.PAN_WINDOWS]
    keys = ["module", "shape", "orientation"]
    original = metrics.set_index(keys).loc[plotted.set_index(keys).index, columns]
    np.testing.assert_array_equal(plotted[columns].to_numpy(), original.to_numpy())
    # Existing combined figure data retain the original vertical heatmap order.
    combined = pd.read_csv(config.SWEEP_DIR / "fig_resilience_data.csv").set_index(keys)
    np.testing.assert_array_equal(plotted[columns].to_numpy(), combined.loc[original.index, columns].to_numpy())
    with pymupdf.open(output.with_suffix(".pdf")) as doc:
        assert len(doc) == 1 and not doc[0].get_images()
        np.testing.assert_allclose([doc[0].rect.width, doc[0].rect.height],
                                  [183 / 25.4 * 72, 95 / 25.4 * 72], atol=1e-3)
    for path, sha in protected.items():
        assert digest(path) == sha, f"A pre-existing export changed: {path}"
    report = dict(source=str(source.relative_to(config.REPO_DIR)), source_sha256=digest(source),
                  width_mm=183, height_mm=95, png_dpi=600, vector_pdf=True,
                  configuration_columns=plotted[keys].to_dict("records"),
                  S_ref_rows_percent=list(resilience.PAN_WINDOWS),
                  all_128_values_match_original=True, transformation="transpose only",
                  M0_shown=False, physical_metrics_recalculated=False,
                  preexisting_exports_unchanged=len(protected))
    (config.SWEEP_DIR / "average_power_windows_horizontal.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Saved {output}.pdf/.png: 16 configuration columns × 8 optical-loss windows; original values preserved.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--horizontal-power-only", action="store_true", help="draw the transposed P_an,S table without M0")
    if parser.parse_args().horizontal_power_only:
        horizontal_power()
        return
    source = config.SWEEP_DIR / "resilience_metrics.csv"
    output = config.CH2_DIR / "fig_resilience"
    protected = {p: digest(p) for p in config.CH2_DIR.rglob("*") if p.is_file()
                 and p.suffix in (".csv", ".tex", ".pdf", ".png")
                 and p not in (output.with_suffix(".pdf"), output.with_suffix(".png"),
                               config.SWEEP_DIR / "fig_resilience_data.csv",
                               config.SWEEP_DIR / "fig_resilience_caption.tex")}
    metrics = pd.read_csv(source)
    metrics = metrics[metrics.module.isin(ORDER)]
    assert resilience.LOW_S_MAX == 2
    keys = [(name, shape, orientation) for shape, orientation in F.CONFIGS for name in ORDER]
    expected = metrics.set_index(["module", "shape", "orientation"]).loc[keys]
    # Read the original plot artists in memory to verify both source figures'
    # mappings and formatting; these plotters do not calculate physical metrics.
    captured = {}

    def capture(fig, path):
        if path.name == "original_m0_check":
            captured["m0"] = np.array([bar.get_height() for bar in fig.axes[0].patches])
        else:
            image = fig.axes[0].images[0]
            captured["power"] = np.asarray(image.get_array()).copy()
            captured["rgba"] = image.to_rgba(captured["power"])
            captured["annotations"] = [t.get_text() for t in fig.axes[0].texts][:128]
        F.plt.close(fig)

    with patch.object(F.plots, "save", side_effect=capture):
        F.initial_amplification(metrics, Path("original_m0_check"))
        F.average_power_windows(metrics, Path("original_power_check"))
    np.testing.assert_array_equal(captured["m0"], expected.M0.to_numpy())
    windows = resilience.PAN_WINDOWS
    cols = [f"P_an_S_0to{window:g}" for window in windows]
    np.testing.assert_array_equal(captured["power"], expected[cols].to_numpy())
    norm = F.Normalize(vmin=.5, vmax=1.)
    np.testing.assert_array_equal(captured["rgba"], F._average_power_colormap()(norm(expected[cols].to_numpy())))
    plotted = F.resilience_overview(metrics, output)
    np.testing.assert_array_equal(plotted.M0, captured["m0"])
    np.testing.assert_array_equal(plotted[cols], captured["power"])
    original_labels = [f"{v:.2f}" if v >= .995 else f"{v:.2f}"[1:] for v in captured["power"].flat]
    assert original_labels == captured["annotations"]
    # Reconcile original rounded table values independently as well.
    original_table = (config.CH2_DIR / "table_resilience_metrics.tex").read_text()
    for row in plotted.itertuples():
        assert f"{row.module} & {row.shape} & {row.orientation} & {row.M0:.2f} &" in original_table
    plotted.to_csv(config.SWEEP_DIR / "fig_resilience_data.csv", index=False)
    caption = (
        r"\caption{\textbf{Electrical resilience under bottom-edge soiling.} "
        r"(a) Initial loss amplification $M_0$, defined by the slope through the origin of relative "
        r"maximum-power loss versus optical loss for $0<S\leq2\%$, shown to one decimal place. "
        r"(b) Average normalised power $P_{\mathrm{an},S}$ over $[0,S_{\mathrm{ref}}]$, relative to the "
        r"uniform-soiling reference. Rows retain the four module and shape--orientation groups. "
        r"Both metrics are dimensionless; colours apply only to $P_{\mathrm{an},S}$. "
        r"All values are taken unchanged from the existing simulations.\label{fig:resilience}}"
    )
    (config.SWEEP_DIR / "fig_resilience_caption.tex").write_text(caption + "\n")
    with pymupdf.open(output.with_suffix(".pdf")) as doc:
        assert len(doc) == 1 and not doc[0].get_images()
        np.testing.assert_allclose([doc[0].rect.width, doc[0].rect.height],
                                  [183 / 25.4 * 72, 145 / 25.4 * 72], atol=1e-3)
    for p, sha in protected.items():
        assert digest(p) == sha, f"A pre-existing scientific export changed: {p}"
    report = dict(source=str(source.relative_to(config.REPO_DIR)), source_sha256=digest(source),
                  row_count=16, S_ref_percent=list(windows), M0_definition_max_S_percent=2,
                  all_16_M0_values_match_original_bar_artists_and_table=True,
                  all_128_heatmap_values_and_annotations_match_original=True,
                  heatmap_colours_and_normalisation_identical=True,
                  physical_metrics_recalculated=False, vector_pdf=True,
                  width_mm=183, height_mm=145, png_dpi=600,
                  preexisting_exports_unchanged=len(protected))
    (config.SWEEP_DIR / "fig_resilience_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(plotted[["module", "shape", "orientation", "M0", "M0_display"]].to_string(index=False), flush=True)
    print(f"Saved {output}.pdf/.png; all 16 M0 and 128 heatmap values verified.", flush=True)


if __name__ == "__main__":
    main()
