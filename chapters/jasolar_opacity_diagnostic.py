"""JA Solar 60 mm: SPICE response and energy-weighted exposure, saved PH15mm run.

Run with MPLBACKEND=Agg. Portrait states follow exactly from the rain-only
mass balance: tau_portrait = tau_landscape ** (L_portrait / L_landscape).
No climate resampling or new weather calculation is performed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from soilspice import config, energy
from soilspice.style import STYLE, save


def main():
    root = config.CH4_DIR / "total_scenarios" / "PH15mm"
    out = root / "data" / "opacity_diagnostic"
    out.mkdir(parents=True, exist_ok=True)
    series = pd.read_csv(root / "data" / "growth_series_reference.csv")
    orientations = ("landscape", "portrait")
    colors = dict(landscape="#2878a5", portrait="#d76a32")
    refs = [("JASOLAR", o, 60., .1, 0.) for o in orientations]
    curves = energy.lut(refs)
    ratio = energy.gains(refs)[1] / energy.gains(refs)[0]
    pd.DataFrame({"tau": energy.TAU, **dict(zip(orientations, 100 * curves))}).to_csv(
        out / "spice_curves.csv", index=False)
    rows, cross, distributions = [], [], []
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True)
        for j, site in enumerate(("Almería", "Doha")):
            s = series[series.site == site]
            y = s.Y.to_numpy()
            effective = y * (1 - s.L_u.to_numpy())
            for i, orientation in enumerate(orientations):
                tau = s.tau.to_numpy() ** (ratio if i else 1)
                idx = np.rint(tau * 200).astype(int)
                loss = curves[i, idx]
                # Chord through the clean and opaque endpoints, on the same LUT grid.
                linear = curves[i, 0] * (1 - energy.TAU[idx])
                hist, bins = np.histogram(tau, bins=np.linspace(0, 1, 21), weights=y / y.sum())
                distributions.extend(dict(site=site, orientation=orientation,
                                          tau_left=a, tau_right=b, energy_share=v)
                                     for a, b, v in zip(bins[:-1], bins[1:], hist))
                rows.append(dict(site=site, orientation=orientation,
                                 tau_mean=np.average(tau, weights=y),
                                 energy_share_tau_below_01=y[tau < .1].sum() / y.sum(),
                                 band_only_pct=100 * np.average(loss, weights=y),
                                 annual_increment_pct=100 * np.sum(effective * loss) / y.sum(),
                                 linear_increment_pct=100 * np.sum(effective * linear) / y.sum(),
                                 nonlinear_increment_pct=100 * np.sum(effective * (loss-linear)) / y.sum()))
                for k, response in enumerate(orientations):
                    cross.append(dict(site=site, distribution=orientation, response=response,
                                      band_only_pct=100*np.average(curves[k, idx], weights=y),
                                      annual_increment_pct=100*np.sum(effective*curves[k, idx])/y.sum()))
                axes[0, j].plot(energy.TAU, 100 * curves[i], color=colors[orientation],
                                label=orientation, lw=1.2)
                axes[1, j].stairs(100 * hist, bins, color=colors[orientation], lw=1.2)
            axes[0, j].set_title(site)
            axes[0, j].set_ylim(bottom=0)
            axes[1, j].set_xlabel(r"Band transmittance $\tau_I$")
            axes[1, j].set_xlim(0, 1)
        axes[0, 0].set_ylabel("SPICE power loss (%)")
        axes[1, 0].set_ylabel("Share of clean energy (%)\nper 0.05 transmittance bin")
        axes[0, 0].legend()
        fig.suptitle("JA Solar, rectangular band 60 mm, f = 0.1; washing at 15 mm/event", fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, .95))
        save(fig, root / "fig_jasolar_opacity_diagnostic")
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "summary.csv", index=False)
    pd.DataFrame(cross).to_csv(out / "crossed_distributions_responses.csv", index=False)
    pd.DataFrame(distributions).to_csv(out / "weighted_distributions.csv", index=False)
    print(summary.round(5).to_string(index=False))
    print(pd.DataFrame(cross).round(5).to_string(index=False))


if __name__ == "__main__":
    main()
