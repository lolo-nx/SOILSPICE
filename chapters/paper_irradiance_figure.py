"""Four-panel SI figures from saved simulations; no model changes or new runs."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from soilspice.style import STYLE, ORDER, module_label, save

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/stc_sensitivity'
GS = [200, 400, 600, 800, 1000]
TAUS = [0.71, 0.44, 0.0]
COLORS = {0.71: '#009E73', 0.44: '#0072B2', 0.0: '#D55E00'}
CONFIGS = [(shape, h, ori, t)
           for shape, h in [('rectangle', 80), ('triangle', 80), ('rectangle', 60)]
           for ori in ['portrait', 'landscape']
           for t in (TAUS if h == 80 else [0.44])]


def values(frame, module, config):
    shape, h, ori, tau = config
    q = frame[(frame.module == module) & (frame['shape'] == shape)
              & (frame.band_height_mm == h) & (frame.orientation == ori)
              & np.isclose(frame.tau_I, tau)].set_index('irradiance_W_m2')
    assert len(q) == 5 and set(q.index) == set(GS)
    return q.loc[GS, 'Delta_L_pp'].to_numpy()


def main():
    frame = pd.read_csv(OUT / 'simulation_results.csv')
    arrays = [np.array([values(frame, m, c) for c in CONFIGS]) for m in ORDER]
    assert sum(a.size for a in arrays) == 280
    assert all(np.allclose(a[:, -1], 0) for a in arrays)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 148/25.4), sharex=True, sharey=True)
        for k, (ax, module, arr) in enumerate(zip(axes.flat, ORDER, arrays)):
            ax.axhspan(-2, 2, color='.95', zorder=0)
            ax.axhline(0, color='.35', lw=.6)
            for config, ys in zip(CONFIGS, arr):
                shape, h, ori, tau = config
                marker = 's' if h == 60 else ('o' if shape == 'rectangle' else '^')
                ax.plot(GS, ys, color=COLORS[tau], ls='-' if ori == 'portrait' else '--',
                        marker=marker, ms=2.8, lw=.85, markerfacecolor='white' if h == 60 else COLORS[tau])
            ax.set_title(module_label(module), pad=6)
            ax.text(-.12, 1.04, 'abcd'[k], transform=ax.transAxes, fontweight='bold', fontsize=8)
            ax.set_xticks(GS)
            ax.set_yticks([-10, -8, -6, -4, -2, 0, 2])
            ax.set_ylim(-11.5, 3.2)
            ax.set_xlim(175, 1025)
            if module == 'SUNPOWER':
                i, j = np.unravel_index(np.argmin(arr), arr.shape)
                ax.annotate(f'{arr[i,j]:.2f} pp', (GS[j], arr[i,j]), xytext=(340, -9.7),
                            fontsize=7, arrowprops=dict(arrowstyle='-', lw=.6, color='.3'))
        handles = [Line2D([], [], color=COLORS[t], lw=1.3, label=rf'$\tau_I = {t:g}$') for t in TAUS]
        handles += [Line2D([], [], color='.25', ls=ls, label=ori) for ori, ls in [('Portrait', '-'), ('Landscape', '--')]]
        handles += [Line2D([], [], color='.25', ls='', marker=marker, ms=3,
                           markerfacecolor='white' if h == 60 else '.25', label=f'{shape}, {h} mm')
                    for shape, h, marker in [('Rectangle', 80, 'o'), ('Triangle', 80, '^'), ('Rectangle', 60, 's')]]
        fig.legend(handles=handles[:5], loc='upper center', bbox_to_anchor=(.54, .995), ncol=5)
        fig.legend(handles=handles[5:], loc='upper center', bbox_to_anchor=(.54, .96), ncol=3)
        fig.supxlabel('Irradiance (W m$^{-2}$)', y=.025, fontsize=7)
        fig.supylabel('Change in relative power loss (percentage points)', x=.015, fontsize=7)
        fig.subplots_adjust(left=.11, right=.98, bottom=.105, top=.85, hspace=.3, wspace=.13)
        save(fig, OUT / 'fig_deviation_from_STC_four_panels', dpi=300)

        fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 166/25.4), sharex=True)
        for k, (ax, module, arr) in enumerate(zip(axes.flat, ORDER, arrays)):
            im = ax.imshow(arr, cmap='RdBu_r', vmin=-11, vmax=11, aspect='auto', interpolation='nearest')
            ax.set_title(module_label(module), pad=7)
            ax.text(-.34, 1.035, 'abcd'[k], transform=ax.transAxes, fontweight='bold', fontsize=8)
            labels = [f'{"Rect." if s == "rectangle" else "Tri."} {h} | {"P" if o == "portrait" else "L"} | {t:g}'
                      for s, h, o, t in CONFIGS]
            ax.set_yticks(range(14), labels=labels, fontsize=5.6)
            ax.set_xticks(range(5), labels=GS)
            ax.tick_params(length=0, pad=3, labelbottom=True)
            ax.set_xticks(np.arange(-.5, 5, 1), minor=True)
            ax.set_yticks(np.arange(-.5, 14, 1), minor=True)
            ax.grid(which='minor', color='#cccccc', linewidth=.35)
            ax.tick_params(which='minor', length=0)
            for boundary in [2.5, 5.5, 8.5, 11.5]:
                ax.axhline(boundary, color='white', lw=.8)
            for i, j in zip(*np.where(np.abs(arr) > 2)):
                ax.text(j, i, f'{arr[i,j]:.2f}', ha='center', va='center', fontsize=5.8,
                        color='white' if np.isclose(arr[i, j], -10.4445999, atol=1e-5) else 'black')
            ax.set_xlabel('Irradiance (W m$^{-2}$)', fontsize=7)
        fig.subplots_adjust(left=.16, right=.98, bottom=.235, top=.91, hspace=.38, wspace=.59)
        cax = fig.add_axes([.31, .15, .5, .014])
        cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-10, -5, -2, 0, 2, 5, 10])
        cb.set_label('Change in relative power loss (percentage points)', fontsize=7)
        cb.outline.set_linewidth(.5)
        save(fig, OUT / 'fig_deviation_from_STC_heatmaps', dpi=300)
    (OUT / 'four_panel_figure_caption.md').write_text(
        'Effect of irradiance on relative power losses at a cell temperature of 25 °C. '
        'Panels show Yingli, SunPower, JA Solar and LONGi. Changes are calculated relative '
        'to the loss at 1000 W/m² for the same configuration. All 14 configurations per module '
        'are included: rectangular and triangular 80 mm bands at transmittances 0.71, 0.44 and 0, '
        'and a rectangular 60 mm band at transmittance 0.44, each in portrait and landscape. '
        'For the curves, the grey region indicates changes within ±2 percentage points and the '
        'horizontal line denotes zero change. For the heatmaps, row labels give shape, height '
        '(mm), orientation (P, portrait; L, landscape) and transmittance. Rect. denotes rectangular '
        'and Tri. denotes triangular. The difference is defined as '
        'ΔL(G) = L(G, 25 °C) − L(1000, 25 °C), in percentage points; negative values indicate '
        'lower relative losses than at the reference irradiance. The final two rows '
        'show the annual reference geometry, a 60 mm rectangle, tested here at τI = 0.44. '
        'Annual calculations use evolving band transmittances rather than a fixed value of 0.44. '
        'Heatmap values are printed '
        'where the absolute change exceeds 2 percentage points. The colour scale is shared '
        'and centred on zero. Of 224 conditions below 1000 W/m², 208 are within 2 percentage '
        'points of the reference. The largest change is −10.44 percentage points for SunPower '
        'in portrait with an opaque triangular 80 mm band at 200 W/m².\n', encoding='utf-8')
    print('Validated 280 points, 14 configurations per module, zero difference at 1000 W/m².')


if __name__ == '__main__':
    main()
