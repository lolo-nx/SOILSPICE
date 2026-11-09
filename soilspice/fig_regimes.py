"""Integrated regime map and atlas, laid out at their final 183 × 170 mm size."""
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from . import plots
from .analysis import REGIME_DI_PCT, REGIME_DV_PCT
from .fig_simulation import CONFIGS, REGIME_EDGE, REGIME_FILL, SHAPE_MARKER, sel
from .style import MODULE_COLOR, ORDER, PRETTY, STYLE, TECH

WIDTH_MM, HEIGHT_MM = 183.0, 170.0
CONTOURS = (-80, -60, -45, -35, -15, -5)


def _rect(x, y, width, height):
    """Figure rectangle specified in mm, measured from the top-left corner."""
    return [x / WIDTH_MM, 1 - (y + height) / HEIGHT_MM,
            width / WIDTH_MM, height / HEIGHT_MM]


def _text(fig, x, y, text, **kwargs):
    return fig.text(x / WIDTH_MM, 1 - y / HEIGHT_MM, text, **kwargs)


def _regime_key(ax, states):
    """Original decision regions and every simulated point, on the full original axes."""
    dv, di = REGIME_DV_PCT, REGIME_DI_PCT
    x0, x1, y0, y1 = -75, 22, -103, 13
    ax.add_patch(Rectangle((dv, y0), x1 - dv, y1 - y0, fc=REGIME_FILL[0], lw=0))
    ax.add_patch(Rectangle((x0, y0), dv - x0, y1 - y0, fc=REGIME_FILL[2], lw=0))
    ax.add_patch(Rectangle((x0, -di), dv - x0, 2 * di,
                           fc=REGIME_FILL[1], ec=REGIME_EDGE[1], lw=.5))
    ax.axvline(dv, color="0.35", lw=.55)
    ax.axhline(0, color="0.4", lw=.35)
    points = 0
    for name in ORDER:
        for shape, orientation in CONFIGS:
            sub = sel(states, name, shape, orientation)
            x, y = sub["Vmp_loss_%"].to_numpy(), sub["Imp_loss_%"].to_numpy()
            if not (np.isfinite(x).all() and np.isfinite(y).all() and
                    (x >= x0).all() and (x <= x1).all() and (y >= y0).all() and (y <= y1).all()):
                raise ValueError("Every saved response point must be finite and within the full map")
            filled = orientation == "landscape"
            line, = ax.plot(x, y, ls="none", marker=SHAPE_MARKER[shape], ms=2.3,
                            mec=MODULE_COLOR[name], mew=0 if filled else .45,
                            mfc=MODULE_COLOR[name] if filled else "none", alpha=.9, zorder=3)
            np.testing.assert_array_equal(line.get_xdata(), x)
            np.testing.assert_array_equal(line.get_ydata(), y)
            points += len(sub)
    ax.set(xlim=(x0, x1), ylim=(y0, y1),
           xlabel=r"$\Delta V_{\mathrm{mp}}/V_{\mathrm{mp,clean}}$ [%]",
           ylabel=r"$\Delta I_{\mathrm{mp}}/I_{\mathrm{mp,clean}}$ [%]")
    ax.set_xticks([-60, -40, -20, 0, 20])
    ax.set_yticks([-100, -80, -60, -40, -20, 0])
    ax.tick_params(pad=2, labelsize=6.3, length=2.5)
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    if points != len(states):
        raise ValueError("The response map must retain every saved point")
    return points


def _definitions(fig):
    """Regime colours and all four shape/orientation combinations for panel a."""
    for regime, title in enumerate(("I, current limited", "II, bypass clamped", "III, mixed")):
        y = 9.5 + 7 * regime
        left, bottom, width, height = _rect(113, y + .2, 3.3, 3.3)
        fig.add_artist(Rectangle((left, bottom), width, height, transform=fig.transFigure,
                                 fc=REGIME_FILL[regime], ec="none", lw=0))
        _text(fig, 119, y + .2 + 3.3 / 2, title, fontsize=7.5, fontweight="bold", va="center")
    _text(fig, 113, 30.5, "Band configuration", fontsize=7, fontweight="bold", va="top")
    for index, (shape, orientation) in enumerate(CONFIGS):
        x = 114 + 34 * (index % 2)
        y = 36.5 + 7 * (index // 2)
        fig.add_artist(Line2D([x / WIDTH_MM], [1 - y / HEIGHT_MM], transform=fig.transFigure,
                             ls="none", marker=SHAPE_MARKER[shape], ms=4.3, mew=.65,
                             mec="0.25", mfc="none" if orientation == "portrait" else "0.25"))
        _text(fig, x + 3, y, f"{shape.capitalize()}, {orientation}", fontsize=6.7, va="center")


def _contour_key(fig, centre_x):
    """Centre the line and its complete label together beneath the atlas."""
    fig.legend(handles=[Line2D([], [], color="0.15", lw=1.1)],
               labels=[r"Contours: $\Delta P_{\mathrm{mp}}/P_{\mathrm{mp,clean}}$ [%]"],
               loc="center", bbox_to_anchor=(centre_x / WIDTH_MM, 1 - 166.5 / HEIGHT_MM),
               bbox_transform=fig.transFigure, fontsize=7, frameon=False,
               borderpad=0, borderaxespad=0, handlelength=2.4, handletextpad=.7)


def _atlas_cell(ax, sub, last_row, first_column):
    grid, acovs, taus = plots._grid(sub, "regime")
    a, t = 100 * np.asarray(acovs), np.asarray(taus)
    mesh = ax.pcolormesh(plots._edges(t), plots._edges(a), np.ma.masked_invalid(grid),
                        cmap=ListedColormap([REGIME_FILL[k] for k in (0, 1, 2)]),
                        norm=Normalize(-0.5, 2.5), edgecolors="white", linewidth=0.4)
    np.testing.assert_array_equal(np.ma.filled(mesh.get_array(), np.nan).reshape(grid.shape), grid)
    loss, _, _ = plots._grid(sub, "Pmp_loss_%")
    levels = [level for level in CONTOURS if np.nanmin(loss) < level < np.nanmax(loss)]
    if levels:
        cs = ax.contour(t, a, loss, levels=levels, colors="0.15", linewidths=0.5,
                        linestyles="solid")
        for label in ax.clabel(cs, inline=True, inline_spacing=1, fontsize=5.3, fmt="%d"):
            label.set_path_effects([pe.withStroke(linewidth=1.4, foreground="white")])
    ax.set_xlim(plots._edges(t)[[0, -1]])
    ax.set_ylim(plots._edges(a)[[0, -1]])
    ax.set_xticks([0, 0.4, 0.8])
    ax.set_yticks([2, 6, 10])
    ax.tick_params(labelsize=6, pad=1.5, length=2,
                   labelleft=first_column, labelbottom=last_row)
    ax.set_ylabel(r"$A_{\mathrm{cov}}$ (%)" if first_column else "", fontsize=6.3, labelpad=.6)
    ax.set_xlabel(r"$\tau_I$" if last_row else "", fontsize=6.3, labelpad=.6)
    return int(np.isfinite(grid).sum())


def draw(states, path):
    """Data-first layout: full response map, concise key and self-keying atlas."""
    plotted_states = 0
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
        _text(fig, 3, 6, "a", fontsize=9, fontweight="bold", va="top")
        response_box = (41, 6, 62, 42.5)
        scatter_points = _regime_key(fig.add_axes(_rect(*response_box)), states)
        _definitions(fig)

        _text(fig, 3, 61, "b", fontsize=9, fontweight="bold", va="top")
        left, top, cell_width, cell_height, gap_x, gap_y = 18.5, 69, 38.5, 20.7, 2.5, 2
        # Module columns retain the layout of the original standalone atlas.
        # Colour frames around the single-line names also key the response scatter.
        for col, name in enumerate(ORDER):
            x = left + col * (cell_width + gap_x)
            centre = x + cell_width / 2
            _text(fig, centre, 63.5, f"{PRETTY[name]} ({TECH[name]})", ha="center", va="center",
                  fontsize=7.5, fontweight="bold",
                  bbox=dict(boxstyle="square,pad=0.3", facecolor="white",
                            edgecolor=MODULE_COLOR[name], linewidth=1.0))
            for row, (shape, orientation) in enumerate(CONFIGS):
                y = top + row * (cell_height + gap_y)
                ax = fig.add_axes(_rect(x, y, cell_width, cell_height))
                plotted_states += _atlas_cell(ax, sel(states, name, shape, orientation),
                                              row == len(CONFIGS) - 1, col == 0)
        # Vertical row labels leave room for wider maps. The same open/filled
        # rectangles and triangles are keyed explicitly beside panel a.
        for row, (shape, orientation) in enumerate(CONFIGS):
            y = top + row * (cell_height + gap_y) + cell_height / 2
            marker_height = 4.3 / 72 * 25.4
            label = _text(fig, 5.5, y - (1 + marker_height) / 2,
                          f"{shape.capitalize()}\n{orientation.capitalize()}", fontsize=6.7,
                          fontweight="bold", ha="center", va="center", rotation=90, multialignment="center")
            label_height = label.get_window_extent(fig.canvas.get_renderer()).height / fig.dpi * 25.4
            marker_y = y + label_height / 2 + .5
            fig.add_artist(Line2D([5.5 / WIDTH_MM], [1 - marker_y / HEIGHT_MM], transform=fig.transFigure,
                                  ls="none", marker=SHAPE_MARKER[shape], ms=4.3, mew=.65,
                                  mec="0.25", mfc="none" if orientation == "portrait" else "0.25"))
        contour_centre = left + (len(ORDER) * cell_width + (len(ORDER) - 1) * gap_x) / 2
        _contour_key(fig, contour_centre)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        assert len(fig.axes[1:]) == 16
        assert [ax.get_ylabel() for ax in fig.axes[1:]] == [r"$A_{\mathrm{cov}}$ (%)"] * 4 + [""] * 12
        assert [ax.get_xlabel() for ax in fig.axes[1:]] == ["", "", "", r"$\tau_I$"] * 4
        assert [bool(ax.get_yticklabels()) for ax in fig.axes[1:]] == [True] * 4 + [False] * 12
        assert [bool(ax.get_xticklabels()) for ax in fig.axes[1:]] == [False, False, False, True] * 4
        text_items = list(fig.texts)
        for ax in fig.axes:
            text_items += [*ax.texts, *ax.get_xticklabels(), *ax.get_yticklabels(), ax.xaxis.label, ax.yaxis.label]
        for legend in fig.legends:
            text_items += legend.get_texts()
        for text in text_items:
            if not text.get_visible() or not text.get_text():
                continue
            box = text.get_window_extent(renderer)
            if box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height:
                raise ValueError(f"Text outside final regime figure: {text.get_text()}")
        if plotted_states != len(states):
            raise ValueError("The atlas must include every saved state of the four modules")
        plots.save(fig, path, dpi=600, bbox_inches=None)
    return dict(width_mm=WIDTH_MM, height_mm=HEIGHT_MM, png_dpi=600,
                atlas_panels=16, saved_states_plotted=plotted_states,
                header="full simulation scatter with regime names and explicit configuration legend",
                scatter_points_plotted=scatter_points, scatter_points_clipped=0,
                scatter_limits_percent=dict(voltage=[-75, 22], current=[-103, 13]),
                atlas_column_modules=list(ORDER),
                atlas_row_configurations=[dict(shape=s, orientation=o) for s, o in CONFIGS],
                point_symbol_key="open/filled rectangles and triangles in panel a legend and atlas row labels",
                configuration_legend_labels=[f"{s.capitalize()}, {o}" for s, o in CONFIGS],
                atlas_row_label_rotation_degrees=90,
                classification_formulas_in_legend=False,
                regime_numerals_in_scatter=False,
                module_colour_key="coloured frames around single-line module and technology labels",
                atlas_configuration_symbols="inline with vertical row labels",
                panel_b_label_y_mm=61,
                coverage_labels_on_first_column_only=True,
                tau_labels_on_bottom_panel_of_each_column=True,
                atlas_tick_labels="left column and bottom row only",
                atlas_cell_size_mm=[cell_width, cell_height],
                atlas_gaps_mm=[gap_x, gap_y],
                response_map_box_mm=list(response_box),
                module_and_technology_font=dict(size_pt=7.5, weight="bold"),
                contour_key_font_pt=7,
                section_titles_removed=True, contour_key="centred below panel b",
                contour_key_y_mm=166.5, contour_key_centre_x_mm=contour_centre,
                voltage_boundary_percent=REGIME_DV_PCT, current_tolerance_percent=REGIME_DI_PCT,
                palette={str(k): REGIME_FILL[k] for k in (0, 1, 2)},
                source_metrics_recalculated=False)
