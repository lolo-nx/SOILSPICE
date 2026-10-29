"""Figures of chapter 4 (energy yield and cleaning costs), in the paper style: fig_band_growth, fig_energy_maps,
fig_design_levers, fig_design_height, fig_sensitivity, fig_cleaning; plus the inverse diagnosis (fig_inverse_chain)
and the cycle-based SI maps (fig_pan_h_maps), using the explicitly selected event washing assumption."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.interpolate import RegularGridInterpolator

from matplotlib.patches import FancyBboxPatch, Rectangle

from . import config, modules
from .energy import AF, ROBOT, cycle_loss, washing_flags, washing_threshold, REFERENCE_TOTAL
from .style import (LETTER, LOSS_CMAP, LOSS_NORM, MODULE_COLOR, ORDER, ORIENTATIONS, PRETTY, SITES as SITES_STYLE,
                    STYLE, LETTER_KW, SITE_COLOR, bold, edges, loss_colorbar, module_label, save)

INK = "black"
SITE = SITE_COLOR


def band_growth(series, path, *, threshold):
    """Band growth updated at rain-event closures (three consecutive dry hours)."""
    summary = []
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(3, 1, figsize=(7.2, 4.6), sharex=True)
        fig.subplots_adjust(left=0.09, right=0.90, bottom=0.12, top=0.84, hspace=0.65)
        ymax = max(30.0, 10 * np.ceil(100 * series.L_band.max() / 10))
        for k, (site, color) in enumerate(SITE.items()):
            s = series[series.site == site]
            dates = pd.to_datetime(s.date)
            ax = axes[k]
            ax.fill_between(dates, 0, 100 * s.L_band, color=color, alpha=0.25, lw=0,
                            step="post", label="band power loss")
            ax.set_ylim(0, ymax)
            ax.set_ylabel("Band loss (%)")
            ax2 = ax.twinx()
            ax2.spines["right"].set_visible(True)
            ax2.plot(dates, s.tau, color=color, lw=0.8, drawstyle="steps-post", label="band transmittance")
            ax2.set_ylim(0, 1.05)
            ax2.set_ylabel(r"Band $\tau_I$")
            washing = washing_flags(s.rain, washing_total=threshold)
            transfer = (s.rain >= 1.0) & ~washing
            event_times = pd.to_datetime(s.end_time)
            ax.plot(event_times[washing], np.full(washing.sum(), ymax * 0.97), marker="v", ls="none", ms=3,
                    color=INK, label="washing rain")
            n_years = dates.dt.year.nunique()
            ax.set_title(f"{site}: {washing.sum() / n_years:.1f} washing events per year, "
                         f"{transfer.sum() / n_years:.1f} transfer events per year",
                         loc="left", fontweight="normal")
            letter(ax, k)
            if k == 0:
                h1, l1 = ax.get_legend_handles_labels()
                h2, l2 = ax2.get_legend_handles_labels()
                fig.legend(h1 + h2, l1 + l2, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.99),
                           columnspacing=1.0)
            summary.append(dict(site=site, years=n_years, washing_events=int(washing.sum()),
                                transfer_events=int(transfer.sum()), washing_per_year=washing.sum() / n_years,
                                transfer_per_year=transfer.sum() / n_years, tau_min=s.tau.min(),
                                band_loss_mean_pct=100 * np.average(s.L_band, weights=s.Y)))
        if threshold != REFERENCE_TOTAL:
            condition = rf"$P_e\geq {threshold:g}$ mm/event"
            fig.text(0.5, 0.035,
                     "Assumed complete washing: " + condition + r"; otherwise transfer if $P_e\geq1$ mm."
                     "\n" + "3 dry hours separate events. Thresholds are not calibrated for bands.",
                     ha="center", fontsize=6)
        save(fig, path)
    return pd.DataFrame(summary)


def letter(ax, k, x=None, y=None, inside=True):
    if inside:
        ax.text(0.03, 0.97, "abcdefgh"[k], transform=ax.transAxes, va="top", **LETTER_KW)
    else:                                       # above the axes, as on the P_mp loss maps
        ax.text(-0.12, 1.02, "abcdefgh"[k], transform=ax.transAxes, va="bottom", **LETTER_KW)




def paper_figures(R, T, S, S16, H1, out, data, *, washing_total=REFERENCE_TOTAL):
    """All the chapter-4 figures of the main text and the Extended Data.
    R: model runs (energy.model_runs), T: maps on the f grid (energy.tau_grid), S: daily series of the reference
    module (energy.growth_series), S16: sensitivity (energy.sensitivity), H1: band height at 1 % loss
    (energy.design_h1). Figures in out, the best strategy per configuration in data. Returns the numbers quoted in
    the text (one row per item)."""
    selected_threshold = washing_threshold(washing_total)
    with plt.rc_context(STYLE):
        R = R.copy(); R["k"] = R.k.astype(int)
        KEY = ["site", "module", "orientation", "h", "f", "margin"]
        G = R[R.market == "ground"]
        man = G[G["mode"] != "robot"]
        B = man.loc[man.groupby(KEY).cost.idxmin()].set_index(KEY)
        B["cost_u_best"] = man.groupby(KEY).cost_u.min()
        B = B.reset_index()
        B["bc"] = 100 * (B.cost - B.cost_u_best) / B.revenue
        B["bc_eur"] = B.cost - B.cost_u_best
        B.to_csv(data / "best_strategy.csv", index=False)
        summary = []

        # ---- band growth at rain-event closures ------------------------------------------------------------------------
        growth = band_growth(S, out / "fig_band_growth", threshold=selected_threshold)
        for r in growth.itertuples(index=False):
            summary.append(dict(item="series", site=r.site, washing_per_yr=r.washing_per_year,
                                transfer_events_per_yr=r.transfer_per_year, tau_min=r.tau_min, Lband_mean=r.band_loss_mean_pct))

        # ---- energy maps against h and f ---------------------------------------------------------------------------------
        T = T.copy()
        T["band"] = 100 * (T.loss - T.loss_u)                    # minus the uniform layer alone (same module)
        med = T.groupby(["site", "h", "f"])[["band", "tau_mean"]].median().reset_index()
        fig = plt.figure(figsize=(7.2, 5.4))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.5, wspace=0.16)
        axs = np.array([fig.add_subplot(gs[0, 0])] + [None, None])
        axs[1] = fig.add_subplot(gs[0, 1], sharey=axs[0]); axs[2] = fig.add_subplot(gs[0, 2], sharey=axs[0])
        bxs = [fig.add_subplot(gs[1, 0])]; bxs += [fig.add_subplot(gs[1, j], sharey=bxs[0]) for j in (1, 2)]
        levels = [0, 0.5, 1, 2, 5, 10, 15, 20, 30, 100]
        cmap = plt.get_cmap("YlOrBr")
        cols = ["#ffffff"] + [cmap(v) for v in np.linspace(0.12, 0.95, len(levels) - 2)]   # white below 0.5 %, as in the paper
        ecmap, enorm = ListedColormap(cols), BoundaryNorm(levels, len(cols))


        def crossing(grid, axis_vals, level):
            """First position along a 1-D profile where it crosses level (None if it does not)."""
            d = grid - level; i = np.flatnonzero(np.sign(d[:-1]) != np.sign(d[1:]))
            return None if not len(i) else float(np.interp(0, [d[i[0]], d[i[0] + 1]], axis_vals[i[0]:i[0] + 2])) \
                if d[i[0]] < d[i[0] + 1] else float(np.interp(0, [d[i[0] + 1], d[i[0]]], axis_vals[i[0]:i[0] + 2][::-1]))


        F0 = 0.1                                                                         # reference f, row 2
        REF = ("JASOLAR", "landscape")   # reference configuration of the maps: closest to the median of the 8 over all (site, h, f)
        ref = T[(T.module == REF[0]) & (T.orientation == REF[1])][["site", "h", "f", "band", "tau_mean"]]
        for k, site in enumerate(SITE):
            ax = axs[k]; q = ref[ref.site == site]
            Hs = np.sort(q.h.unique()); Fs = np.sort(q.f.unique())
            Z = q.pivot(index="h", columns="f", values="band").loc[Hs, Fs].to_numpy()
            Tm = q.pivot(index="h", columns="f", values="tau_mean").loc[Hs, Fs].to_numpy()
            lf = np.log10(Fs)
            xi = np.linspace(lf[0], lf[-1], 200); yi = np.linspace(Hs[0], Hs[-1], 200); XX, YY = np.meshgrid(xi, yi)
            zi = RegularGridInterpolator((Hs, lf), Z)((YY, XX)); ti = RegularGridInterpolator((Hs, lf), Tm)((YY, XX))
            im = ax.contourf(XX, YY, zi, levels=levels, cmap=ecmap, norm=enorm)        # linear interpolation in (h, log f)
            # labels placed by hand so that they do not overlap: tau lines labelled near the top, loss lines on the right
            cs = ax.contour(XX, YY, ti, levels=[0.2, 0.5, 0.8], colors="#1f4e79", linewidths=0.7, linestyles="--")
            row = int(np.argmin(np.abs(yi - 90)))
            man = []; y = None
            for v in cs.levels:                          # near the top, or on the right edge if the line leaves earlier
                x = crossing(ti[row], xi, v)
                if x is None:
                    y = crossing(ti[:, -8], yi, v); x = xi[-8] if y is not None else None
                man += [] if x is None else [(x, 90 if y is None or crossing(ti[row], xi, v) is not None else y)]
                if man and man[-1][0] < xi[0] + 0.12:      # too close to the left edge: label lower, inside the axes
                    yl = crossing(ti[:, 12], yi, v)
                    if yl is not None and yl > yi[0] + 5:
                        man[-1] = (xi[12], yl)
                    else:
                        man.pop()                          # line along the edge: left unlabelled (same level as next panel)
                y = None
            ax.clabel(cs, fmt=lambda v: f"$\\bar\\tau$={v:.1f}", fontsize=5.5, inline_spacing=1, manual=man)
            c2 = ax.contour(XX, YY, zi, levels=[1, 5], colors=INK, linewidths=0.6)
            col = int(np.argmin(np.abs(xi - np.log10(0.6))))
            man = [(np.log10(0.6), y) for v in c2.levels if (y := crossing(zi[:, col], yi, v)) is not None]
            ax.clabel(c2, fmt="%g %%", fontsize=5.5, inline_spacing=1, manual=man)
            FF, HH = np.meshgrid(lf, Hs); ax.scatter(FF, HH, s=0.6, color=INK, alpha=0.5, lw=0)    # simulated points
            ax.axvline(np.log10(F0), color=INK, lw=0.5, ls=":")
            ft = [0.01, 0.03, 0.1, 0.3, 1.0]; ax.set_xticks(np.log10(ft)); ax.set_xticklabels([f"{v:g}" for v in ft])
            ax.set_title(site); ax.set_xlabel("Trapped fraction $f$")
            if k == 0:
                ax.set_ylabel("Band height $h$ (mm)")
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            letter(ax, k, inside=False)
            # row 2: each module and orientation at f = F0 (the dotted line of the map)
            bx = bxs[k]
            for mname in ORDER:
                for o, ls in (("portrait", "-"), ("landscape", "--")):
                    z = T[(T.site == site) & (T.module == mname) & (T.orientation == o) & np.isclose(T.f, F0)].sort_values("h")
                    isref = (mname, o) == REF
                    bx.plot(z.h, z.band, color=MODULE_COLOR[mname], ls=ls, lw=2.0 if isref else 0.9, zorder=3 if isref else 2)
            bx.set_xlim(10, 100); bx.set_xlabel("Band height $h$ (mm)")
            bx.set_title(site)
            if k == 0:
                bx.set_ylabel("Energy lost to the band\n(% per year)")
            else:
                plt.setp(bx.get_yticklabels(), visible=False)
            letter(bx, k + 3)
        bxs[0].set_ylim(bottom=0)
        from matplotlib.lines import Line2D  # noqa: E402
        hd = [Line2D([], [], color=MODULE_COLOR[mname], lw=1.0, label=PRETTY[mname]) for mname in ORDER]
        hd += [Line2D([], [], color="0.3", ls="-", label="portrait"), Line2D([], [], color="0.3", ls="--", label="landscape"),
               Line2D([], [], color=MODULE_COLOR[REF[0]], ls="--", lw=2.0, label="reference (maps)")]
        fig.legend(handles=hd, loc="upper center", ncol=7, bbox_to_anchor=(0.5, 0.475), handlelength=1.6, columnspacing=1.0)
        cb = fig.colorbar(im, ax=list(axs), orientation="horizontal", location="top", fraction=0.05, pad=0.13, aspect=60,
                          spacing="uniform", ticks=levels[:-1])
        cb.ax.set_xticklabels([f"{v:g}" for v in levels[:-1]])
        cb.set_label("Energy lost to the band (% per year)", labelpad=2)
        cb.outline.set_linewidth(0.6)
        save(fig, out / "fig_energy_maps")
        for site in SITE:
            for h in (40.0, 60.0, 80.0):
                q = T[(T.site == site) & (T.h == h)]
                for f in (0.03, 0.1, 0.3):
                    z = q[np.isclose(q.f, f)]
                    summary.append(dict(item="energy", site=site, h=h, f=f, band_med=z.band.median(), band_min=z.band.min(),
                                        band_max=z.band.max(), tau_med=z.tau_mean.median(), tau_min=z.tau_mean.min(),
                                        tau_max=z.tau_mean.max()))

        # ---- design levers ------------------------------------------------------------------------------------------------
        p = B.pivot_table(index=["site", "module", "orientation", "h", "f"], columns="margin", values="bc_eur")
        Lslope = np.array([modules.MODULES[m].height if o == "portrait" else modules.MODULES[m].width for (_, m, o, _, _) in p.index])
        V = pd.DataFrame({"val10": AF * (p[0.0] - p[10.0]), "glass": 100 * 10 / Lslope}).reset_index()
        net = pd.DataFrame({m: AF * (p[0.0] - p[m]) - 100 * m / Lslope for m in [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]})
        V["m_opt"] = net.idxmax(axis=1).to_numpy()
        fig = plt.figure(figsize=(7.2, 2.5))
        outer = fig.add_gridspec(1, 2, width_ratios=[1.15, 3.2], wspace=0.32)
        inner = outer[1].subgridspec(1, 3, wspace=0.1)
        axs = [fig.add_subplot(outer[0])]
        for i in range(3):
            axs.append(fig.add_subplot(inner[i], sharey=axs[1] if i else None))
        ax = axs[0]
        for site, col in SITE.items():
            q = V[V.site == site]
            m = q[q.f == 0.1].groupby("h").val10.median(); lo = q[q.f == 0.03].groupby("h").val10.median()
            hi = q[q.f == 0.3].groupby("h").val10.median()
            ax.fill_between(m.index, lo, hi, color=col, alpha=0.15, lw=0)
            ax.plot(m.index, m, color=col, label=site)
        from matplotlib.patches import Patch  # noqa: E402
        ax.plot([], [], color="0.4", label="$f$ = 0.1")
        ax.add_artist(ax.legend(loc="upper right", fontsize=6))
        ax.legend(handles=[Patch(color="0.4", alpha=0.25, lw=0, label="$f$ = 0.03–0.3")], loc="center right", fontsize=6)
        ax.set_xlim(10, 100); ax.set_ylim(0, 16)
        ax.set_xlabel("Band height $h$ (mm)"); ax.set_ylabel("Value of +10 mm margin\nover 25 years (€ kWp$^{-1}$)")
        letter(ax, 0, x=-0.24)
        O = B[B.margin == 0].pivot_table(index=["site", "module", "h", "f"], columns="orientation", values="bc").reset_index()
        O["d"] = O.portrait - O.landscape
        for k, site in enumerate(SITE):
            ax = axs[k + 1]
            for mname in ORDER:
                q = O[(O.site == site) & (O.module == mname)]
                m = q[q.f == 0.1].set_index("h").d
                ax.plot(m.index, m, color=MODULE_COLOR[mname], label=PRETTY[mname])
            ax.axhline(0, color=INK, lw=0.5)
            ax.set_xlim(10, 100); ax.set_ylim(-0.7, 0.9); ax.set_title(site, color=SITE[site])
            ax.set_xlabel("Band height $h$ (mm)")
            ax.text(98, 0.85, "landscape cheaper", fontsize=5.5, va="top", ha="right", color="0.35")
            ax.text(98, -0.65, "portrait cheaper", fontsize=5.5, va="bottom", ha="right", color="0.35")
            if k == 0:
                ax.set_ylabel("Portrait $-$ landscape\nband cost (pp of revenue)")
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            if k == 2:
                fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=4, bbox_to_anchor=(0.62, 0.98), handlelength=1.2)
            letter(ax, k + 1, x=-0.33 if k == 0 else -0.08)
        save(fig, out / "fig_design_levers")
        for site in SITE:
            for h in (40.0, 60.0, 80.0):
                q = V[(V.site == site) & (V.h == h)]
                r = dict(item="margin", site=site, h=h, glass_min=V.glass.min(), glass_max=V.glass.max())
                for f in (0.03, 0.1, 0.3):
                    z = q[np.isclose(q.f, f)]; r[f"val10_f{f}"] = z.val10.median(); r[f"val10min_f{f}"] = z.val10.min()
                    r[f"val10max_f{f}"] = z.val10.max(); r[f"mopt_f{f}"] = z.m_opt.median()
                summary.append(r)
                for mname in ORDER:
                    z = O[(O.site == site) & (O.module == mname) & (O.h == h)]
                    summary.append(dict(item="orientation", site=site, h=h, module=mname,
                                        **{f"d_f{f}": z[np.isclose(z.f, f)].d.iloc[0] for f in (0.03, 0.1, 0.3)}))
            for mname in ORDER:                                                   # height where the order reverses
                z = O[(O.site == site) & (O.module == mname) & (O.f == 0.1)].sort_values("h")
                hh, dd = z.h.to_numpy(), z.d.to_numpy(); i = np.flatnonzero((dd[:-1] < 0) & (dd[1:] >= 0))
                summary.append(dict(item="crossover", site=site, module=mname,
                                    h_cross=np.interp(0, dd[i[0]:i[0] + 2], hh[i[0]:i[0] + 2]) if len(i) else np.nan))

        # ---- lever 1: band height at 1 % loss, rectangle (a) and triangle (b) --------------------------------------------
        fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), sharey=True, gridspec_kw=dict(wspace=0.12))
        for k, (ax, shape) in enumerate(zip(axs, ("rectangle", "triangle"))):
            for (m, o), g in H1[H1["shape"] == shape].groupby(["module", "orientation"]):
                g = g.sort_values("tau_I"); c = MODULE_COLOR[m]; por = o == "portrait"
                ax.plot(g.tau_I, g.h1, color=c, ls="-" if por else "--", marker="o", ms=2.5, mfc=c if por else "white",
                        mew=0.7)
                ax.plot([-0.012], [g.d.iloc[0]], marker=">", ms=3.5, color=c, mfc=c if por else "white", mew=0.7,
                        clip_on=False)
            ax.set_xlim(0, 0.9); ax.set_ylim(0, 100)
            ax.set_xlabel("Band transmittance $\\tau_I$")
            ax.set_title("Rectangular band" if shape == "rectangle" else "Triangular band")
            letter(ax, k, inside=False)
        axs[0].set_ylabel("Band height at 1 % loss\n(mm, from the frame edge)")
        hd = [Line2D([], [], color=MODULE_COLOR[m], label=PRETTY[m]) for m in ORDER]
        hd += [Line2D([], [], color="0.4", ls="-", marker="o", ms=2.5, label="portrait"),
               Line2D([], [], color="0.4", ls="--", marker="o", ms=2.5, mfc="white", label="landscape"),
               Line2D([], [], color="0.4", ls="none", marker=">", ms=3.5, label="frame-to-cell distance $d$")]
        fig.legend(handles=hd, loc="lower center", ncol=7, bbox_to_anchor=(0.5, 0.98))
        save(fig, out / "fig_design_height")

        # ---- sensitivity (Extended Data), from econ_new_massband_sens.py -------------------------------------------------
        PAR = [("f", "Trapped fraction $f$", "0.03 / 0.3"),
               ("rate", "Soiling rate", "×0.5 / ×2"), ("alpha", "Loss per dust mass", "×0.75 / ×1.25"),
               ("light", "Event forming the band", "0.5 / 2 mm"), ("cost", "Cleaning cost", "×0.5 / ×2")]
        fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.2), sharey=True, gridspec_kw=dict(hspace=0.6, wspace=0.1))
        for r, (col, lab) in enumerate([("band_loss", "Energy lost to the band, rain only (% per year)"),
                                         ("band_cost", "Extra cost with cleaning when it pays (% of revenue)")]):
            for k, s_ in enumerate(SITE):
                ax = axs[r, k]; q = S16[S16.site == s_].set_index("case")[col]; ref = q["reference"]
                for i_, (p_, name, rng) in enumerate(PAR):
                    ax.barh(i_, q[f"{p_} low"] - ref, left=ref, height=0.6, color="#5aa3b5", **{"edgecolor": "black", "lw": 0.3})
                    ax.barh(i_, q[f"{p_} high"] - ref, left=ref, height=0.6, color="#c04a2b", **{"edgecolor": "black", "lw": 0.3})
                ax.axvline(ref, color=INK, lw=0.8)
                ax.set_yticks(range(len(PAR))); ax.set_yticklabels([f"{n} ({g})" for _, n, g in PAR], fontsize=6.5)
                ax.invert_yaxis(); ax.set_xlim(0, None); ax.tick_params(axis="y", length=0)
                if r == 0:
                    ax.set_title(s_)
                if k == 1:
                    ax.set_xlabel(lab)
                letter(ax, 3 * r + k, inside=False)
        fig.legend(handles=[Patch(facecolor="#5aa3b5", edgecolor="black", lw=0.3, label="low value"),
                            Patch(facecolor="#c04a2b", edgecolor="black", lw=0.3, label="high value"),
                            Line2D([], [], color=INK, lw=0.8, label="reference")],
                   loc="lower center", ncol=3, bbox_to_anchor=(0.55, 0.98))
        fig.text(0.55, -0.045, f"Scenario {selected_threshold:g} " + "mm/event" + "; complete washing assumed; threshold not calibrated for bands.",
                 ha="center", va="top", fontsize=6)
        save(fig, out / "fig_sensitivity")

        # ---- cleaning -------------------------------------------------------------------------------------------------------
        N = G[(G["mode"] == "none") & (G.margin == 0)].assign(bc=lambda x: 100 * (x.cost - x.cost_u) / x.revenue)
        Tg = G[(G["mode"] == "trigger") & (G.margin == 0)]
        fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(wspace=0.4, width_ratios=[1, 1, 1.3]))
        ax = axs[0]
        for site, col in SITE.items():
            n = N[(N.site == site) & (N.f == 0.1)].groupby("h").bc.median()
            b = B[(B.site == site) & (B.f == 0.1) & (B.margin == 0)].groupby("h").bc
            ax.plot(n.index, n, color=col, ls="--", lw=0.9)
            ax.fill_between(b.min().index, b.min(), b.max(), color=col, alpha=0.18, lw=0)
            ax.plot(b.median().index, b.median(), color=col, label=site)
        ax.set_ylim(0, 20); ax.set_xlim(10, 100)
        ax.set_xlabel("Band height $h$ (mm)"); ax.set_ylabel("Band cost (% of revenue)")
        ax.legend(handles=[Line2D([], [], color="0.4", ls="--", label="rain only"), Line2D([], [], color="0.4", label="cheapest cleaning")],
                  loc="upper left", bbox_to_anchor=(0.0, 0.88), fontsize=6)
        letter(ax, 0, x=-0.28)
        ax = axs[1]
        for site, col in SITE.items():
            q = Tg[(Tg.site == site) & (Tg.f == 0.1)]
            c = q.groupby("h").cleans
            ax.fill_between(c.min().index, c.min(), c.max(), color=col, alpha=0.18, lw=0)
            ax.plot(c.median().index, c.median(), color=col)
            ax.axhline(q.cleans_u.median(), color=col, ls=":", lw=0.9)
        ax.plot([], [], color="0.4", label="with the band"); ax.plot([], [], color="0.4", ls=":", label="without band")
        ax.set_xlim(10, 100); ax.set_ylim(0, 15)
        ax.set_xlabel("Band height $h$ (mm)"); ax.set_ylabel("Cleanings per year (trigger rule)")
        ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.77), fontsize=6)
        letter(ax, 1, x=-0.25)
        ax = axs[2]
        h0 = 60.0
        strategies = [("Rain\nonly", lambda q: q[q["mode"] == "none"]),
                      ("Best fixed\ninterval", lambda q: q[q["mode"] == "interval"]),
                      ("Trigger\nrule", lambda q: q[q["mode"] == "trigger"]),
                      ("Robot", lambda q: q[q["mode"] == "robot"])]
        w = 0.26
        for i, (site, col) in enumerate(SITE.items()):
            q0 = G[(G.site == site) & (G.f == 0.1) & (G.h == h0) & (G.margin == 0)]
            for j, (lab, sel) in enumerate(strategies):
                q = sel(q0)
                if lab == "Robot":
                    # model_runs excludes robot investment; include it exactly once here, in both counterfactuals.
                    cb_ = 100 * (q.cost + np.mean(ROBOT)) / q.revenue; cu_ = 100 * (q.cost_u + np.mean(ROBOT)) / q.revenue
                    err = 100 * (ROBOT[1] - ROBOT[0]) / 2 / q.revenue.iloc[0]
                else:
                    per = q.groupby(["module", "orientation"])
                    cb_ = per.apply(lambda t: (100 * t.cost / t.revenue).min()); cu_ = per.apply(lambda t: (100 * t.cost_u / t.revenue).min())
                    err = None
                x = j + (i - 1) * w
                ax.bar(x, np.median(cb_), w * 0.9, color=col, label=site if j == 0 else None, yerr=err, error_kw=dict(lw=0.6, capsize=1.5))
                ax.plot([x - w * 0.4, x + w * 0.4], [np.median(cu_)] * 2, color=INK, lw=0.9, label="without band" if (i, j) == (0, 0) else None)
                summary.append(dict(item="strategy", site=site, strategy=lab.replace("\n", " "), cost_band=np.median(cb_),
                                    cost_noband=np.median(cu_)))
        ax.set_xticks(range(len(strategies))); ax.set_xticklabels([s for s, _ in strategies], fontsize=6); ax.tick_params(axis="x", length=0)
        ax.set_ylabel("Soiling cost (% of revenue)"); ax.set_ylim(0, 40)
        hh, ll = ax.get_legend_handles_labels()
        ax.legend(hh[:1] + [Line2D([], [], color=INK, marker="|", ls="none", ms=6, mew=0.6)], ll[:1] + ["robot cost range"],
                  loc="upper right", fontsize=6)
        letter(ax, 2, x=-0.17)
        fig.legend(handles=[Line2D([], [], color=c, lw=1.4, label=s_) for s_, c in SITE.items()]
                   + [Patch(color="0.4", alpha=0.3, lw=0, label="range over 4 modules × 2 orientations")],
                   loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.98))
        save(fig, out / "fig_cleaning")
        for site in SITE:
            for h in (40.0, 60.0, 80.0):
                for f in (0.03, 0.1, 0.3):
                    q = Tg[(Tg.site == site) & (Tg.h == h) & np.isclose(Tg.f, f)]
                    b = B[(B.site == site) & (B.h == h) & np.isclose(B.f, f) & (B.margin == 0)]
                    n = N[(N.site == site) & (N.h == h) & np.isclose(N.f, f)]
                    summary.append(dict(item="cleaning", site=site, h=h, f=f, cleans=q.cleans.median(), cleans_u=q.cleans_u.median(),
                                        bc_best=b.bc.median(), bc_best_min=b.bc.min(), bc_best_max=b.bc.max(), bc_none=n.bc.median(),
                                        best_modes="; ".join(sorted(set(b["mode"] + b.k.astype(str))))))
        rt = R[(R.market == "rooftop") & (R["mode"] == "trigger") & (R.margin == 0)]
        for site in SITE:
            for h in (40.0, 60.0, 80.0):
                q = rt[(rt.site == site) & (rt.h == h) & (rt.f == 0.1)]
                summary.append(dict(item="rooftop", site=site, h=h, f=0.1, cleans=q.cleans.median(), cleans_u=q.cleans_u.median(),
                                    bc_best=(100 * (q.cost - q.cost_u) / q.revenue).median()))
    return pd.DataFrame(summary)


# ---- inverse diagnosis --------------------------------------------------------------------------------------------------
M_INV, O_INV, H0_INV, T0_INV = "SUNPOWER", "landscape", 80, 0.7
BLUE, PURPLE, NAVY, BROWN, GREEN = "#1f77b4", "#7b3f99", "#1f3d6b", "#a0522d", "#2e7d32"


def arrow(fig, x0, x1, y, label=None):
    fig.add_artist(plt.Annotation("", xy=(x1, y), xytext=(x0, y), xycoords="figure fraction",
                                  arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.25")))
    if label:
        fig.text((x0 + x1) / 2, y + 0.035, label, ha="center", va="bottom", fontsize=7)


def bracket(fig, x0, x1, y, text, color):
    fig.add_artist(plt.Line2D([x0, x0, x1, x1], [y + 0.02, y, y, y + 0.02], transform=fig.transFigure,
                              color=color, lw=1.0))
    fig.text((x0 + x1) / 2, y - 0.02, text, ha="center", va="top", color=color, style="italic", fontsize=7)


def inverse_chain(path, iv_csv, e_cycle, CP):
    """From a measured I-V curve back to the band and the cleaning decision (SunPower landscape, 80 mm,
    tau_I = 0.7 as the unknown band). iv_csv: SPICE I-V curves used as a synthetic measurement; e_cycle: mean clean
    energy between two washing rains per site [kWh/kWp] (energy.washing_cycles); CP: one cleaning in kWh/kWp."""
    M, O, H0, T0 = M_INV, O_INV, H0_INV, T0_INV
    sw = pd.read_csv(config.EQUAL_HEIGHT_CSV)
    g = sw[(sw.module == M) & (sw.orientation == O) & (sw["shape"] == "rectangle")]
    H = np.sort(g.height_mm.unique()).astype(float); TAU = np.sort(g.tau_I.unique())
    P = g.pivot_table(index="height_mm", columns="tau_I", values="Pmp_loss_%").to_numpy()
    I = g.pivot_table(index="height_mm", columns="tau_I", values="Imp_loss_%").to_numpy()
    V = g.pivot_table(index="height_mm", columns="tau_I", values="Vmp_loss_%").to_numpy()
    i0, j0 = list(H).index(H0), list(TAU).index(T0)
    dP, dI, dV = P[i0, j0], I[i0, j0], V[i0, j0]
    D = np.hypot(I - dI, V - dV)
    ib, jb = np.unravel_index(np.argmin(D), D.shape)

    # optical loss S = A_cov (1 - tau) in % of the cell area (as in the paper); for a band within the bottom row,
    # the I-V depends on S only
    NROW = 100
    A = g.groupby("height_mm").A_cov.first().reindex(H).to_numpy()
    S = (NROW * A)[:, None] * (1 - TAU[None, :])
    fam = list(zip(*np.where(D < 1.0)))                              # bands that fit the signature
    S_est = np.mean([S[i, j] for i, j in fam])

    # every (h, tau) with the same S: tau = 1 - S / A_cov(h), from the lowest (opaque) band upwards
    FH = [(h, 1 - S_est / (NROW * a), i) for i, (h, a) in enumerate(zip(H, A)) if NROW * a >= S_est]
    saved_all = [cycle_loss(TAU, np.abs(P[i]), max(tf, 0.0), "opt") for _, tf, i in FH] + \
                [cycle_loss(TAU, np.abs(P[i]), TAU[j], "opt") for i, j in fam]      # incl. the exact grid bands
    s_lo, s_hi = min(saved_all), max(saved_all)
    hc = np.linspace(0, H[-1], 2000); fc = np.interp(hc, H, NROW * A)
    h_min = hc[np.argmax(fc >= S_est)]
    th = {n: 100 * CP / e_cycle[n] for n in ["Neuchâtel", "Almería", "Doha"]}


    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(14, 4.3))
        L = [0.03, 0.225, 0.42, 0.665, 0.855]; W = [0.15, 0.145, 0.15, 0.15, 0.135]; B, Ht = 0.33, 0.52
        ARR = [(0.183, 0.198), (0.373, 0.383), (0.622, 0.632), (0.818, 0.85)]

        # (a) I-V curves simulated with SPICE (ngspice, same netlist), used as a synthetic measurement
        a = fig.add_axes([L[0], B, W[0], Ht])
        ivc = pd.read_csv(iv_csv)
        for col, c, ls, lab in [("I_clean", BLUE, "-", "clean"),
                                ("I_h80_t0.7", PURPLE, "-", "soiled module")]:
            a.plot(ivc.V, ivc[col], color=c, lw=1.4 if ls == "-" else 1.6, ls=ls, label=lab)
            if ls == "-":
                k = np.argmax(ivc.V * ivc[col]); a.plot(ivc.V[k], ivc[col][k], "o", color=c, ms=4.5)
        a.set_xlim(0, 72); a.set_ylim(0, 7.2)
        a.set_xlabel("Voltage [V]"); a.set_ylabel("Current [A]")
        a.legend(loc="lower left", fontsize=6.3, frameon=False, handlelength=1.6)
        a.set_title("a   Measured I–V curve", loc="left", pad=8)

        # (b) signature
        b = fig.add_axes([L[1], B, W[1], Ht])
        names = [r"$\Delta P_\mathrm{mp}$", r"$\Delta I_\mathrm{mp}$", r"$\Delta V_\mathrm{mp}$"]
        vals = [dP, dI, dV]
        b.barh(names, vals, color=["0.55", PURPLE, BLUE], edgecolor="black", lw=0.4, height=0.55)
        for k, x in enumerate(vals):
            b.text(x + (0.6 if x >= 0 else -0.6), k, f"{x:+.1f} %", va="center", ha="left" if x >= 0 else "right", fontsize=7)
        b.axvline(0, color="black", lw=0.6); b.set_xlim(-20, 12); b.invert_yaxis()
        b.set_xlabel("Change vs clean [%]")
        b.text(0.03, 0.04, "bypass step at ~45 V,\nMPP above it", transform=b.transAxes, ha="left", va="bottom", fontsize=7, color="0.3")
        b.set_title("b   Signature", loc="left", pad=8)

        # (c) library search, with iso-S lines
        c = fig.add_axes([L[2], B, W[2], Ht])
        m = c.pcolormesh(np.r_[TAU - 0.05, TAU[-1] + 0.05], np.r_[H[0] - 2.5, (H[1:] + H[:-1]) / 2, H[-1] + 10], D,
                         cmap="Greys_r", vmin=0, vmax=20, shading="flat")
        hh = np.linspace(H[0] - 2.5, H[-1] + 10, 300); tt = np.linspace(-0.05, 0.95, 300)
        SS = np.interp(hh, H, NROW * A)[:, None] * (1 - tt[None, :])
        cs = c.contour(tt, hh, SS, levels=[1, 3, 5, 7], colors="#5fa8d3", linewidths=0.7, linestyles="--")
        c.clabel(cs, fmt="S=%.0f %%", fontsize=5.5, inline_spacing=2)
        c.contour(tt, hh, SS, levels=[S_est], colors=PURPLE, linewidths=1.3)
        ii, jj = np.array(fam).T
        c.plot(TAU[jj], H[ii], "o", mfc="none", mec=PURPLE, mew=1.2, ms=9, label=f"misfit < 1 %  (S ≈ {S_est:.1f} %)")
        c.plot(TAU[jb], H[ib], marker="*", color="#f2b705", mec="black", mew=0.6, ms=12)
        c.legend(loc="upper left", fontsize=6, handletextpad=0.3, borderpad=0.3)
        c.set_xlabel(r"$\tau$"); c.set_ylabel(r"$h$ [mm]"); c.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
        c.set_xlim(-0.05, 0.95); c.set_ylim(H[0] - 2.5, H[-1] + 10)
        cb = fig.colorbar(m, ax=c, fraction=0.05, pad=0.03); cb.set_label("Misfit [%]", fontsize=7); cb.ax.tick_params(labelsize=6)
        c.set_title("c   SPICE library search", loc="left", pad=8)

        # (d) the loss depends on S only
        d = fig.add_axes([L[3], B, W[3], Ht])
        k = A[:, None].repeat(len(TAU), 1) > 0
        d.scatter(S[k], np.abs(P[k]), s=9, c=np.repeat(H[:, None], len(TAU), 1)[k], cmap="Greys", vmin=0, vmax=110,
                  edgecolors="black", linewidths=0.3, zorder=2, label="all simulated bands")
        d.axvline(S_est, color=PURPLE, lw=1.0, ls="--"); d.axhline(-dP, color=PURPLE, lw=1.0, ls="--")
        d.plot([S[i, j] for i, j in fam], [abs(P[i, j]) for i, j in fam], "o", mfc="none", mec=PURPLE, mew=1.2, ms=8, zorder=3)
        d.text(S_est + 0.25, 1.0, f"S ≈ {S_est:.1f} %", color=PURPLE, fontsize=7)
        d.text(6.6, np.abs(P).max() - 1.2, "bypass cap", ha="center", va="top", fontsize=6.5, color="0.3")
        d.set_xlim(0, 9); d.set_ylim(0, None)
        d.set_xlabel(r"Optical loss $S = A_\mathrm{cov}(1-\tau)$ [%]"); d.set_ylabel(r"$P_\mathrm{mp}$ loss [%]")
        d.legend(loc="upper left", fontsize=6, handletextpad=0.2, borderpad=0.3, markerscale=1.3)
        d.set_title("d   What the I–V measures", loc="left", pad=8)

        # (e) all bands with this S, and the cleaning decision
        e = fig.add_axes([L[4], B + 0.3, W[4], Ht - 0.3])
        hh2 = np.linspace(h_min, H[-1], 200)
        e.plot(hh2, 1 - S_est / np.interp(hh2, H, NROW * A), color=PURPLE, lw=1.3)
        e.plot([h for h, _, _ in FH], [max(tf, 0) for tf, in [(x[1],) for x in FH]], "o", mfc="white", mec=PURPLE, ms=4)
        e.axvspan(0, h_min, color="0.85", lw=0)
        e.text(h_min / 2, 0.5, "no band\ngives S", ha="center", va="center", fontsize=6, color="0.35")
        e.set_xlim(0, H[-1] + 5); e.set_ylim(0, 1)
        e.set_xlabel(r"$h$ [mm]", labelpad=1); e.set_ylabel(r"$\tau$")
        e.set_title(f"e   All bands with S ≈ {S_est:.1f} %", loc="left", pad=8)
        e2 = fig.add_axes([L[4], B - 0.03, W[4], 0.2]); e2.set_axis_off(); e2.set_xlim(0, 1); e2.set_ylim(0, 1)
        e2.add_patch(FancyBboxPatch((0.0, 0.0), 1.0, 1.0, boxstyle="round,pad=0.01", fc="#f4f4f4", ec="0.6", lw=0.6))
        lines = [rf"h ≥ {h_min:.0f} mm; taller $\Rightarrow$ lighter",
                 f"Saved by one cleaning: {s_lo:.1f}–{s_hi:.1f} %"] + \
                [f"{n}: {'pays' if s_lo > t else ('does not pay' if s_hi < t else 'unsure')} (> {t:.1f} %)"
                 for n, t in th.items()]
        e2.text(0.04, 0.93, "\n".join(lines), va="top", fontsize=6.4, linespacing=1.45)

        for x0, x1 in ARR:
            arrow(fig, x0, x1, B + Ht / 2)
        CAP = ["IV tracer on a soiled module\n(SunPower, landscape) vs clean\nreference; the band is unknown",
               "relative change of $P_\\mathrm{mp}$, $I_\\mathrm{mp}$, $V_\\mathrm{mp}$\nand bypass step (yes / no)",
               "misfit to every simulated band;\n$\\star$ closest, ○ within 1 point:\nthey lie on one iso-S line",
               "all bands collapse on one curve:\nthe I–V gives S,\nnot h and τ separately",
               "τ = 1 − S / A_cov(h): every band on the curve\nfits the I–V; same decision for all"]
        for x, w, tx in zip(L, [W[0], W[1], W[2] + 0.03, W[3], W[4]], CAP):
            fig.text(x + w / 2, 0.205 if x != L[4] else 0.25, tx, ha="center", va="top", fontsize=7, linespacing=1.3)
        bracket(fig, L[0], L[0] + W[0], 0.07, "IV tracer (here: SPICE curve)", BROWN)
        bracket(fig, L[1], L[1] + W[1], 0.07, "Python", NAVY)
        bracket(fig, L[2], L[3] + W[3], 0.07, "SPICE library (LTspice, precomputed)", BROWN)
        bracket(fig, L[4], L[4] + W[4], 0.07, "Python + cleaning map", NAVY)
        save(fig, path)


# ---- SI: energy lost while a band grows in height (cycle-based model) --------------------------------------------------
def pan_h_maps(t, cp_ground, e_cycle, path):
    """1 - P_an,h over (tau_I, h_max) for every module, with the break-even contour of one cleaning at each site.
    t: energy.pan_h(); cp_ground: one cleaning in kWh/kWp; e_cycle: energy between two washing rains per site."""
    LINES = [(n, 100 * cp_ground / e_cycle[n], ls) for n, _, ls, _ in SITES_STYLE]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.0), sharex=True, sharey=True)
        for i, o in enumerate(ORIENTATIONS):
            for c, (ax, m) in enumerate(zip(axes[i], ORDER)):
                g = t[(t.module == m) & (t.orientation == o)].pivot_table(index="h_max_mm", columns="tau_I", values="cycle_loss_pct")
                tau, h, z = g.columns.to_numpy(), g.index.to_numpy(), g.to_numpy()
                ax.pcolormesh(tau, h, z, cmap=LOSS_CMAP, norm=LOSS_NORM, shading="nearest")
                for name, lvl, ls in LINES:
                    ax.contour(tau, h, z, levels=[lvl], colors="k", linewidths=0.9, linestyles=ls)
                ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
                k = "abcdefgh"[i * 4 + c]
                ax.text(-0.02, 1.02, f"{k}  {module_label(m)}" if i == 0 else k, transform=ax.transAxes, fontweight="bold",
                        va="bottom", fontsize=7)
                if i == 1:
                    ax.set_xlabel(r"$\tau_I$")
            axes[i, 0].set_ylabel(bold(o.capitalize()) + "\n" + r"$h_{max}$ [mm]", labelpad=2)
        loss_colorbar(fig, axes.ravel().tolist(), r"Energy lost over the cycle, band growing in height, $1 - P_{an,h}$ [%]")
        hd = [Line2D([], [], color="k", ls=ls, lw=0.9, label=f"Cleaning pays above: {n} ({lvl:.1f} %)") for n, lvl, ls in LINES]
        fig.legend(handles=hd, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handlelength=2.2)
        save(fig, path)



# ---- rectangle against triangle -------------------------------------------------------------------------------------------
def band_shape(T, path, f0=0.1):
    """Energy lost to the band (rain only) against its height, rectangular (solid) and triangular (dashed) band, at
    the three sites: median of the 8 configurations (line) and their range (shaded). T: energy.tau_grid of both
    shapes (column shape), at least f = f0."""
    T = T[np.isclose(T.f, f0)].assign(band=lambda x: 100 * (x.loss - x.loss_u))
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(wspace=0.28))   # own y scale per site
        for k, (ax, site) in enumerate(zip(axs, SITE)):
            for shape, ls in (("rectangle", "-"), ("triangle", "--")):
                q = T[(T.site == site) & (T["shape"] == shape)].groupby("h").band
                ax.fill_between(q.min().index, q.min(), q.max(), color=SITE[site], alpha=0.15 if ls == "-" else 0.08,
                                lw=0)
                ax.plot(q.median().index, q.median(), color=SITE[site], ls=ls, lw=1.1)
            ax.set_xlim(10, 100); ax.set_ylim(bottom=0); ax.set_title(site); ax.set_xlabel("Band height $h$ (mm)")
            letter(ax, k, inside=False)
        axs[0].set_ylabel("Energy lost to the band\n(% per year, rain only)")
        hd = [Line2D([], [], color="0.3", ls="-", label="rectangular band"),
              Line2D([], [], color="0.3", ls="--", label="triangular band"),
              Patch(color="0.3", alpha=0.2, lw=0, label="range over 4 modules × 2 orientations")]
        fig.legend(handles=hd, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.98))
        save(fig, path)
