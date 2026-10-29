"""Threshold sensitivity from saved chapter-4 runs, without changing the energy model."""
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from . import config, energy
from .style import MODULE_COLOR, ORDER, ORIENTATIONS, PRETTY, STYLE, save

THRESHOLDS = (5., 10., 15., 20., 25.)
KEY = ["site", "market", "module", "orientation", "h", "f", "margin"]
CONFIG = ["site", "market", "module", "orientation", "h_mm", "f"]
PAIR = ["site", "market", "module", "h_mm", "f"]
EPS_EUR = 1e-9
EPS_PCT = 1e-8
MARGIN_MM = 10.
MODULE_PRICE_EUR_KWP = 100.  # Existing figure convention: 0.10 EUR/Wp.
REF_HEIGHT, REF_FRACTION = 60., .1


def states(values, tolerance):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("nonfinite value in a conclusion")
    return np.where(values > tolerance, "positive", np.where(values < -tolerance, "negative", "zero"))


def winner(frame, cost, frequency):
    """Independent minimum for one counterfactual, with co-optimal frequency ranges."""
    group = frame.groupby(KEY, sort=False)
    selected = frame.loc[group[cost].idxmin()].set_index(KEY)
    tied = frame[np.abs(frame[cost] - group[cost].transform("min")) <= EPS_EUR]
    ties = tied.groupby(KEY).agg(cooptimal_policy_count=(cost, "size"),
                                cooptimal_frequency_min=(frequency, "min"),
                                cooptimal_frequency_max=(frequency, "max"))
    return selected.join(ties)


def summarize_runs(runs, threshold, *, slope_length=energy.slope_length):
    """All h/f configurations, both saved markets; reference margin 0 versus +10 mm."""
    r = runs[runs.margin.isin((0., MARGIN_MM)) & (runs["mode"] != "robot")].copy()
    if r.duplicated(KEY + ["mode", "k"]).any():
        raise ValueError("duplicate saved policy/configuration")
    if r.select_dtypes("number").isna().any().any():
        raise ValueError("missing saved economic inputs")
    for column in ("Y", "revenue"):
        if not np.allclose(r.groupby(KEY)[column].max(), r.groupby(KEY)[column].min(), atol=1e-10):
            raise ValueError(f"{column} changes between policies")
    none = r[r["mode"] == "none"].set_index(KEY)
    band = winner(r, "cost", "cleans")
    uniform = winner(r, "cost_u", "cleans_u")
    manual = winner(r[r["mode"] != "none"], "cost", "cleans")
    manual_u = winner(r[r["mode"] != "none"], "cost_u", "cleans_u")
    if not all(frame.index.equals(band.index) for frame in (uniform, manual, manual_u)):
        raise ValueError("counterfactual configuration grids differ")
    none = none.reindex(band.index)
    out = pd.DataFrame(index=band.index)
    out["R_H_mm_event"] = threshold
    out["clean_energy_kwh_kwp_year"] = none.Y
    out["clean_revenue_eur_kwp_year"] = none.revenue
    out["energy_none_total_loss_pct"] = 100 * none.loss
    out["energy_none_uniform_loss_pct"] = 100 * none.loss_u
    out["energy_none_band_extra_loss_pct"] = 100 * (none.loss - none.loss_u)
    out["energy_none_band_extra_kwh_kwp_year"] = (none.loss - none.loss_u) * none.Y
    for prefix, frame, c, loss, frequency in (
        ("none_total", none, "cost", "loss", "cleans"),
        ("none_uniform", none, "cost_u", "loss_u", "cleans_u"),
        ("opt_total", band, "cost", "loss", "cleans"),
        ("opt_uniform", uniform, "cost_u", "loss_u", "cleans_u"),
        ("manual_total", manual, "cost", "loss", "cleans"),
        ("manual_uniform", manual_u, "cost_u", "loss_u", "cleans_u"),
    ):
        out[f"cost_{prefix}_eur_kwp_year"] = frame[c]
        out[f"cost_{prefix}_pct_clean_revenue"] = 100 * frame[c] / none.revenue
        if not prefix.startswith("none"):
            out[f"{prefix}_policy"] = frame["mode"]
            out[f"{prefix}_interval_days"] = frame.k
            out[f"{prefix}_cleanings_per_year"] = frame[frequency]
            out[f"{prefix}_energy_loss_pct"] = 100 * frame[loss]
            for col in ("cooptimal_policy_count", "cooptimal_frequency_min", "cooptimal_frequency_max"):
                out[f"{prefix}_{col}"] = frame[col]
    out["cost_none_band_extra_eur_kwp_year"] = none.cost - none.cost_u
    out["cost_opt_band_extra_eur_kwp_year"] = band.cost - uniform.cost_u
    for mode in ("none", "opt"):
        out[f"cost_{mode}_band_extra_pct_clean_revenue"] = 100 * out[f"cost_{mode}_band_extra_eur_kwp_year"] / none.revenue
    out["manual_gain_eur_kwp_year"] = none.cost - manual.cost
    out["manual_uniform_gain_eur_kwp_year"] = none.cost_u - manual_u.cost_u
    out["selected_gain_eur_kwp_year"] = none.cost - band.cost
    out["band_extra_cost_reduction_by_optimization_eur_kwp_year"] = out.cost_none_band_extra_eur_kwp_year - out.cost_opt_band_extra_eur_kwp_year
    for col in ("manual_gain", "manual_uniform_gain", "selected_gain", "band_extra_cost_reduction_by_optimization"):
        out[f"{col}_pct_clean_revenue"] = 100 * out[f"{col}_eur_kwp_year"] / none.revenue
    # Uniform and band policies are optimised again independently at each margin.
    base = out.xs(0., level="margin").copy()
    extra = out.xs(MARGIN_MM, level="margin").reindex(base.index)
    np.testing.assert_allclose(base.cost_opt_uniform_eur_kwp_year, extra.cost_opt_uniform_eur_kwp_year, atol=1e-10)
    lengths = np.array([slope_length(m, o) for _, _, m, o, _, _ in base.index])
    if (lengths <= 0).any():
        raise ValueError("invalid slope length")
    base["margin10_capex_eur_kwp"] = MODULE_PRICE_EUR_KWP * MARGIN_MM / lengths
    base["margin10_savings_eur_kwp_year"] = base.cost_opt_band_extra_eur_kwp_year - extra.cost_opt_band_extra_eur_kwp_year
    base["margin10_annualized_capex_eur_kwp_year"] = base.margin10_capex_eur_kwp / energy.AF
    base["margin10_net_eur_kwp_year"] = base.margin10_savings_eur_kwp_year - base.margin10_annualized_capex_eur_kwp_year
    base["margin10_net_pct_clean_revenue"] = 100 * base.margin10_net_eur_kwp_year / base.clean_revenue_eur_kwp_year
    base["margin10_npv25_eur_kwp"] = energy.AF * base.margin10_savings_eur_kwp_year - base.margin10_capex_eur_kwp
    base["margin10_opt_total_cost_eur_kwp_year"] = extra.cost_opt_total_eur_kwp_year
    base["margin10_opt_policy"] = extra.opt_total_policy
    base["margin10_opt_interval_days"] = extra.opt_total_interval_days
    base["margin10_opt_cleanings_per_year"] = extra.opt_total_cleanings_per_year
    base["margin10_state"] = states(base.margin10_net_eur_kwp_year, EPS_EUR)
    base["manual_gain_state"] = states(base.manual_gain_eur_kwp_year, EPS_EUR)
    base["selected_policy_uses_manual_cleaning"] = base.opt_total_cleanings_per_year > EPS_EUR
    base["shape"] = "rectangle"
    base["base_margin_mm"] = 0.
    return base.reset_index().rename(columns={"h": "h_mm"})


def washing_summary(intervals):
    """Use complete inter-washing gaps across years; retain both censored boundaries."""
    intervals = intervals.rename(columns={"P_H_mm": "R_H_mm_event"}).copy()
    for col in ("start", "end"):
        intervals[col] = pd.to_datetime(intervals[col], utc=True)
    rows = []
    for (site, threshold), q in intervals.groupby(["site", "R_H_mm_event"], sort=False):
        q = q.sort_values("start")
        if q.left_censored.sum() != 1 or q.right_censored.sum() != 1 or not q.left_censored.iloc[0] or not q.right_censored.iloc[-1]:
            raise ValueError("expected exactly one initial and final censored interval")
        if not np.array_equal(q.end.iloc[:-1].to_numpy(), q.start.iloc[1:].to_numpy()):
            raise ValueError("washing interval gaps/overlaps")
        start, stop = q.start.iloc[0], q.end.iloc[-1]
        nominal_stop = pd.Timestamp("2026-01-01", tz="UTC")
        # Correctly labelled source records end at 23:00 on the last day; do not invent the missing hour.
        if start != pd.Timestamp("2016-01-01", tz="UTC") or stop not in (nominal_stop, nominal_stop - pd.Timedelta(hours=1)):
            raise ValueError("expected the 2016–2025 period, allowing the final missing source hour")
        np.testing.assert_allclose(q.days.sum(), (stop - start).total_seconds() / 86400)
        complete = q[~q.left_censored & ~q.right_censored]
        longest = complete.loc[complete.days.idxmax()] if len(complete) else None
        left, right = q[q.left_censored].iloc[0], q[q.right_censored].iloc[0]
        rows.append(dict(site=site, R_H_mm_event=threshold, washing_events=len(q) - 1,
                         washing_events_per_calendar_year=(len(q) - 1) / 10,
                         max_complete_washing_gap_days=longest.days if longest is not None else np.nan,
                         max_complete_gap_start=longest.start if longest is not None else pd.NaT,
                         max_complete_gap_end=longest.end if longest is not None else pd.NaT,
                         initial_censored_start=left.start, initial_censored_end=left.end, initial_censored_days=left.days,
                         final_censored_start=right.start, final_censored_end=right.end, final_censored_days=right.days,
                         max_observed_including_censored_days=q.days.max(),
                         economic_years_365p25=q.days.sum() / 365.25, calendar_year_count=10,
                         study_start_utc=start, study_stop_utc=stop,
                         missing_trailing_source_hours=(nominal_stop - stop) / pd.Timedelta(hours=1)))
    return pd.DataFrame(rows), intervals


def orientation_comparison(summary):
    keys = PAIR + ["R_H_mm_event"]
    p = summary[summary.orientation == "portrait"].set_index(keys)
    l = summary[summary.orientation == "landscape"].set_index(keys).reindex(p.index)
    rows = pd.DataFrame(index=p.index)
    for label, column in (("energy", "energy_none_band_extra_loss_pct"),
                          ("cost", "cost_opt_band_extra_pct_clean_revenue")):
        rows[f"{label}_portrait_minus_landscape_pp"] = p[column] - l[column]
        sign = states(rows[f"{label}_portrait_minus_landscape_pp"], EPS_PCT)
        rows[f"{label}_preferred_orientation"] = np.where(sign == "positive", "landscape", np.where(sign == "negative", "portrait", "tie"))
    return rows.reset_index()


def classify_across_thresholds(summary, orientations):
    """Preserve signs/ties at every threshold; report discrete adjacent brackets, no interpolation."""
    robust, transitions = [], []
    for identifiers, q in summary.groupby(CONFIG, sort=False):
        q = q.sort_values("R_H_mm_event")
        if tuple(q.R_H_mm_event) != THRESHOLDS:
            raise ValueError("robustness requires all five thresholds exactly once")
        ids = dict(zip(CONFIG, identifiers))
        row = dict(**ids)
        for label, column, tolerance in (("margin10", "margin10_net_eur_kwp_year", EPS_EUR),
                                          ("manual_cleaning", "manual_gain_eur_kwp_year", EPS_EUR)):
            values = q[column].to_numpy()
            labels = states(values, tolerance)
            row[f"{label}_all_profitable"] = bool((labels == "positive").all())
            row[f"{label}_never_profitable"] = bool((labels != "positive").all())
            row[f"{label}_sign_stable"] = len(set(labels)) == 1
            row[f"{label}_states_5_10_15_20_25"] = ";".join(labels)
            row[f"{label}_min_eur_kwp_year"], row[f"{label}_max_eur_kwp_year"] = values.min(), values.max()
            row[f"{label}_profitable_thresholds_mm"] = ";".join(f"{t:g}" for t, s in zip(THRESHOLDS, labels) if s == "positive")
            for j in range(len(q) - 1):
                if labels[j] != labels[j + 1]:
                    transitions.append(dict(**ids, conclusion=label, R_H_low_mm=THRESHOLDS[j], R_H_high_mm=THRESHOLDS[j + 1],
                                            state_low=labels[j], state_high=labels[j + 1],
                                            value_low=values[j], value_high=values[j + 1], unit="EUR/kWp/year"))
        robust.append(row)
    orientation_robust = []
    for identifiers, q in orientations.groupby(PAIR, sort=False):
        q = q.sort_values("R_H_mm_event")
        if tuple(q.R_H_mm_event) != THRESHOLDS:
            raise ValueError("orientation comparison requires all five thresholds")
        ids = dict(zip(PAIR, identifiers))
        row = dict(**ids)
        for label in ("energy", "cost"):
            labels = q[f"{label}_preferred_orientation"].to_numpy()
            values = q[f"{label}_portrait_minus_landscape_pp"].to_numpy()
            row[f"{label}_ranking_stable"] = len(set(labels)) == 1
            row[f"{label}_rankings_5_10_15_20_25"] = ";".join(labels)
            row[f"{label}_all_landscape_better"] = bool((labels == "landscape").all())
            row[f"{label}_all_portrait_better"] = bool((labels == "portrait").all())
            for j in range(len(q) - 1):
                if labels[j] != labels[j + 1]:
                    transitions.append(dict(**ids, orientation="portrait vs landscape", conclusion=f"orientation_{label}",
                                            R_H_low_mm=THRESHOLDS[j], R_H_high_mm=THRESHOLDS[j + 1],
                                            state_low=labels[j], state_high=labels[j + 1], value_low=values[j], value_high=values[j + 1],
                                            unit="pp clean yield" if label == "energy" else "pp clean revenue"))
        orientation_robust.append(row)
    orientation_robust = pd.DataFrame(orientation_robust)
    robust = pd.DataFrame(robust).merge(orientation_robust, on=PAIR, validate="many_to_one")
    return robust, orientation_robust, pd.DataFrame(transitions)


def reference(frame):
    return frame[(frame.h_mm == REF_HEIGHT) & np.isclose(frame.f, REF_FRACTION)].copy()


def plot_summary(summary, path, market):
    q = reference(summary)
    q = q[q.market == market]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8))
        fig.subplots_adjust(left=.11, right=.99, bottom=.15, top=.80, hspace=.40, wspace=.37)
        for col, site in enumerate(energy.SITES):
            g = q[q.site == site]
            for module in ORDER:
                for orientation, ls, marker in (("portrait", "-", "o"), ("landscape", "--", "^")):
                    z = g[(g.module == module) & (g.orientation == orientation)].sort_values("R_H_mm_event")
                    for row, value in enumerate(("energy_none_band_extra_loss_pct", "margin10_npv25_eur_kwp")):
                        axes[row, col].plot(z.R_H_mm_event, z[value], color=MODULE_COLOR[module], ls=ls,
                                            marker=marker, ms=2.5, lw=.9, mew=.3)
            axes[0, col].set_title(site)
            for row in (0, 1):
                ax = axes[row, col]
                ax.axhline(0, color="0.25", lw=.5, zorder=1)
                ax.set_xticks(THRESHOLDS)
                ax.tick_params(pad=1.5)
                ax.text(.02, .98, "abcdef"[row * 3 + col], transform=ax.transAxes, va="top", fontweight="bold", fontsize=8)
            axes[0, col].set_ylim(bottom=0)
            axes[1, col].set_xlabel(r"Washing threshold $R_H$ [mm/event]")
        axes[0, 0].set_ylabel("Additional annual band loss\n[% of clean yield]")
        axes[1, 0].set_ylabel("Net value of +10 mm margin\nover 25 years [€ kWp$^{-1}$]")
        hd = [Line2D([], [], color=MODULE_COLOR[m], lw=1, label=PRETTY[m]) for m in ORDER]
        fig.legend(handles=hd, loc="upper center", bbox_to_anchor=(.55, .955), ncol=4)
        hd = [Line2D([], [], color="0.3", ls=ls, marker=mk, ms=3, label=o.capitalize())
              for o, ls, mk in (("portrait", "-", "o"), ("landscape", "--", "^"))]
        fig.legend(handles=hd, loc="upper center", bbox_to_anchor=(.55, .91), ncol=2)
        fig.suptitle(f"{market.capitalize()} market · rectangular band h = 60 mm, f = 0.1", fontsize=8, y=.99)
        fig.text(.55, .065, "Each line is one module/orientation; site panels have independent vertical scales.", ha="center", fontsize=6)
        fig.text(.55, .04, "Top: no manual cleaning. Bottom: independently optimised policies, minus added-glass cost; 5% discount, 25 years.", ha="center", fontsize=5.5)
        fig.text(.55, .015, "Complete washing assumed; thresholds are not calibrated for bands. Illustrative reference: R_H = 15 mm/event.", ha="center", fontsize=5.5)
        save(fig, path)


def load_and_audit(root):
    """Require identical stored assumptions and candidate policies before comparing thresholds."""
    tables, inputs, reference_settings, grid, invariant_y = [], [], None, None, None
    columns = ["washing_family", "washing_rule", "dry_hours", "transfer_event_rain_mm", "years", "climate_csvs",
               "hourly_time_convention", "daily_cleaning_boundaries", "simulation_source_sha256",
               "alpha", "deposition_g_m2_day", "spice_sweep_sha256", "event_action_time"]
    for threshold in THRESHOLDS:
        folder = root / "total_scenarios" / f"PH{threshold:g}mm" / "data"
        settings = json.loads((folder / "model_settings.json").read_text())
        if settings["washing_event_total_mm"] != threshold or settings["washing_intensity_mm_h"] is not None:
            raise ValueError("incorrect washing scenario metadata")
        if settings["washing_family"] != "total" or settings["washing_rule"] != "P_e >= P_H only":
            raise ValueError("expected event-total washing without an intensity override")
        assumptions = {c: settings[c] for c in columns}
        if reference_settings is None:
            reference_settings = assumptions
        elif assumptions != reference_settings:
            raise ValueError("other model assumptions differ between scenarios")
        path = folder / "model_runs.csv"
        runs = pd.read_csv(path)
        current_grid = runs.set_index(KEY + ["mode", "k"]).index.sort_values()
        current_y = runs.groupby(KEY)[["Y", "revenue"]].first().sort_index()
        if grid is None:
            grid, invariant_y = current_grid, current_y
        else:
            pd.testing.assert_index_equal(grid, current_grid)
            pd.testing.assert_frame_equal(invariant_y, current_y, rtol=1e-12, atol=1e-12)
        # Check the existing cost includes lost revenue and the actual manual cleaning charges.
        for market, (price, charge) in energy.MARKETS.items():
            z = runs[(runs.market == market) & (runs["mode"] != "robot")]
            np.testing.assert_allclose(z.cost, price * z.Y * z.loss + charge * z.cleans, atol=1e-10)
            np.testing.assert_allclose(z.cost_u, price * z.Y * z.loss_u + charge * z.cleans_u, atol=1e-10)
        tables.append(summarize_runs(runs, threshold))
        # Reconcile the margin calculation against the independently selected optima
        # saved by fig_energy.paper_figures (ground market).
        saved_best = pd.read_csv(folder / "best_strategy.csv")
        best_pivot = saved_best.pivot(index=["site", "module", "orientation", "h", "f"], columns="margin", values="bc_eur")
        result = tables[-1][tables[-1].market == "ground"].rename(columns={"h_mm": "h"}).set_index(best_pivot.index.names)
        expected = energy.AF * (best_pivot[0.] - best_pivot[10.])
        np.testing.assert_allclose(result.margin10_npv25_eur_kwp.reindex(expected.index) + result.margin10_capex_eur_kwp.reindex(expected.index), expected, atol=1e-10)
        inputs.append(dict(threshold_mm_event=threshold, file=str(path.relative_to(config.REPO_DIR)),
                           sha256=hashlib.sha256(path.read_bytes()).hexdigest(), rows=len(runs)))
    return pd.concat(tables, ignore_index=True), dict(identical_assumptions=True, identical_candidate_policy_grids=True,
                 identical_clean_production_and_prices=True, saved_costs_reconcile=True, margin_reconciles_with_existing_figures=True,
                 assumptions=reference_settings, inputs=inputs)
