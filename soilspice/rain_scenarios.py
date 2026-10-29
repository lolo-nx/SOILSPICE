"""Event-total washing: illustrative R_H=15 mm reference and 5/10/20/25 mm sensitivity."""
import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config, energy, rain
from .style import ORDER, ORIENTATIONS, SITE_COLOR, STYLE, save

ASSUMPTION = "Complete washing is assumed; thresholds are not calibrated for bottom-edge bands."
REF_HEIGHT = 60.0
REF_FRACTION = 0.1


def simulation_provenance():
    """Record the timestamp convention and calculation sources used by these exports."""
    return dict(hourly_time_convention="Open-Meteo hour ends shifted -1 h after PV calculation; interval starts in UTC",
                daily_cleaning_boundaries="UTC midnight",
                simulation_source_sha256={name: hashlib.sha256((config.REPO_DIR / "soilspice" / name).read_bytes()).hexdigest()
                                          for name in ("energy.py", "rain.py")})


def options(threshold):
    return {"washing_total": energy.washing_threshold(threshold)}


def scenario_tag(threshold):
    energy.washing_threshold(threshold)
    return "PH" + f"{threshold:g}".replace(".", "p") + "mm"


def write_scenario_settings(threshold, folder):
    """Metadata for a specifically requested scenario, including single-scenario CLI runs."""
    options(threshold)
    metadata = dict(washing_family="total",
                    scenario_role="illustrative reference" if threshold == energy.REFERENCE_TOTAL else "sensitivity",
                    washing_intensity_mm_h=None,
                    washing_event_total_mm=threshold,
                    washing_rule="P_e >= P_H only",
                    washing_assumption=ASSUMPTION, transfer_event_rain_mm=energy.LIGHT,
                    selected_reference_threshold=energy.REFERENCE_TOTAL, dry_hours=rain.DRY_HOURS, years=energy.YEARS,
                    intensity_definition="max(hourly precipitation [mm] / 1 h); not event mean",
                    event_action_time="one update at event closure, after 3 consecutive dry hours",
                    climate_csvs=energy.FILES, alpha=energy.ALPHA,
                    deposition_g_m2_day={site: energy.deposition(site) for site in energy.SITES},
                    spice_sweep_sha256=hashlib.sha256(config.EQUAL_HEIGHT_CSV.read_bytes()).hexdigest())
    metadata.update(simulation_provenance())
    target = folder / "data" / "model_settings.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")


def finish_economics_exports():
    """Verify the five retained event-total scenarios."""
    figures = ("fig_band_growth", "fig_energy_maps", "fig_design_levers", "fig_design_height",
               "fig_sensitivity", "fig_cleaning", "fig_band_shape", "fig_inverse_chain", "fig_pan_h_maps_SI")
    thresholds = energy.TOTAL_THRESHOLDS
    for threshold in thresholds:
        folder = config.CH4_DIR / "total_scenarios" / scenario_tag(threshold)
        for stem in figures:
            for extension in (".pdf", ".png"):
                if not (folder / (stem + extension)).is_file():
                    raise ValueError(f"Incomplete scenario export: {folder / (stem + extension)}")
        write_scenario_settings(threshold, folder)
    target = config.ECON_DIR / "rain_total_settings.json"
    metadata = json.loads(target.read_text())
    metadata["economics_scenarios_fully_generated"] = list(thresholds)
    metadata["economics_outputs"] = "total_scenarios/; illustrative reference PH15mm/"
    target.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")


def event_intervals(events, start, stop, threshold):
    """Real closure-to-closure gaps, plus explicitly censored study-boundary intervals."""
    selected = events[events.complete & (events.end > start) & (events.end <= stop)
                      & energy.washing_flags(events.rain, **options(threshold))].sort_values("end")
    times = pd.DatetimeIndex(selected.end)
    boundaries = pd.DatetimeIndex([start, *times, stop])
    intervals = pd.DataFrame({"start": boundaries[:-1], "end": boundaries[1:],
                              "days": np.diff(boundaries.asi8) / pd.Timedelta(days=1).value})
    intervals["left_censored"] = np.arange(len(intervals)) == 0
    intervals["right_censored"] = np.arange(len(intervals)) == len(intervals) - 1
    full = intervals[~intervals.left_censored & ~intervals.right_censored]
    longest = full.loc[full.days.idxmax()] if len(full) else None
    return dict(washing_events=len(times),
                max_between_washings_days=longest.days if longest is not None else np.nan,
                longest_gap_start=longest.start if longest is not None else pd.NaT,
                longest_gap_end=longest.end if longest is not None else pd.NaT,
                max_observed_without_washing_days=intervals.days.max(),
                initial_censored_days=intervals.days.iloc[0], final_censored_days=intervals.days.iloc[-1]), intervals


def reference_metrics(threshold):
    """Annual outputs for h=60 mm, f=0.1, rectangular bands, eight configurations, no manual cleaning.
    Transmittance is weighted by clean hourly production, as in energy.simulate.
    Losses use the same clean-production denominator; additional band loss subtracts the uniform-only case.
    """
    rain_options = options(threshold)
    threshold_column = "P_H_mm"
    combos = [(m, o, REF_HEIGHT, REF_FRACTION, 0.0) for m in ORDER for o in ORIENTATIONS]
    annual, details, interval_tables = [], [], []
    for site in energy.SITES:
        d, y, totals, kw = energy.simulation_inputs(site, combos, energy.deposition(site), **rain_options)
        lost, _, tau_mean, tracked = energy.simulate(
            y, totals, energy.deposition(site), energy.gains(combos), energy.lut(combos),
            track=np.arange(len(combos)), **kw)
        tau, lb, lu = np.asarray(tracked).transpose(1, 0, 2)
        losses = y * (1 - (1 - lb) * (1 - lu))
        extra = y * lb * (1 - lu)
        np.testing.assert_allclose(losses.sum(axis=0), lost, atol=1e-9)
        np.testing.assert_allclose((tau * y).sum(axis=0) / y.sum(axis=0), tau_mean, atol=1e-12)
        years = np.unique(d.index.year)
        start, stop = d.index[0], d.end_time.iloc[-1]
        event_info, intervals = event_intervals(energy.site_events(site), start, stop, threshold)
        interval_tables.append(intervals.assign(site=site, **{threshold_column: threshold}))
        washing = energy.washing_flags(d.rain, **rain_options)
        assert washing.sum() == event_info["washing_events"]
        for j, (module, orientation, h, f, margin) in enumerate(combos):
            common = dict(site=site, **{threshold_column: threshold}, module=module, orientation=orientation,
                          h_mm=h, f=f, margin_mm=margin, shape="rectangle", manual_cleaning="none")
            details.append(dict(**common, **event_info, washing_events_per_year=int(washing.sum()) / len(years),
                                tau_band_energy_weighted=tau_mean[j],
                                clean_energy_kwh_kwp_year=y[:, j].sum() / len(years),
                                total_energy_loss_pct=100 * losses[:, j].sum() / y[:, j].sum(),
                                additional_band_energy_loss_pct=100 * extra[:, j].sum() / y[:, j].sum(),
                                total_energy_lost_kwh_kwp_year=losses[:, j].sum() / len(years),
                                additional_band_energy_lost_kwh_kwp_year=extra[:, j].sum() / len(years)))
            for year in years:
                mask = d.index.year == year
                clean = y[mask, j].sum()
                closures = d.end_time.dt.year.to_numpy() == year
                annual.append(dict(**common, year=int(year), washing_events=int((washing & closures).sum()),
                                   tau_band_energy_weighted=np.sum(tau[mask, j] * y[mask, j]) / clean,
                                   clean_energy_kwh_kwp=clean,
                                   total_energy_loss_pct=100 * losses[mask, j].sum() / clean,
                                   additional_band_energy_loss_pct=100 * extra[mask, j].sum() / clean,
                                   total_energy_lost_kwh_kwp=losses[mask, j].sum(),
                                   additional_band_energy_lost_kwh_kwp=extra[mask, j].sum()))
    return pd.DataFrame(details), pd.DataFrame(annual), pd.concat(interval_tables, ignore_index=True)


def comparison_figure(summary, path):
    metrics = [("washing_events_per_year", "Assumed washing events/year"),
               ("max_observed_without_washing_days", "Longest observed period\nwithout washing [days]"),
               ("tau_band_energy_weighted", "Production-weighted band\n" + r"transmittance $\bar\tau_I$"),
               ("total_energy_loss_pct", "Annual total energy loss\n[% of clean yield]")]
    thresholds = energy.TOTAL_THRESHOLDS
    threshold_column = "P_H_mm"
    colors = {5: "#009E73", 10: "#0072B2", 15: "#9467BD", 20: "#E69F00", 25: "#D55E00"}
    n = len(thresholds)
    width = 0.8 / n
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
        for k, (ax, (metric, label)) in enumerate(zip(axes.flat, metrics)):
            for j, threshold in enumerate(thresholds):
                q = summary[summary[threshold_column] == threshold].set_index("site").loc[list(energy.SITES)]
                x = np.arange(len(energy.SITES)) + (j - (n - 1) / 2) * width
                ax.bar(x, q[metric], width * 0.95, color=colors[threshold], label=rf"$R_H={threshold:g}$ mm/event")
                if metric == "max_observed_without_washing_days":
                    ax.plot(x, q.max_between_washings_days, "o", ms=3, mfc="white", mec="black", mew=0.5)
                lo, hi = metric + "_min", metric + "_max"
                if lo in q:
                    ax.errorbar(x, q[metric], yerr=[q[metric] - q[lo], q[hi] - q[metric]],
                                fmt="none", color="black", lw=0.5, capsize=1.5)
            ax.set_xticks(range(len(energy.SITES)), energy.SITES)
            ax.set_ylabel(label)
            ax.text(0.015, 0.98, "abcd"[k], transform=ax.transAxes, va="top", fontweight="bold", fontsize=8)
            if metric == "tau_band_energy_weighted":
                ax.set_ylim(0, 1.08)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, ncol=min(3, n),
                   loc="upper center", bbox_to_anchor=(0.5, 0.99))
        fig.text(0.5, 0.075, "Rectangular band: h = 60 mm, f = 0.1; no manual cleaning. "
                 "Median and range over 4 modules × 2 orientations.", ha="center", fontsize=6)
        fig.text(0.5, 0.055, "Waiting periods include censored study boundaries; open circles: longest complete closure-to-closure interval.",
                 ha="center", fontsize=5.5)
        fig.text(0.5, 0.035, ASSUMPTION, ha="center", fontsize=6)
        fig.tight_layout(rect=(0, 0.11, 1, 0.91), h_pad=2.5, w_pad=1.8)
        save(fig, path)


def growth_comparison(series, path):
    thresholds = energy.TOTAL_THRESHOLDS
    threshold_column = "P_H_mm"
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(len(energy.SITES), len(thresholds), figsize=(7.2, 4.6), sharex=True, sharey=True)
        for r, site in enumerate(energy.SITES):
            for c, threshold in enumerate(thresholds):
                q = series[(series.site == site) & (series[threshold_column] == threshold)]
                ax = axes[r, c]
                ax.plot(pd.to_datetime(q.date), q.tau, color=SITE_COLOR[site], lw=0.5, drawstyle="steps-post")
                ax.set_ylim(0, 1.05)
                if r == 0:
                    ax.set_title(rf"$R_H={threshold:g}$" + "\nmm/event")
                if c == 0:
                    ax.set_ylabel(site + "\n" + r"Band $\tau_I$")
                ax.tick_params(pad=1.5)
                ax.set_xticks(pd.to_datetime(["2016-01-01", "2020-01-01", "2024-01-01"]), ["2016", "2020", "2024"])
        fig.suptitle("JA Solar, landscape · h = 60 mm, f = 0.1 · no manual cleaning", fontsize=8)
        fig.text(0.5, 0.02, ASSUMPTION, ha="center", fontsize=6)
        fig.tight_layout(rect=(0, 0.05, 1, 0.95), h_pad=0.8, w_pad=0.8)
        save(fig, path)


def run_comparison():
    """Run the five event-total thresholds with identical inputs."""
    options(1)
    thresholds = energy.TOTAL_THRESHOLDS
    threshold_column = "P_H_mm"
    symbol, unit = ("R_H", "mm/event")
    all_details, all_annual, all_intervals, all_series, all_maps, event_tables = [], [], [], [], [], []
    raw_events = []
    for site in energy.SITES:
        e = energy.site_events(site).reset_index().assign(site=site)
        e["P_e_mm"], e["I_max_e_mm_h"] = e.rain, e.rain_max
        raw_events.append(e)
    catalogue = pd.concat(raw_events, ignore_index=True)
    catalogue.to_csv(config.ECON_DIR / "rain_events.csv", index=False)
    for threshold in thresholds:
        print(f"Comparing {symbol}={threshold:g} {unit}", flush=True)
        folder = config.CH4_DIR / "total_scenarios" / scenario_tag(threshold)
        data = folder / "data"
        data.mkdir(parents=True, exist_ok=True)
        write_scenario_settings(threshold, folder)
        details, annual, intervals = reference_metrics(threshold)
        series = energy.growth_series(**options(threshold)).assign(**{threshold_column: threshold})
        maps = energy.tau_grid(**options(threshold)).assign(**{threshold_column: threshold})
        for name, frame in (("reference_metrics", details), ("annual_reference_metrics", annual),
                            ("washing_intervals", intervals), ("growth_series_reference", series), ("energy_maps", maps)):
            exported = frame.drop(columns=threshold_column) if name in ("growth_series_reference", "energy_maps") else frame
            exported.to_csv(data / f"{name}.csv", index=False)
        # The saved energy maps and tracked metrics must represent exactly the same model.
        check = maps[(maps.h == REF_HEIGHT) & np.isclose(maps.f, REF_FRACTION)].merge(details, on=["site", "module", "orientation"])
        np.testing.assert_allclose(100 * check.loss, check.total_energy_loss_pct, atol=1e-10)
        np.testing.assert_allclose(100 * (check.loss - check.loss_u), check.additional_band_energy_loss_pct, atol=1e-10)
        np.testing.assert_allclose(check.tau_mean, check.tau_band_energy_weighted, atol=1e-12)
        e = catalogue.copy()
        washing = energy.washing_flags(e.rain, **options(threshold))
        e["rain_regime"] = np.where(~e.complete, "incomplete", np.where(washing, "washing",
                                   np.where(e.rain >= energy.LIGHT, "transfer", "negligible")))
        event_tables.append(e.assign(**{threshold_column: threshold}))
        all_details.append(details); all_annual.append(annual); all_intervals.append(intervals)
        all_series.append(series); all_maps.append(maps)
    details, annual = pd.concat(all_details, ignore_index=True), pd.concat(all_annual, ignore_index=True)
    intervals, series = pd.concat(all_intervals, ignore_index=True), pd.concat(all_series, ignore_index=True)
    summary_rows = []
    vary = ["tau_band_energy_weighted", "total_energy_loss_pct", "additional_band_energy_loss_pct",
            "total_energy_lost_kwh_kwp_year", "additional_band_energy_lost_kwh_kwp_year"]
    for (site, threshold), q in details.groupby(["site", threshold_column], sort=False):
        row = dict(site=site, **{threshold_column: threshold}, h_mm=REF_HEIGHT, f=REF_FRACTION,
                   washing_events=int(q.washing_events.iloc[0]), washing_events_per_year=q.washing_events_per_year.iloc[0],
                   max_between_washings_days=q.max_between_washings_days.iloc[0],
                   longest_gap_start=q.longest_gap_start.iloc[0], longest_gap_end=q.longest_gap_end.iloc[0],
                   max_observed_without_washing_days=q.max_observed_without_washing_days.iloc[0])
        for column in vary:
            row.update({column: q[column].median(), column + "_min": q[column].min(), column + "_max": q[column].max()})
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    for name, frame in ((f"rain_total_summary", summary), (f"rain_total_reference_metrics", details),
                        (f"rain_total_annual_metrics", annual), (f"rain_total_washing_intervals", intervals),
                        (f"rain_total_growth_series", series), (f"rain_total_energy_maps", pd.concat(all_maps, ignore_index=True)),
                        (f"rain_events_total_scenarios", pd.concat(event_tables, ignore_index=True))):
        frame.to_csv(config.ECON_DIR / f"{name}.csv", index=False)
    settings = dict(washing_family="total", scenario_role="reference plus sensitivity",
                    washing_rule="P_e >= P_H only",
                    washing_assumption=ASSUMPTION, thresholds=list(thresholds), threshold_units=unit,
                    selected_reference_threshold=energy.REFERENCE_TOTAL, transfer_event_rain_mm=energy.LIGHT,
                    event_definition="3 consecutive zero-rain hours separate events", dry_hours=rain.DRY_HOURS,
                    intensity_definition="max(hourly precipitation in mm / 1 hour); not event mean intensity",
                    event_action_priority=["assumed full washing from event total", "otherwise transfer if P_e >= 1 mm", "otherwise no effect"],
                    event_action_time="once, at closure after 3 dry hours; incomplete events have no action",
                    years=energy.YEARS, timezone="UTC", h_mm=REF_HEIGHT, f=REF_FRACTION,
                    manual_cleaning="none for scenario comparison", tau_weight="clean hourly production",
                    loss_denominator="clean hourly production; positive loss; annual percentages are yield weighted",
                    summary_aggregation="median and range of 4 modules x 2 orientations, detailed values exported",
                    gap_definition="closure to next closure, across years; incomplete boundary intervals separately censored",
                    annualization="2016–2025; event counts divided by 10; energy loss weighted over available hourly yields",
                    source_coverage_note="The source files lack the hour ending 2026-01-01 00:00 UTC; no synthetic weather added.",
                    climate_csvs=energy.FILES, source_sha256={site: hashlib.sha256((config.CLIMATE_DIR / f).read_bytes()).hexdigest()
                                                            for site, f in energy.FILES.items()},
                    deposition_g_m2_day={site: energy.deposition(site) for site in energy.SITES}, alpha=energy.ALPHA,
                    spice_sweep=str(config.EQUAL_HEIGHT_CSV.relative_to(config.REPO_DIR)),
                    spice_sweep_sha256=hashlib.sha256(config.EQUAL_HEIGHT_CSV.read_bytes()).hexdigest())
    settings.update(simulation_provenance())
    (config.ECON_DIR / f"rain_total_settings.json").write_text(json.dumps(settings, indent=2, ensure_ascii=False) + "\n")
    comparison_figure(summary, config.CH4_DIR / "fig_rain_total_comparison")
    growth_comparison(series, config.CH4_DIR / "fig_band_growth_total_comparison")
    print(summary[["site", threshold_column, "washing_events_per_year", "max_between_washings_days",
                   "tau_band_energy_weighted", "total_energy_loss_pct", "additional_band_energy_loss_pct"]].round(4).to_string(index=False), flush=True)
    return summary
