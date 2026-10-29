"""LaTeX tables (booktabs) written from DataFrames."""
from pathlib import Path

import pandas as pd


def to_latex(df, path, caption, label, decimals=2, decimal_comma=False, column_names=None):
    """Write `df` as a booktabs table. Numbers are rounded to `decimals`;
    missing values become '---'."""
    def fmt(v):
        if isinstance(v, float):
            if pd.isna(v):
                return "---"
            s = f"{v:.{decimals}f}"
            return s.replace(".", ",") if decimal_comma else s
        return str(v).replace("_", r"\_").replace("%", r"\%")

    names = column_names or [str(c).replace("_", r"\_").replace("%", r"\%") for c in df.columns]
    lines = [r"\begin{table}[htbp]", r"\centering", rf"\caption{{{caption}}}", rf"\label{{{label}}}",
             r"\begin{tabular}{@{}" + "l" * len(df.columns) + r"@{}}", r"\toprule",
             " & ".join(names) + r" \\", r"\midrule"]
    lines += [" & ".join(fmt(v) for v in row) + r" \\" for row in df.itertuples(index=False)]
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  saved {path}")


SITE_NOTES = {"Neuchâtel": "a", "Almería": "b", "Doha": "c"}   # footnotes of the soiling rates (sources in the caption)


def sites_table(sites, bound, path, *, washing_total=None):
    """Table sites of the paper (LaTeX, same layout as the manuscript): ERA5 grid point, yield, precipitation, rainy
    days, soiling rate, deposition D and L_min (event total >= R_H; complete washing assumed).
    sites: energy.site_table(); bound: energy.cleaning_bound()."""
    from . import energy as E
    threshold = E.washing_threshold(E.REFERENCE_TOTAL if washing_total is None else washing_total)
    wash_column = "washing_event_total_mm"
    condition = rf"$P_e\geq {threshold:g}$~mm/event"
    s = sites.set_index("site")
    lmin = bound[bound[wash_column] == threshold].pivot(index="site", columns="market", values="L_min")
    names = list(s.index)

    def deg(v, pos, neg):
        return f"{abs(v):.2f}$^\\circ${pos if v >= 0 else neg}"

    rows = [("Latitude", [deg(s.loc[n, "Latitude"], "N", "S") for n in names]),
            ("Longitude", [deg(s.loc[n, "Longitude"], "E", "W") for n in names]),
            ("Yield [kWh~kWp$^{-1}$~yr$^{-1}$]", [f"{s.loc[n, 'Yield [kWh/kWp/yr]']:.0f}" for n in names]),
            ("Precipitation [mm~yr$^{-1}$]", [f"{s.loc[n, 'Precipitation [mm/yr]']:.0f}" for n in names]),
            ("Days $\\geq$ 1~mm [yr$^{-1}$]", [f"{s.loc[n, 'Days >= 1 mm [/yr]']:.0f}" for n in names]),
            ("Washing events [yr$^{-1}$]",
             [f"{s.loc[n, 'Washing events [/yr]']:.1f}" for n in names]),
            ("Soiling rate [\\%~day$^{-1}$]", [f"{s.loc[n, 'Soiling rate [%/day]']:.2f}\\textsuperscript{{{SITE_NOTES.get(n, '')}}}"
                                              for n in names]),
            ("$D$ [g~m$^{-2}$~day$^{-1}$]", [f"{v:.2f}" if (v := s.loc[n, 'D [g/m2/day]']) >= 0.1 else f"{v:.3f}" for n in names]),
            None,
            ("$L_{\\mathrm{min}}$, ground [\\%]", [f"{lmin.loc[n, 'ground']:.1f}" for n in names]),
            ("$L_{\\mathrm{min}}$, rooftop [\\%]", [f"{lmin.loc[n, 'rooftop']:.1f}" for n in names])]
    tex = {"Neuchâtel": "Neuch\\^atel", "Almería": "Almer\\'ia"}
    out = ["\\begin{table}[!b]",
           "\\caption{Sites (ERA5 grid points, 2016--2025). Soiling rates are those measured on uniformly soiled modules "
           "during dry periods at the site or the nearest site with published data; $D$ follows from the rate and the "
           "mass--transmittance relation of \\cite{boyleNaturalSoilingPhotovoltaic2015}. $L_{\\mathrm{min}}$ is the "
           "smallest band loss for which one cleaning can pay (Table~\\ref{tab:clean}). "
           f"Complete washing is assumed when {condition}; "
           "this threshold is not calibrated for bands. Events end after 3 consecutive dry hours."
           "\\label{tab:sites}}",
           "\\definecolor{tabgrey}{gray}{0.92}",
           "\\begin{tabularx}{\\columnwidth}{X" + ">{\\centering\\arraybackslash}p{1.45cm}" * len(names) + "}",
           "\\toprule",
           " & " + " & ".join(f"\\textbf{{{tex.get(n, n)}}}" for n in names) + " \\\\",
           "\\midrule"]
    k = 0
    for r in rows:
        if r is None:
            out.append("\\midrule"); continue
        out.append(("\\rowcolor{tabgrey} " if k % 2 == 0 else "") + r[0] + " & " + " & ".join(r[1]) + " \\\\")
        k += 1
    out += ["\\bottomrule", "\\end{tabularx}",
            "\\par\\smallskip{\\footnotesize \\textsuperscript{a} Karlsruhe, 16$^\\circ$ tilt \\cite{wagnerVogel2019}. "
            "\\textsuperscript{b} Granada, 30$^\\circ$, summer dry season \\cite{micheliEconomicsSeasonalPhotovoltaic2021}. "
            "\\textsuperscript{c} Doha, 22$^\\circ$, multi-year mean \\cite{figgis2018}.}",
            "\\end{table}"]
    Path(path).write_text("\n".join(out) + "\n")


def band_params_table(path, cite_karlsruhe="WorldEnergyMeteorology", cite_doha="ahziInvestigationPVSoiling", *,
                      washing_total=None):
    """Table of the parameters of the band growth model (LaTeX, short version: no 'read in' column).
    Values are read from soilspice.energy so that the table follows the code."""
    from . import energy as E
    threshold = E.washing_threshold(E.REFERENCE_TOTAL if washing_total is None else washing_total)
    unit = "mm/event"
    wash_name = "Washing event total $R_H$"
    total_range = ", ".join(f"{value:g}" for value in E.TOTAL_THRESHOLDS)
    pg, cg = E.MARKETS["ground"]; pr, cr = E.MARKETS["rooftop"]
    rows = [
        ("Soiling rate $k$, Neuch\\^atel", f"{E.RATE['Neuchâtel']:.2f}\\%~day$^{{-1}}$", "$\\times 0.5$, $\\times 2$",
         f"Karlsruhe, dry periods \\cite{{{cite_karlsruhe}}}"),
        ("Soiling rate $k$, Almer\\'ia", f"{E.RATE['Almería']:.2f}\\%~day$^{{-1}}$", "$\\times 0.5$, $\\times 2$",
         "Granada, summer dry season \\cite{micheliEconomicsSeasonalPhotovoltaic2021}"),
        ("Soiling rate $k$, Doha", f"{E.RATE['Doha']:.2f}\\%~day$^{{-1}}$", "$\\times 0.5$, $\\times 2$",
         f"Doha test facility, multi-year mean \\cite{{{cite_doha}}}"),
        ("Loss per unit mass", "4.1\\% per g~m$^{-2}$", "--", "Linear up to 2~g~m$^{-2}$ \\cite{boyleNaturalSoilingPhotovoltaic2015}"),
        ("Attenuation $\\alpha$", f"{E.ALPHA:g}~m$^{{2}}$~g$^{{-1}}$", "$\\times 0.75$, $\\times 1.25$", "Beer--Lambert, initial slope above"),
        ("Transfer event rain $R_L$", f"{E.LIGHT:g}~mm", "0.5, 2", "Scenario assumption; cumulative rain per event"),
        (wash_name, f"{threshold:g}~{unit}", total_range,
         "Complete washing assumed; thresholds not calibrated for bands"),
        ("Event separation", "3 dry hours", "fixed", "An event ends after 3 consecutive hours with zero rain"),
        ("Trapped fraction $f$", "0.1", "0.03, 0.3", "Not measured (maps 0.01--1)"),
        ("Band height $h$", "10--100~mm", "--", "5--10~cm observed \\cite{yanYao2018}"),
        ("Band shape", "rectangle, triangle", "--", "Along the frame \\cite{cuiCharacterizationSoilingBands2021a}"),
        ("Electricity price $p$", f"{pg:g} ({pr:g})~\\euro~kWh$^{{-1}}$", "--", "Ground (rooftop) \\cite{ilseTechnoEconomicAssessmentSoiling2019a}"),
        ("Cleaning cost $C$", f"{cg:g} ({cr:g})~\\euro~kWp$^{{-1}}$", "$\\times 0.5$, $\\times 2$",
         "Ground (rooftop) \\cite{ilseTechnoEconomicAssessmentSoiling2019a}"),
        ("Cleaning robot", "2.4--8.2~\\euro~m$^{-2}$", "range",
         "Investment over 10 years \\cite{ilseTechnoEconomicAssessmentSoiling2019a}"),
        ("Discount rate, lifetime", "5\\%, 25 years", "--", "\\cite{ilseTechnoEconomicAssessmentSoiling2019a}"),
        ("Module price", "0.10~\\euro~Wp$^{-1}$", "--", "Assumption (margin cost)"),
    ]
    out = ["\\begin{table*}[!t]",
           "\\caption{Parameters of the band growth model. Sensitivity: low and high values of the one-at-a-time "
           "analysis (Extended Data Fig.~\\ref{fig:ed-sens}).\\label{tab:bandparams}}",
           "\\definecolor{tabgrey}{gray}{0.92}",
           "\\begin{tabularx}{\\textwidth}{p{3.6cm}p{3.0cm}p{2.3cm}X}",
           "\\toprule",
           "\\textbf{Parameter} & \\textbf{Value} & \\textbf{Sensitivity} & \\textbf{Source} \\\\",
           "\\midrule"]
    for k, r in enumerate(rows):
        out.append(("\\rowcolor{tabgrey} " if k % 2 == 0 else "") + " & ".join(r) + " \\\\")
    out += ["\\bottomrule", "\\end{tabularx}", "\\end{table*}"]
    Path(path).write_text("\n".join(out) + "\n")
