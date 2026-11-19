"""Refresh PH15mm robot economics from saved losses; no physical resimulation.

Main Fig. 8g: explicit linear-remaining-life credit. Export a named no-residual
sensitivity and both sets of full cash flows, using the same figure style.
Run: python chapters/update_robot_economics.py
Then: python Nature-builder-economics.py
"""
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from soilspice import config, fig_energy, tables
from soilspice.robot_economics import (DISCOUNT_RATE, EVALUATION_YEARS,
    ROBOT_INVESTMENT_EUR_KWP, ROBOT_LIFETIME_YEARS, ROBOT_MAINTENANCE_FRACTION,
    RESIDUAL_MODES, annuity_factor, robot_cashflows, robot_ownership_cost, robot_strategy_costs)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assert_physics_unchanged():
    current = (config.REPO_DIR/'soilspice/energy.py').read_text(encoding='utf-8')
    previous = subprocess.check_output(['git','show','HEAD:soilspice/energy.py'], cwd=config.REPO_DIR).decode('utf-8')
    old = {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(previous).body
           if isinstance(node, ast.FunctionDef)}
    new = {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(current).body
           if isinstance(node, ast.FunctionDef)}
    # model_runs only receives a bookkeeping docstring change; simulation,
    # weather/yield, LUT, deposition and cleaning decisions must be identical.
    changed = [name for name in old if name != 'model_runs' and old[name] != new[name]]
    assert not changed, changed
    a = next(n for n in ast.parse(previous).body if isinstance(n,ast.FunctionDef) and n.name=='model_runs')
    b = next(n for n in ast.parse(current).body if isinstance(n,ast.FunctionDef) and n.name=='model_runs')
    assert ast.dump(ast.Module(body=a.body[1:],type_ignores=[])) == ast.dump(ast.Module(body=b.body[1:],type_ignores=[]))


def main():
    scenario = config.CH4_DIR/'total_scenarios/PH15mm'
    data = scenario/'data'
    out = data/'robot_economics'
    out.mkdir(parents=True,exist_ok=True)
    protected = {p:sha(p) for p in data.rglob('*') if p.is_file() and out not in p.parents}
    for p in scenario.glob('*.pdf'):
        if p.stem not in ('fig_cleaning','fig_cleaning_with_residual','fig_cleaning_without_residual',
                           'fig_design_cleaning_panel','fig_design_cleaning_panel_with_residual',
                           'fig_design_cleaning_panel_without_residual'):
            protected[p]=sha(p)
            png=p.with_suffix('.png')
            if png.exists(): protected[png]=sha(png)
    # No changes to electrical data, netlists, or climate inputs are permitted.
    for directory in (config.NETLIST_DIR,config.CLIMATE_DIR,config.EQUAL_HEIGHT_DIR):
        protected.update({p:sha(p) for p in directory.rglob('*') if p.is_file()})
    # This economic parameter table is intentionally refreshed below.
    assert_physics_unchanged()
    runs = pd.read_csv(data/'model_runs.csv')
    best = pd.read_csv(data/'best_strategy.csv')
    original_numbers = pd.read_csv(scenario/'numbers_in_text.csv')
    reference = runs[(runs.market=='ground') & (runs['mode']=='robot')
                     & np.isclose(runs.h,60) & np.isclose(runs.f,.1) & np.isclose(runs.margin,0)]
    assert len(reference)==24 and reference.groupby('site').size().eq(8).all()
    summaries, ownership, cashflows, config_costs, plot_summaries = [], [], [], [], {}
    for residual in RESIDUAL_MODES:
        costs = robot_strategy_costs(reference,residual_mode=residual)
        config_costs.append(costs)
        for bound,investment in zip(('low','mid','high'),
                (ROBOT_INVESTMENT_EUR_KWP[0],np.mean(ROBOT_INVESTMENT_EUR_KWP),ROBOT_INVESTMENT_EUR_KWP[1])):
            ownership.append(dict(bound=bound,**robot_ownership_cost(investment,residual_mode=residual)))
            cashflows.append(robot_cashflows(investment,residual_mode=residual).assign(
                residual_mode=residual,investment_bound=bound))
        for site,q in costs.groupby('site',sort=False):
            row=dict(site=site,residual_mode=residual,module_orientation_count=len(q),h_mm=60,f=.1)
            for case in ('band','noband'):
                for bound in ('low','mid','high'):
                    for suffix in ('pct_clean_revenue','total_eur_kwp_year','total_pv_eur_kwp'):
                        column=(f'{case}_{bound}_pct_clean_revenue' if suffix=='pct_clean_revenue'
                                else f'{case}_total_{bound}_{suffix.removeprefix("total_")}')
                        row[f'{case}_{bound}_{suffix}']=float(q[column].median())
            # Reconstruct the actual previous model, including the original
            # midpoint and 10-year undiscounted annualisation, for comparison.
            old_annual=np.mean(ROBOT_INVESTMENT_EUR_KWP)/10
            row['previous_band_mid_pct_clean_revenue']=float((100*(q.cost+old_annual)/q.revenue).median())
            row['band_effect_mid_pp']=float((100*(q.cost-q.cost_u)/q.revenue).median())
            summaries.append(row)
        name='fig_cleaning_without_residual' if residual=='none' else 'fig_cleaning'
        plot_summaries[residual]=fig_energy.cleaning_figures(runs,best,scenario,
                                            robot_residual_mode=residual,output_name=name)
        plot_summaries[residual].to_csv(out/f'cleaning_numbers_{residual}.csv',index=False)

    summary=pd.DataFrame(summaries)
    summary.to_csv(out/'robot_cost_summary.csv',index=False)
    pd.DataFrame(ownership).to_csv(out/'ownership_costs.csv',index=False)
    pd.concat(cashflows,ignore_index=True).to_csv(out/'discounted_cashflows.csv',index=False)
    pd.concat(config_costs,ignore_index=True).to_csv(out/'per_configuration_costs.csv',index=False)
    updated=plot_summaries['linear_remaining_life']
    # All three manual/no-cleaning strategy bars remain numerically identical.
    manual=updated[(updated.item=='strategy') & (updated.strategy!='Robot')]
    for _,row in manual.iterrows():
        old=original_numbers[(original_numbers.item=='strategy') & (original_numbers.site==row.site)
                             & (original_numbers.strategy==row.strategy)].iloc[0]
        np.testing.assert_allclose([row.cost_band,row.cost_noband],[old.cost_band,old.cost_noband],rtol=1e-12,atol=1e-12)
    # Replace only the three Robot rows in the paper's saved numbers.
    numbers=original_numbers.copy()
    robots=updated[(updated.item=='strategy') & (updated.strategy=='Robot')]
    for _,row in robots.iterrows():
        mask=(numbers.item=='strategy') & (numbers.strategy=='Robot') & (numbers.site==row.site)
        assert mask.sum()==1
        for column in ('cost_band','cost_noband','robot_residual_mode','cost_band_low','cost_band_high',
                       'cost_noband_low','cost_noband_high'):
            if column not in numbers:
                numbers[column]='' if column=='robot_residual_mode' else np.nan
            numbers.loc[mask,column]=row[column]
    nonrobot=~((numbers.item=='strategy') & (numbers.strategy=='Robot'))
    pd.testing.assert_frame_equal(numbers.loc[nonrobot,original_numbers.columns],
                                  original_numbers.loc[nonrobot],check_exact=True)
    numbers.to_csv(scenario/'numbers_in_text.csv',index=False)
    tables.band_params_table(scenario/'table_band_parameters.tex',washing_total=15.)

    settings=dict(horizon_years=EVALUATION_YEARS,discount_rate=DISCOUNT_RATE,
        robot_lifetime_years=ROBOT_LIFETIME_YEARS,purchase_years=[0,12,24],
        initial_investment_eur_kwp=ROBOT_INVESTMENT_EUR_KWP,
        maintenance_fraction_of_initial_investment_per_year=ROBOT_MAINTENANCE_FRACTION,
        maintenance_payment_years=list(range(1,26)),annuity_factor=annuity_factor(),
        main_figure_residual_mode='linear_remaining_life',alternative_residual_mode='none',
        alternative_terminal_fraction=11/12,terminal_credit_year=25,
        discounting_convention='Purchases at t=0,12,24; annual recurring costs/revenue at t=1..25',
        annual_energy_cost_convention='Saved 10-year mean annual losses held stationary over 25-year economic horizon',
        inflation_escalation_degradation='Unchanged: no added escalation or degradation',
        robot_energy_water_or_extra_operating_costs='No new operating assumptions: only specified 5% maintenance added',
        reference=dict(market='ground',h_mm=60,f=.1,margin_mm=0),
        investment_uncertainty='Input range endpoints, not statistical confidence intervals',
        no_band_markers='Same ownership EAC as band case; independently calculated uniform-layer energy loss',
        model_runs_cost_contract='Robot cost/cost_u are energy-loss costs only; no ownership already included')
    (out/'settings.json').write_text(json.dumps(settings,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    write_report(summary,updated,out,settings)
    assert all(sha(p)==value for p,value in protected.items()),'A protected physical input/output changed'
    audit=dict(physical_functions_unchanged=True,saved_model_runs_unchanged=True,
        protected_files_sha256={str(p.relative_to(config.REPO_DIR)):value for p,value in protected.items()},
        manual_strategy_bar_values_unchanged=True,nonrobot_numbers_unchanged=True,
        main_case_residual_mode='linear_remaining_life',
        ownership_pv_equals_eac_times_annuity=True,
        notes='Calculations checked by independent tests in tests/test_robot_economics.py')
    (out/'verification.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(summary[['site','residual_mode','band_low_pct_clean_revenue','band_mid_pct_clean_revenue',
                   'band_high_pct_clean_revenue','noband_mid_pct_clean_revenue']].to_string(index=False))


def write_report(summary,plots,out,settings):
    order=('Almería','Doha','Neuchâtel')
    def case(mode,site):
        return summary[(summary.residual_mode==mode)&(summary.site==site)].iloc[0]
    def vals(row):
        return f'{row.band_mid_pct_clean_revenue:.2f}% ({row.band_low_pct_clean_revenue:.2f}–{row.band_high_pct_clean_revenue:.2f}%)'
    results=('For daily robotic cleaning, the median total annual soiling costs are '+
        ', '.join(f'{vals(case("linear_remaining_life",site))} in {site}' for site in order)+
        ' of clean annual revenue (60 mm rectangular band, f=0.1, no additional glass margin, '
        'ground-market prices, PH15 mm scenario). Costs include remaining production losses, '
        'robot purchases/replacements and maintenance. Investment ranges are propagated at '
        'configuration level before taking the median over four modules and two orientations. '
        'The main case includes an explicit linear remaining-life credit in year 25. '
        'Without a terminal residual credit, the corresponding totals are '+
        ', '.join(f'{vals(case("none",site))} in {site}' for site in order)+'.')
    methods=('All cleaning strategies are compared over a 25-year economic horizon at a 5% discount rate. '
        'The annual costs and clean revenues are held constant at their simulated ten-year mean values, '
        'as in the existing model. Robot initial investment is 2.4–8.2 EUR/m², converted using the '
        'existing 4.6 m²/kWp factor to 11.04–37.72 EUR/kWp. Robots are purchased at years 0, 12 and 24, '
        'with unchanged real purchase cost at each replacement. Maintenance and repair cost 5% of '
        'initial investment annually, with year-end payments in years 1–25. The present value '
        'of these costs is divided by the 25-year annuity factor (14.09394456604475) to obtain '
        'an equivalent annual ownership cost. This is added once to the simulated residual '
        'production-loss cost in both with-band and no-band cases. Manual cleaning fees are '
        'not applied to the daily robot. The main case uses a linear remaining-life '
        'credit of 11/12 of the year-24 robot purchase '
        'cost at year 25, discounted to year zero. This is an assumed service-value credit, '
        'not an observed or guaranteed resale value. Zero residual credit is evaluated '
        'separately as a sensitivity. No residual credit is assigned to the '
        'robots exhausted at replacement. No new energy, water, escalation or degradation '
        'assumptions are introduced.')
    caption=('g, Total soiling costs expressed as equivalent annual costs divided by clean annual '
        'revenue, for a 60 mm rectangular band, f=0.1 and zero added glass margin in the PH15 mm '
        'scenario. Bars show medians over four module types and two orientations; robot bar '
        'centres use the midpoint investment (24.38 EUR/kWp). Robot error bars span the '
        '11.04–37.72 EUR/kWp initial-investment range, propagated through the 25-year '
        'discounted calculation including purchases at years 0, 12 and 24 and annual '
        'maintenance/repair of 5% of initial investment. A 5% discount rate and a terminal '
        'remaining-life credit of 11/12 of the year-24 purchase cost at year 25 are used '
        'in the main figure. Black horizontal markers show the '
        'corresponding no-band costs, including the same robot ownership cost. The investment '
        'range is not a statistical confidence interval. Other cleaning-strategy assumptions '
        'and panels a–f are unchanged; a figure without terminal residual credit '
        'is provided separately.')
    (out/'manuscript_text.md').write_text('## Results\n\n'+results+'\n\n## Methods\n\n'+methods+
        '\n\n## Figure 8 caption: panel g\n\n'+caption+'\n',encoding='utf-8')
    for name,text in (('results',results),('methods',methods),('figure8g_caption',caption)):
        tex=text.replace('%',r'\%').replace('–','--').replace('m²',r'm$^2$')
        (out/f'{name}.tex').write_text(tex+'\n',encoding='utf-8')
    lines=['# PH15 mm: robotic cleaning economic revision','','## What changed','',
        'Previous ownership charge: investment / 10, without discounting, replacement timing '
        'or maintenance. New charge: discounted 25-year ownership cost converted to an equivalent '
        'annual cost, plus unchanged saved annual lost revenue.','',
        'Main figure: **linear remaining-life residual credit**, as selected by the user. '
        'The named alternative assumes zero residual value.','',
        'For initial investment C, r=0.05 and A=sum((1+r)^(-t), t=1..25):','',
        '    PV_no_residual = C [1 + (1+r)^(-12) + (1+r)^(-24)] + 0.05 C A',
        '    PV_with_residual = PV_no_residual - (11/12) C (1+r)^(-25)',
        '    EAC = PV_with_residual / A  (main case)',
        '    annual total cost (% clean revenue) = 100 (saved energy-loss cost + EAC) / clean annual revenue','',
        'For manual strategies, multiplying annual costs and revenue by the same annuity factor '
        'leaves the plotted cost/revenue ratio unchanged. Comparing their existing annual '
        'costs to robot EAC is therefore consistent. The 2016–2025 meteorological record '
        'still supplies mean annual quantities; it is not a separate 10-year economic horizon.','',
        '## Exact plotted robot values','',
        '| Site | Residual assumption | Band low (%) | Band midpoint (%) | Band high (%) | No-band midpoint (%) |',
        '|---|---|---:|---:|---:|---:|']
    for _,r in summary.iterrows():
        lines.append(f'| {r.site} | {r.residual_mode} | {r.band_low_pct_clean_revenue:.9f} | '
            f'{r.band_mid_pct_clean_revenue:.9f} | {r.band_high_pct_clean_revenue:.9f} | {r.noband_mid_pct_clean_revenue:.9f} |')
    lines+=['','## Impact on other panels and conclusions','',
        'Only robot costs, their investment-range error bars and no-band robot markers in '
        'panel g change. Panels a–f do not use robot ownership costs: the manual optima '
        'explicitly exclude robots. The cleaning frequencies, physical energy losses, '
        'transmittances, band growth, threshold sensitivities and electrical calculations '
        'are unchanged. Other rainfall scenarios are not re-exported; a future full run '
        'uses the corrected robot ownership accounting.','',
        'A common ownership charge cancels exactly in the paired band-minus-no-band robot '
        'difference for each configuration. Site-level differences of medians need not '
        'equal the median of paired differences. All uncertainties here are investment '
        'bounds, separate from module/orientation variability.','']
    for site in order:
        r=case('linear_remaining_life',site)
        q=plots[(plots.item=='strategy') & (plots.site==site)]
        manual=q[q.strategy.isin(['Best fixed interval','Trigger rule'])]
        best=manual.loc[manual.cost_band.idxmin()]
        rain=q[q.strategy=='Rain only'].iloc[0]
        lines.append(f'- {site}: old robot total {r.previous_band_mid_pct_clean_revenue:.9f}%; '
            f'new robot total {r.band_mid_pct_clean_revenue:.9f}%; '
            f'best displayed manual strategy {best.strategy}: {best.cost_band:.9f}%; '
            f'rain only {rain.cost_band:.9f}%. '
            f'Robot cheaper than the best displayed manual strategy at midpoint: '
            f'{bool(r.band_mid_pct_clean_revenue<best.cost_band)}. '
            f'At low and high robot investment respectively: '
            f'{bool(r.band_low_pct_clean_revenue<best.cost_band)}, '
            f'{bool(r.band_high_pct_clean_revenue<best.cost_band)}. '
            f'Robot band-only effect (median paired difference) {r.band_effect_mid_pp:.9f} percentage points.')
    lines+=['','The robot almost removes the band-specific loss, but its total economic '
        'cost includes a substantial ownership charge. Avoid interpreting robot total '
        'soiling cost as a band-only loss or claiming one strategy is optimal across '
        'the whole investment range from midpoint bars alone.','',
        'See manuscript_text.md and the .tex snippets for ready-to-use Results, Methods '
        'and caption text; CSVs retain full precision. All physical inputs and saved '
        'model runs are protected by hashes; numerical cash-flow tests are independent '
        'of the figure implementation.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
