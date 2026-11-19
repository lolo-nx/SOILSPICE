"""Compact SI tables derived from the irradiance sensitivity CSVs.

Writes booktabs LaTeX snippets and vector PDF/PNG previews. No simulations.
"""
from pathlib import Path
import hashlib
import json
import sys
import textwrap

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from soilspice import config
from soilspice.style import ORDER, STYLE, module_label, save


def latex_table(path,title,label,headers,rows,note):
    columns='l'+'r'*(len(headers)-1)
    def escape(value):
        return str(value).replace('%',r'\%')
    lines=[r'% Requires \usepackage{booktabs}. Table numbering is set by the SI document.',
           r'\begin{table}[t]',r'\centering',r'\caption{\textbf{'+title+r'}}',
           r'\label{'+label+r'}',r'\small',r'\setlength{\tabcolsep}{6pt}',
           r'\begin{tabular}{'+columns+r'}',r'\toprule']
    lines += [' & '.join(r'\shortstack{'+escape(h).replace('\n',r'\\')+'}' for h in headers)+r' \\',r'\midrule']
    for row in rows:
        lines.append(' & '.join(escape(v) for v in row)+r' \\')
    lines += [r'\bottomrule',r'\end{tabular}',r'\par\smallskip',
              r'\begin{minipage}{\linewidth}',r'\footnotesize '+escape(note),
              r'\end{minipage}',r'\end{table}']
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')


def preview(path,title,headers,rows,widths,note):
    lines=textwrap.wrap(note,width=141)
    end=19+6*len(rows)
    height=end+7+3*len(lines)+3
    edges=[3]
    for width in widths:
        edges.append(edges[-1]+177*width)
    with plt.rc_context(STYLE):
        fig=plt.figure(figsize=(183/25.4,height/25.4))
        ax=fig.add_axes([0,0,1,1])
        ax.set(xlim=(0,183),ylim=(height,0))
        ax.axis('off')
        ax.text(3,4.7,title,fontsize=9,fontweight='bold',va='center')
        for y,lw in ((10,.8),(19,.45),(end,.8)):
            ax.plot([3,180],[y,y],color='black',lw=lw)
        for col,header in enumerate(headers):
            x=edges[col] if col==0 else (edges[col]+edges[col+1])/2
            ax.text(x,14.5,header,ha='left' if col==0 else 'center',va='center',fontsize=7.4)
        for index,row in enumerate(rows):
            for col,value in enumerate(row):
                ax.text(edges[col] if col==0 else edges[col+1]-2,22+6*index,value,
                        ha='left' if col==0 else 'right',va='center',fontsize=8)
        for index,line in enumerate(lines):
            ax.text(3,end+6+3*index,line,fontsize=6.8,va='center')
        save(fig,path,dpi=300)


def main():
    source=config.OUTPUT_DIR/'stc_sensitivity'
    out=source/'paper_tables'
    out.mkdir(exist_ok=True)
    power_path=source/'simulation_results.csv'
    annual_path=source/'annual_irradiance/annual_results_per_configuration.csv'
    power=pd.read_csv(power_path)
    power=power[(power.status=='ok')&(power.irradiance_W_m2<1000)]
    rows=[]
    power_data=[]
    for name in ORDER:
        q=power[power.module==name]
        n=int((q.Delta_L_pp.abs()<=2).sum())
        maximum=float(q.Delta_L_pp.abs().max())
        rows.append([module_label(name),f'{n} / {len(q)}',f'{maximum:.2f}'])
        power_data.append(dict(module=name,conditions_within_2pp=n,conditions=len(q),max_abs_change_pp=maximum))
    total=int((power.Delta_L_pp.abs()<=2).sum())
    rows.append(['All modules',f'{total} / {len(power)}',f'{power.Delta_L_pp.abs().max():.2f}'])
    headers=['Module','Conditions within\n2 pp of STC','Largest absolute\nchange (pp)']
    note=('Changes are relative to the loss at 1000 W/m² and 25 °C for the same band configuration. '
          'Each module has 56 conditions at 200, 400, 600 and 800 W/m². '
          'Bands are 80 mm rectangles and triangles with transmittance 0.71, 0.44 or 0 in both orientations. '
          'A 60 mm rectangle with transmittance 0.44 is also included in both orientations. '
          'pp denotes percentage points.')
    title='Sensitivity of relative power losses to irradiance'
    latex_table(out/'table_S5_power_irradiance.tex',title,'tab:irradiance-power',headers,rows,
                note.replace('W/m²',r'W m$^{-2}$').replace('°C',r'$^\circ$C'))
    preview(out/'table_S5_power_irradiance',title,headers,rows,[.42,.29,.29],note)
    pd.DataFrame(power_data).to_csv(out/'table_S5_power_irradiance.csv',index=False)

    annual=pd.read_csv(annual_path)
    annual=annual[annual.irradiance_basis=='after_uniform_layer']
    rows=[]
    annual_data=[]
    for site in ('Neuchâtel','Almería','Doha'):
        q=annual[annual.site==site]
        assert len(q)==8
        original=float(q.original_band_loss_pct.median())
        corrected=float(q.corrected_band_loss_pct.median())
        delta=corrected-original
        largest=float(q.annual_difference_pp.abs().max())
        rows.append([site,f'{original:.3f}',f'{corrected:.3f}',f'{delta:+.3f}',f'{largest:.3f}'])
        annual_data.append(dict(site=site,original_median_pct=original,corrected_median_pct=corrected,
                                change_in_median_pp=delta,largest_abs_individual_change_pp=largest))
    headers=['Site','STC estimate\n(%)','With irradiance\ndependence (%)',
             'Change in\nmedian (pp)','Largest individual\nchange (pp)']
    note=('Values are additional annual energy losses due to the band, relative to clean annual production. '
          'Medians cover four modules and two orientations. The largest individual change is the maximum '
          'absolute change among these eight combinations. The scenario uses a 60 mm rectangular band, '
          'f = 0.1, no manual cleaning and complete band removal after rainfall events of at least 15 mm. '
          'The correction isolates irradiance dependence at 25 °C. pp denotes percentage points.')
    title='Effect of irradiance on additional annual band energy losses'
    latex_table(out/'table_S6_annual_irradiance.tex',title,'tab:irradiance-annual',headers,rows,
                note.replace('°C',r'$^\circ$C'))
    preview(out/'table_S6_annual_irradiance',title,headers,rows,[.19,.18,.23,.18,.22],note)
    pd.DataFrame(annual_data).to_csv(out/'table_S6_annual_irradiance.csv',index=False)
    (out/'README.md').write_text(
        'Use Table S6 for the annual energy paragraph. Table S5 is a compact optional summary '
        'of the instantaneous power test accompanying Fig. S6. The numbering is suggested '
        'and follows the surrounding manuscript discussion.\n\n'
        'LaTeX snippets require booktabs. The PDF and PNG previews use the same rounded '
        'values. CSV files retain full precision. Horizontal rules only, no vertical '
        'borders or coloured cells. All statistics are derived from the simulation CSVs.\n\n'
        'The change in the site median is the corrected median minus the original median. '
        'It is not the median of the paired changes. The final column is an absolute '
        'maximum, not a confidence interval.\n',encoding='utf-8')
    (out/'sources.json').write_text(json.dumps({str(p.relative_to(config.REPO_DIR)):
        hashlib.sha256(p.read_bytes()).hexdigest() for p in (power_path,annual_path)},indent=2)+'\n')
    print(pd.DataFrame(annual_data).to_string(index=False))


if __name__=='__main__':
    main()
