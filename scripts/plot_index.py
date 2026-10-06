"""公开实验指数审计图；兼容范围不是统计置信区间。"""
import os
from pathlib import Path
cache=Path(__file__).resolve().parents[1]/'.cache'
cache.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(cache))
os.environ.setdefault('XDG_CACHE_HOME',str(cache))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import json
import pandas as pd
import numpy as np
from PIL import Image

plt.rcParams.update({
    'figure.dpi': 600.1, 'savefig.dpi': 600.1, 'font.family': 'Arial',
    'font.size': 10, 'axes.titlesize': 14, 'axes.labelsize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 10,
    'axes.linewidth': .5, 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.facecolor': 'white',
    'savefig.bbox': 'tight', 'savefig.pad_inches': .1,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'qa/figures'; OUT.mkdir(exist_ok=True)
COLORS=['#E64B35','#4DBBD5','#00A087','#3C5488','#7E6163','#F39B7F','#8491B4']
all_data=pd.read_csv(ROOT/'qa/all_7467_hierarchical_v04.csv')
representative=pd.read_csv(ROOT/'qa/representative_v04.csv')
sensitivity=pd.read_csv(ROOT/'qa/weight_sensitivity_v04.csv')
audit=json.loads((ROOT/'qa/index_audit_v04.json').read_text())
chosen=['R175H','R273H','Y220C','G245S','R249S','R282W']
fig=plt.figure(figsize=(7.2,10.2))
gs=fig.add_gridspec(3,1,height_ratios=[1.1,1,1.4],hspace=.68)
ax=fig.add_subplot(gs[0])
presets=['function_priority','equal','structure_priority']
for i,preset in enumerate(presets):
    for y,q in enumerate(chosen):
        row=sensitivity.loc[(sensitivity['query']==q)&(sensitivity['preset']==preset)].iloc[0]
        yy=y+(i-1)*.2
        ax.plot([row.total_compatible_low,row.total_compatible_high],[yy,yy],color=COLORS[i],lw=1.5)
        ax.plot([row.total_compatible_low,row.total_compatible_high],[yy,yy],'|',color=COLORS[i],ms=5,mew=1)
ax.set_yticks(range(len(chosen)),chosen);ax.invert_yaxis()
ax.set_xlim(0,100);ax.set_xlabel('Compatible index range (0–100)')
ax.set_title('A  Weight sensitivity of selected variants',loc='left')
ax.legend([Line2D([],[],color=COLORS[i],lw=1.5) for i in range(3)],
          ['S/F/D = 25/50/25','33/33/33','50/25/25'],
          ncol=1,frameon=False,loc='upper left',fontsize=9)

ax=fig.add_subplot(gs[1])
coverage=all_data.groupby('coverage').size()
labels=[f'{v*100:.1f}%' for v in coverage.index]
bars=ax.bar(range(len(coverage)),coverage.values,color=COLORS[:len(coverage)],width=.65)
for b,count in zip(bars,coverage.values):ax.text(b.get_x()+b.get_width()/2,b.get_height()+70,f'{count:,}',ha='center',fontsize=10)
ax.set_xticks(range(len(coverage)),labels)
ax.set_ylabel('Protein substitutions (N)')
ax.set_xlabel('Weighted quantitative coverage')
ax.set_ylim(0,max(coverage.values)*1.22)
ax.set_title('B  Quantitative coverage: N = 7,467',loc='left')

ax=fig.add_subplot(gs[2])
reportable=all_data.loc[all_data.reportable].sort_values(['total_compatible_low','total_compatible_high','query']).reset_index(drop=True)
for y,row in reportable.iterrows():
    color=COLORS[0] if row.structure_coverage>0 else COLORS[1]
    ax.plot([row.total_compatible_low,row.total_compatible_high],[y,y],color=color,lw=.45,alpha=.55)
    ax.plot(row.total_compatible_low,y,'|',color=color,ms=1.8,mew=.45)
    ax.plot(row.total_compatible_high,y,'|',color=color,ms=1.8,mew=.45)
ax.set_xlim(0,100);ax.set_ylim(-10,len(reportable)+10)
ax.set_xlabel('Default compatible index range (0–100)')
ax.set_ylabel('Variants ordered by lower bound')
ax.set_title(f'C  Display-eligible variants: N = {len(reportable):,}',loc='left')
ax.legend([Line2D([],[],color=COLORS[0],lw=1.5),Line2D([],[],color=COLORS[1],lw=1.5)],
          ['At least one structure axis measured','Structure axes unavailable'],
          frameon=False,loc='upper left',fontsize=9)
fig.text(.01,.005,'Ranges retain missing values and source ambiguity; they are not confidence intervals.\n'
         'Default index = 0.25 S + 0.50 F + 0.25 D; display rule: coverage ≥60% and ≥2 functional families.\n'
         'This descriptive research index is uncalibrated for pathogenicity, treatment response or clinical outcomes.',
         fontsize=8,ha='left',va='bottom')
stem=OUT/'Fig_v04_hierarchical_index_audit'
fig.savefig(stem.with_suffix('.png'),dpi=600.1,bbox_inches='tight')
fig.savefig(stem.with_suffix('.pdf'),bbox_inches='tight')
fig.savefig(stem.with_suffix('.svg'),bbox_inches='tight')
plt.close(fig)
img=Image.open(stem.with_suffix('.png'));dpi=img.info.get('dpi',(0,0))
assert min(dpi)>=600, dpi
assert audit['variants']==len(all_data)==7467
assert audit['reportable']==len(reportable)==1124
assert coverage.sum()==7467
for suffix in ['.png','.pdf']:
    assert stem.with_suffix(suffix).stat().st_size>=50*1024
qc={'png_dpi':dpi,'png_pixels':img.size,'png_bytes':stem.with_suffix('.png').stat().st_size,
    'pdf_bytes':stem.with_suffix('.pdf').stat().st_size,'all_variants':len(all_data),
    'reportable':len(reportable),'snapshot':audit['snapshot_id'],
    'ranges':'Missing/source compatibility envelope, not CI','rasterized':False}
(OUT/'figure_qc_v04.json').write_text(json.dumps(qc,indent=2)+'\n')
print(json.dumps(qc,indent=2))
