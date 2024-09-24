import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import matplotlib as mpl
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os
from scipy.stats import spearmanr


gwas_session_ls = \
    ['GCST004599',
    'GCST004600',
    'GCST004601',
    'GCST004602',
    'GCST004603',
    'GCST004604',
    'GCST004605',
    'GCST004606',
    'GCST004607',
    'GCST004608',
    'GCST004609',
    'GCST004610',
    'GCST004611',
    'GCST004613',
    'GCST004614',
    'GCST004612',
    'GCST004615',
    'GCST004616',
    'GCST004617',
    'GCST004618',
    'GCST004619',
    'GCST004620',
    'GCST004621',
    'GCST004622',
    'GCST004623',
    'GCST004624',
    'GCST004625',
    'GCST004626',
    'GCST004627',
    'GCST004628',
    'GCST004629',
    'GCST004630',
    'GCST004631',
    'GCST004632',
    'GCST004633',
    'GCST004634']

path = '/Users/phoebel/github/locuscompare2_private/Astle_blood_trait'
coloc = pd.read_csv(f'{path}/coloc_total_withna.tsv',sep='\t',index_col=0)
ecaviar = pd.read_csv(f'{path}/ecaviar_total_withna.tsv',sep='\t',index_col=0)
fastenloc = pd.read_csv(f'{path}/fastenlocv3.1_gwas_loci_total_withna_20240923.tsv',sep='\t',index_col=0)
smr = pd.read_csv(f'{path}/smr_total_withna.tsv',sep='\t',index_col=0)
predixcan = pd.read_csv(f'{path}/predixcan_total_withna.tsv',sep='\t',index_col=0)
fusion = pd.read_csv(f'{path}/fusion_total_withna.tsv',sep='\t',index_col=0)

coloc_gene_trait_ls = []
coloc_result = []
ecaviar_gene_trait_ls = []
ecaviar_result = []
fastenloc_gene_trait_ls = []
fastenloc_result = []
smr_gene_trait_ls = []
smr_result = []
predixcan_gene_trait_ls = []
predixcan_result = []
fusion_gene_trait_ls = []
fusion_result = []

for ix, row in coloc.iterrows():
  for trait in gwas_session_ls:
    coloc_gene_trait_ls.append(f"{ix}_{trait}")
    coloc_result.append(row[trait])
    
for ix, row in ecaviar.iterrows():
  for trait in gwas_session_ls:
    ecaviar_gene_trait_ls.append(f"{ix}_{trait}")
    ecaviar_result.append(row[trait])
    
for ix, row in smr.iterrows():
  for trait in gwas_session_ls:
    smr_gene_trait_ls.append(f"{ix}_{trait}")
    smr_result.append(row[trait])
    
for ix, row in fastenloc.iterrows():
  for trait in gwas_session_ls:
    fastenloc_gene_trait_ls.append(f"{ix}_{trait}")
    fastenloc_result.append(row[trait])
    
for ix, row in predixcan.iterrows():
  for trait in gwas_session_ls:
    predixcan_gene_trait_ls.append(f"{ix}_{trait}")
    predixcan_result.append(row[trait])
    
for ix, row in fusion.iterrows():
  for trait in gwas_session_ls:
    fusion_gene_trait_ls.append(f"{ix}_{trait}")
    fusion_result.append(row[trait])
    
    
coloc_integrate = pd.DataFrame({'ix':coloc_gene_trait_ls,'coloc':coloc_result})
ecaviar_integrate = pd.DataFrame({'ix':ecaviar_gene_trait_ls,'ecaviar':ecaviar_result})
fastenloc_integrate = pd.DataFrame({'ix':fastenloc_gene_trait_ls,'fastenloc':fastenloc_result})
smr_integrate = pd.DataFrame({'ix':smr_gene_trait_ls,'smr':smr_result})
predixcan_integrate = pd.DataFrame({'ix':predixcan_gene_trait_ls,'predixcan':predixcan_result})
fusion_integrate = pd.DataFrame({'ix':fusion_gene_trait_ls,'fusion':fusion_result})

coloc_integrate.index = coloc_integrate['ix']
ecaviar_integrate.index = ecaviar_integrate['ix']
fastenloc_integrate.index = fastenloc_integrate['ix']
smr_integrate.index = smr_integrate['ix']
predixcan_integrate.index = predixcan_integrate['ix']
fusion_integrate.index = fusion_integrate['ix']


integrate_df = pd.concat([coloc_integrate['coloc'],ecaviar_integrate['ecaviar'],fastenloc_integrate['fastenloc'],smr_integrate['smr'],predixcan_integrate['predixcan'],fusion_integrate['fusion']],axis=1)


integrate_df[['smr','predixcan','fusion']] = -np.log10(integrate_df[['smr','predixcan','fusion']])

tmp = integrate_df['predixcan'].replace([np.inf, -np.inf], -1)
integrate_df['predixcan'] = integrate_df['predixcan'].replace([-1], np.max(tmp))

tmp = integrate_df['fusion'].replace([np.inf, -np.inf], -1)
integrate_df['fusion'] = integrate_df['fusion'].replace([-1], np.max(tmp))

#plt.figure()
#plt.title("prob -log10(pval) Heatmap")
#plt.rcParams['font.size'] = '20'
#sns.set(font_scale=1.7)
#sns.clustermap(integrate_df.corr(method='spearman'), annot=True, cmap="crest",annot_kws={"size": 20})
#plt.show()

integrate_df.to_csv('/Users/phoebel/github/locuscompare2_private/Astle_blood_trait/Astle_result_integrate_df_prob_log10p_20240923.tsv',sep='\t')


integrate_df = pd.read_csv('/Users/phoebel/github/locuscompare2_private/Astle_blood_trait/Astle_result_integrate_df_prob_log10p_20240923.tsv',sep='\t')

def hide_current_axis(*args, **kwds):
    plt.gca().set_visible(False)
    
# 数据处理
integrate_dropna = integrate_df.dropna()
tmp = integrate_dropna.replace([np.inf, -np.inf], np.nan)
tmp = tmp.iloc[:, 1:]
integrate_dropna = integrate_dropna.replace([np.inf, -np.inf], np.max(tmp))
integrate_dropna.columns = ['ix', 'COLOC', 'eCAVIAR', 'fastENLOC', 'SMR', 'PrediXcan', 'FUSION']

# 相关性函数
def corrfunc(x, y, **kwds):
    cmap = kwds['cmap']
    norm = kwds['norm']
    ax = plt.gca()
    ax.tick_params(bottom=False, top=False, left=False, right=False)
    sns.despine(ax=ax, bottom=True, top=True, left=True, right=True)
    r, _ = spearmanr(x, y)
    facecolor = cmap(norm(r))
    lightness = (max(facecolor[:3]) + min(facecolor[:3]))
    ax.grid(False)
    ax.scatter(np.max(x)/2, np.max(y)/2, s=r*10000, color=facecolor, edgecolors='silver')
    ax.annotate(r'$\rho$' + f"={r:.2f}", xy=(.5, .5), xycoords=ax.transAxes, weight="bold",
                color='black' if lightness < 0.7 else 'black', size=25, ha='center', va='center')

# 创建 PairGrid
g = sns.PairGrid(integrate_dropna[['eCAVIAR', 'COLOC', 'fastENLOC', 'SMR', 'PrediXcan', 'FUSION']], height=1.8, aspect=0.95)
g.map_lower(plt.scatter, s=40, marker='.', color='black', alpha=0.6)
g.map_diag(sns.histplot, bins=25)
g.map_upper(corrfunc, cmap=plt.get_cmap('PuOr'), norm=plt.Normalize(vmin=0, vmax=1))

# 添加颜色条
#cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=plt.get_cmap('PuOr'), norm=plt.Normalize(vmin=0, vmax=1)), ax=g.axes, orientation='vertical', label='Correlation ρ', pad=0.01)
#cbar.ax.tick_params(labelsize=15)

# 手动设置颜色条位置
#cbar.ax.set_position([0.95, 0.1, 0.03, 0.8])  # [左, 下, 宽, 高]

# 设置坐标轴标签
for ax in g.axes[:, 0]:
    ax.set_ylabel('')
for ax in g.axes[-1, :]:
    ax.set_xlabel('')

# 设置样式和图例
sns.set_style("whitegrid", {'axes.grid': False})
g = g.add_legend(fontsize='15')

# 设置坐标轴限制
g.axes.flat[0].set_xlim(0, 1)
g.axes.flat[0].set_ylim(0, 1)

#plt.show()

g.savefig(f'/Users/phoebel/github/locuscompare2_private/figure/mainfigure/correlation_scatter_combine_gene_trait_pairs_20240923.png',format='png', dpi=900, bbox_inches='tight')

plt.figure()
sns.set(font_scale=1.5)
sns.clustermap(integrate_df.iloc[:,1:].dropna().corr(method='spearman'), annot=True, cmap="crest",annot_kws={"size": 20})
plt.savefig(f'/Users/phoebel/github/locuscompare2_private/Astle_blood_trait/heatmap_dropna_spearman_20240923.png',format='png', dpi=300, bbox_inches='tight')