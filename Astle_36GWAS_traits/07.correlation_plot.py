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

# coloc_df = pd.DataFrame()
# ecaviar_df = pd.DataFrame()
# fastenloc_df = pd.DataFrame()
# smr_df = pd.DataFrame()
# predixcan_df = pd.DataFrame()
# fusion_df = pd.DataFrame()

# for gwas_session in gwas_session_ls:
#     result_dir = f"/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/{gwas_session}/"
#     files = os.listdir(result_dir)

#     for file in files:
#         if file[:3] == 'smr':
#             smr = pd.read_csv(os.path.join(result_dir, file), sep='\t')
#         elif file[:5] == 'coloc':
#             coloc = pd.read_csv(os.path.join(result_dir, file), sep='\t')
#         elif file[:9] == 'predixcan':
#             predixcan = pd.read_csv(os.path.join(result_dir, file), sep='\t')
#         elif file[:6] == 'fusion':
#             fusion = pd.read_csv(os.path.join(result_dir, file), sep='\t')
#         elif file[:7] == 'ecaviar':
#             ecaviar = pd.read_csv(os.path.join(result_dir, file), sep='\t')
#         elif file[:9] == 'fastenloc':
#             fastenloc = pd.read_csv(os.path.join(result_dir, file), sep='\t')

#     coloc.index = coloc['gene_id']
#     smr.index = smr['gene_id']
#     fusion.index = fusion['gene_id']
#     predixcan.index = predixcan['gene_id']
#     ecaviar.index = ecaviar['gene_id']
#     fastenloc.index = fastenloc['gene_id']

#     coloc = coloc.sort_values('overall_H4',ascending=False)
#     coloc = coloc.drop_duplicates('gene_id',keep='first')
#     ecaviar = ecaviar.sort_values('clpp',ascending=False)
#     ecaviar = ecaviar.drop_duplicates('gene_id')

#     coloc_df = pd.concat([coloc_df, coloc['overall_H4']], axis=1)
#     ecaviar_df = pd.concat([ecaviar_df, ecaviar['clpp']], axis=1)
#     fastenloc_df = pd.concat([fastenloc_df, fastenloc['GLCP']], axis=1)
#     smr_df = pd.concat([smr_df, smr['p_SMR']], axis=1)
#     predixcan_df = pd.concat([predixcan_df, predixcan['pvalue']], axis=1)
#     fusion_df = pd.concat([fusion_df, fusion['TWAS.P']], axis=1)


# coloc_df.columns = gwas_session_ls
# ecaviar_df.columns = gwas_session_ls
# fastenloc_df.columns = gwas_session_ls
# smr_df.columns = gwas_session_ls
# predixcan_df.columns = gwas_session_ls
# fusion_df.columns = gwas_session_ls

# coloc_df.to_csv('/home/liulab/codes/locuscompare2_private/Astle_blood_trait/coloc_total_withna.tsv',sep='\t',index=True)
# ecaviar_df.to_csv('/home/liulab/codes/locuscompare2_private/Astle_blood_trait/ecaviar_total_withna.tsv',sep='\t',index=True)
# fastenloc_df.to_csv('/home/liulab/codes/locuscompare2_private/Astle_blood_trait/fastenloc_total_withna.tsv',sep='\t',index=True)
# smr_df.to_csv('/home/liulab/codes/locuscompare2_private/Astle_blood_trait/smr_total_withna.tsv',sep='\t',index=True)
# predixcan_df.to_csv('/home/liulab/codes/locuscompare2_private/Astle_blood_trait/predixcan_total_withna.tsv',sep='\t',index=True)
# fusion_df.to_csv('/home/liulab/codes/locuscompare2_private/Astle_blood_trait/fusion_total_withna.tsv',sep='\t',index=True)


path = '/Users/phoebel/github/locuscompare2_private/Astle_blood_trait'
coloc = pd.read_csv(f'{path}/coloc_total_withna.tsv',sep='\t',index_col=0)
ecaviar = pd.read_csv(f'{path}/ecaviar_total_withna.tsv',sep='\t',index_col=0)
fastenloc = pd.read_csv(f'{path}/fastenloc_total_withna.tsv',sep='\t',index_col=0)
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

plt.figure()
plt.title("prob -log10(pval) Heatmap")
plt.rcParams['font.size'] = '20'
sns.set(font_scale=1.7)
sns.clustermap(integrate_df.corr(), annot=True, cmap="crest",annot_kws={"size": 20})
plt.show()

integrate_df.to_csv('/Users/phoebel/github/locuscompare2_private/Astle_blood_trait/Astle_result_integrate_df_prob_log10p.tsv',sep='\t')


integrate_df = pd.read_csv('/Users/phoebel/github/locuscompare2_private/Astle_blood_trait/Astle_result_integrate_df_prob_log10p.tsv',sep='\t')

def hide_current_axis(*args, **kwds):
    plt.gca().set_visible(False)
    
integrate_dropna = integrate_df.dropna()
tmp = integrate_dropna.replace([np.inf, -np.inf], np.nan)
tmp = tmp.iloc[:,1:]
integrate_dropna = integrate_dropna.replace([np.inf, -np.inf], np.max(tmp))
integrate_dropna.columns = ['ix','COLOC','eCAVIAR','fastENLOC','SMR','PrediXcan','FUSION']
def corrfunc(x, y, **kwds):
    cmap = kwds['cmap']
    norm = kwds['norm']
    ax = plt.gca()
    ax.tick_params(bottom=False, top=False, left=False, right=False)
    sns.despine(ax=ax, bottom=True, top=True, left=True, right=True)
    r, _ = pearsonr(x, y)
    facecolor = cmap(norm(r))
    print(facecolor)
    #ax.set_facecolor(facecolor)
    lightness = (max(facecolor[:3]) + min(facecolor[:3]) )
    ax.grid(False)

    ax.scatter(np.max(x)/2,np.max(y)/2, s=r*10000,color=facecolor,edgecolors='silver')
    #ax.scatter(.5, .5, s=r*10000,color=facecolor,edgecolors='silver')
    #ax.add_patch(plt.Circle((np.max(x)/2,np.max(y)/2), 0.2, color=facecolor))
    ax.annotate(f"r={r:.2f}", xy=(.5, .5), xycoords=ax.transAxes, weight="bold",
                color='white' if lightness < 0.7 else 'black', size=25, ha='center', va='center')
    #ax.annotate(f"{r:.2f}", xy=(.3, .47), xycoords=ax.transAxes)    

order = ['ecaviar', 'coloc', 'fastenloc', 'smr', 'predixcan', 'fusion']
lims = [(0, 1), (0, 1),(0, 1),(0, 100), (0, 250), (0, 250)]
#tick_inc = [0.5, 0.5, 1, 0.5]
mpl.rcParams["axes.labelsize"] = 25
g = sns.PairGrid(integrate_dropna[['eCAVIAR', 'COLOC', 'fastENLOC', 'SMR', 'PrediXcan', 'FUSION']],height=1.8)
g.map_lower(plt.scatter, s=40, marker='.', color='black', alpha=0.6).add_legend(fontsize= '15')
g.map_diag(sns.histplot,bins=25)
g.map_upper(corrfunc, cmap=plt.get_cmap('PuOr'), norm=plt.Normalize(vmin=0, vmax=1))
#g.fig.subplots_adjust(wspace=0.1, hspace=0.1) # equal spacing in both directions
sns.set_style("whitegrid", {'axes.grid' : False})
g = g.add_legend(fontsize='15')

#loc = matplotlib.ticker.MultipleLocator(0.5)
g.axes.flat[0].set_xlim(0,1)
g.axes.flat[0].set_ylim(0,1)
#plt.show()

g.savefig(f'/Users/phoebel/github/locuscompare2_private/figure/mainfigure/correlation_scatter_combine_gene_trait_pairs.png',format='png', dpi=300, bbox_inches='tight')
