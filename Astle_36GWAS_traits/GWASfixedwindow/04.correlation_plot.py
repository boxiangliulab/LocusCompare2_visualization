import glob
import pandas as pd
import os
import scipy.stats as ss
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path

gwas_session_ls = \
    ['GCST004599','GCST004600','GCST004601','GCST004602','GCST004603','GCST004604',
     'GCST004605','GCST004606','GCST004607','GCST004608','GCST004609','GCST004610',
     'GCST004611','GCST004613','GCST004614','GCST004612','GCST004615','GCST004616',
     'GCST004617','GCST004618','GCST004619','GCST004620','GCST004621','GCST004622',
     'GCST004623','GCST004624','GCST004625','GCST004626','GCST004627','GCST004628',
     'GCST004629','GCST004630','GCST004631','GCST004632','GCST004633','GCST004634']


coloc_df = pd.DataFrame()
fastenloc_df = pd.DataFrame()
ecaviar_df = pd.DataFrame()
smr_df = pd.DataFrame()
predixcan_df = pd.DataFrame()
fusion_df = pd.DataFrame()


for gwas_session in gwas_session_ls:
    coloc_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/coloc/analyzed/coloc*')
    if len(coloc_paths) != 0:
        coloc = pd.read_csv(coloc_paths[0], sep='\t')
        coloc = coloc.sort_values('overall_H4',ascending=False)
        #coloc = coloc.drop_duplicates(['gene_id', 'gwas_lead_snp'],keep='first')
        #coloc.index = coloc['gene_id'] + coloc['gwas_lead_snp']+gwas_session
        coloc = coloc.drop_duplicates(['gene_id'],keep='first')
        coloc.index = coloc['gene_id'] + gwas_session
    
    fastenloc_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/fastenloc/analyzed/fastenloc_outp*')
    fastenloc = pd.read_csv(fastenloc_paths[0], sep='\t')
    fastenloc = fastenloc.sort_values('GRCP',ascending=False)
    #fastenloc = fastenloc.drop_duplicates(['gene_id', 'lead_variant'],keep='first')
    #fastenloc.index = fastenloc['gene_id'] + fastenloc['lead_variant']+gwas_session
    fastenloc = fastenloc.drop_duplicates(['gene_id'],keep='first')
    fastenloc.index = fastenloc['gene_id'] + gwas_session
    
    ecaviar_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/ecaviar/analyzed/ecaviar*')
    ecaviar = pd.read_csv(ecaviar_paths[0], sep='\t')
    ecaviar = ecaviar.sort_values('clpp',ascending=False)
    #ecaviar = ecaviar.drop_duplicates(['gene_id', 'lead_variant'],keep='first')
    #ecaviar.index = ecaviar['gene_id'] + ecaviar['lead_variant']+gwas_session
    ecaviar = ecaviar.drop_duplicates(['gene_id'],keep='first')
    ecaviar.index = ecaviar['gene_id'] + gwas_session
    
    smr_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/smr/analyzed/smr*')
    smr = pd.read_csv(smr_paths[0], sep='\t')
    smr.index = smr['gene_id']+gwas_session
    smr = smr.sort_values('p_SMR',ascending=True)
    smr = smr.drop_duplicates(['gene_id', 'topSNP'],keep='first')
    smr.rename(columns={'p_SMR': 'p_smr', 'b_SMR':'b_smr', 'se_SMR':'se_smr'}, inplace=True)

    predixcan_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/predixcan/analyzed/predixcan*')
    predixcan = pd.read_csv(predixcan_paths[0], sep='\t')
    predixcan.index = predixcan['gene_id']+gwas_session
    predixcan = predixcan.sort_values('pvalue',ascending=True)
    predixcan = predixcan.drop_duplicates(['gene_id'],keep='first')
    predixcan.rename(columns={'pvalue': 'p_predixcan', 'zscore':'z_predixcan', 'effect_size':'b_predixcan'}, inplace=True)
    
    
    fusion_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/fusion/analyzed/fusion*')
    fusion = pd.read_csv(fusion_paths[0], sep='\t')
    fusion.index = fusion['gene_id']+gwas_session
    fusion = fusion.sort_values('TWAS.P',ascending=True)
    fusion = fusion.drop_duplicates(['gene_id'],keep='first')
    fusion.rename(columns={'TWAS.P': 'p_fusion', 'TWAS.Z':'z_fusion'}, inplace=True)
    
    
    
    coloc_df = pd.concat([coloc_df, coloc[['overall_H4','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf']]], axis=0)
    fastenloc_df = pd.concat([fastenloc_df, fastenloc['GRCP']], axis=0)
    ecaviar_df = pd.concat([ecaviar_df, ecaviar['clpp']], axis=0)
    smr_df = pd.concat([smr_df, smr[['p_smr','b_smr','se_smr','p_HEIDI']]], axis=0)
    predixcan_df = pd.concat([predixcan_df, predixcan[['p_predixcan','z_predixcan','b_predixcan']]], axis=0)
    fusion_df = pd.concat([fusion_df, fusion[['p_fusion','z_fusion']]], axis=0)
    
    
coloc_df = coloc_df.sort_values('overall_H4',ascending=False)
fastenloc_df = fastenloc_df.sort_values('GRCP',ascending=False)
ecaviar_df = ecaviar_df.sort_values('clpp',ascending=False)
smr_df = smr_df.sort_values('p_smr',ascending=True)
predixcan_df = predixcan_df.sort_values('p_predixcan',ascending=True)
fusion_df = fusion_df.sort_values('p_fusion',ascending=True)
    
coloc_df = coloc_df[~coloc_df.index.duplicated(keep='first')]
fastenloc_df = fastenloc_df[~fastenloc_df.index.duplicated(keep='first')]
ecaviar_df = ecaviar_df[~ecaviar_df.index.duplicated(keep='first')]
smr_df = smr_df[~smr_df.index.duplicated(keep='first')]
predixcan_df = predixcan_df[~predixcan_df.index.duplicated(keep='first')]
fusion_df = fusion_df[~fusion_df.index.duplicated(keep='first')]

df = pd.concat([coloc_df['overall_H4'], fastenloc_df['GRCP'], ecaviar_df['clpp'],fusion_df['p_fusion'],predixcan_df['p_predixcan'], smr_df['p_smr']],join='inner',axis=1)
df.columns = ['coloc','fastenloc','ecaviar','fusion','predixcan','smr']
df['smr'] = -np.log10(df['smr'])
df['fusion'] = -np.log10(df['fusion'])
df['predixcan'] = -np.log10(df['predixcan'])

tmp = df['predixcan'].replace([np.inf, -np.inf], -1)
df['predixcan'] = tmp.replace([-1], np.max(tmp))

tmp = df['fusion'].replace([np.inf, -np.inf], -1)
df['fusion'] = tmp.replace([-1], np.max(tmp))

#fastenloc_df.columns = gwas_session_ls
coloc_df['ix'] = coloc_df.index
coloc_df['gene_id'] = coloc_df['ix'].apply(lambda x: str(x).split('GCST')[0])
coloc_df['gwas'] = coloc_df['ix'].apply(lambda x: 'GCST'+str(x).split('GCST')[1])
fastenloc_df['ix'] = fastenloc_df.index
fastenloc_df['gene_id'] = fastenloc_df['ix'].apply(lambda x: str(x).split('GCST')[0])
fastenloc_df['gwas'] = fastenloc_df['ix'].apply(lambda x: 'GCST'+str(x).split('GCST')[1])
ecaviar_df['ix'] = ecaviar_df.index
ecaviar_df['gene_id'] = ecaviar_df['ix'].apply(lambda x: str(x).split('GCST')[0])
ecaviar_df['gwas'] = ecaviar_df['ix'].apply(lambda x: 'GCST'+str(x).split('GCST')[1])
fusion_df['ix'] = fusion_df.index
fusion_df['gene_id'] = fusion_df['ix'].apply(lambda x: str(x).split('GCST')[0])
fusion_df['gwas'] = fusion_df['ix'].apply(lambda x: 'GCST'+str(x).split('GCST')[1])
predixcan_df['ix'] = predixcan_df.index
predixcan_df=predixcan_df.dropna()
predixcan_df['gene_id'] = predixcan_df['ix'].apply(lambda x: str(x).split('GCST')[0])
predixcan_df['gwas'] = predixcan_df['ix'].apply(lambda x: 'GCST'+str(x).split('GCST')[1])
smr_df['ix'] = smr_df.index
smr_df['gene_id'] = smr_df['ix'].apply(lambda x: str(x).split('GCST')[0])
smr_df['gwas'] = smr_df['ix'].apply(lambda x: 'GCST'+str(x).split('GCST')[1])



output_dir = '/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window'
Path(output_dir).mkdir(exist_ok=True, parents=True)
coloc_df.to_csv(os.path.join(output_dir, 'coloc_df.tsv.gz'),sep='\t',compression='gzip', index=False)
fastenloc_df.to_csv(os.path.join(output_dir, 'fastenloc_df.tsv.gz'),sep='\t',compression='gzip', index=False)
ecaviar_df.to_csv(os.path.join(output_dir, 'ecaviar_df.tsv.gz'),sep='\t',compression='gzip', index=False)
fusion_df.to_csv(os.path.join(output_dir, 'fusion_df.tsv.gz'),sep='\t',compression='gzip', index=False)
predixcan_df.to_csv(os.path.join(output_dir, 'predixcan_df.tsv.gz'),sep='\t',compression='gzip', index=False)
smr_df.to_csv(os.path.join(output_dir, 'smr_df.tsv.gz'),sep='\t',compression='gzip', index=False)


df.to_csv(os.path.join(output_dir, 'integrate_df_prob_log10p_20250218.tsv'),sep='\t')

integrate_df = pd.read_csv(os.path.join(output_dir, 'integrate_df_prob_log10p_20250218.tsv'),sep='\t')

def hide_current_axis(*args, **kwds):
    ax = plt.gca()
    ax.set_xticks([])  # 移除 X 轴刻度
    ax.set_yticks([])  # 移除 Y 轴刻度
    ax.set_xticklabels([])  # 隐藏 X 轴刻度数值
    ax.set_yticklabels([])  # 隐藏 Y 轴刻度数值
    ax = plt.gca().set_visible(False)
    sns.despine(ax=ax, bottom=True, top=True, left=True, right=True)
    
# 数据处理
integrate_dropna = integrate_df.dropna()
tmp = integrate_dropna.replace([np.inf, -np.inf], np.nan)
tmp = tmp.iloc[:, 1:]
integrate_dropna = integrate_dropna.replace([np.inf, -np.inf], np.max(tmp))
integrate_dropna.columns = ['ix', 'COLOC', 'fastENLOC', 'eCAVIAR', 'FUSION', 'PrediXcan', 'SMR']

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
    #ax.grid(False)
    ax.scatter(np.max(x)/2, np.max(y)/2, s=r*10000, color=facecolor, edgecolors='silver')
    ax.annotate(f"{r:.2f}", xy=(.5, .5), xycoords=ax.transAxes, weight="bold",
                color='black' if lightness < 0.7 else 'black', size=25, ha='center', va='center')
    ax.set_xticklabels([])  # 隐藏 x 轴刻度数值
    ax.set_yticklabels([])  # 隐藏 y 轴刻度数值
    
# 创建 PairGrid
g = sns.PairGrid(integrate_dropna[[ 'eCAVIAR', 'fastENLOC', 'COLOC', 'SMR', 'PrediXcan', 'FUSION']], height=1.8, aspect=0.95, corner=True)
g.map_lower(corrfunc, cmap=plt.get_cmap('PuOr'), norm=plt.Normalize(vmin=0, vmax=1))


g = g.add_legend(fontsize='25')
g.map_diag(hide_current_axis)

#g.axes.flat[0].set_xlim(0, 1)
#g.axes.flat[0].set_ylim(0, 1)
#sns.set_style("whitegrid", {'axes.grid': False})

g.savefig(os.path.join(output_dir, 'correlation_heatmap_20250320.pdf'),format='pdf', bbox_inches='tight')
plt.close()

plt.figure()
sns.set(font_scale=1.5)
sns.clustermap(integrate_df.iloc[:,1:].dropna().corr(method='spearman'), annot=True, cmap="crest",annot_kws={"size": 20})
plt.savefig(os.path.join(output_dir, 'heatmap_dropna_spearman_20250320.pdf'),format='pdf', bbox_inches='tight')


