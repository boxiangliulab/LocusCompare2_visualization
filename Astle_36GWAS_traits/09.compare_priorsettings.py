import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

def scatterplot(df, tool1=None, tool2=None, tool1metric=None, tool2metric=None, tool1level=None, tool2level=None, tool1prior=None, tool2prior=None):

    spearmanr_coefficient, spearmanr_pval = spearmanr(df[tool1metric], df[tool2metric]) 
    pearsonr_coefficient, pearsonr_pval = pearsonr(df[tool1metric], df[tool2metric]) 
    
    y = np.linspace(min(df[tool1metric]), max(df[tool1metric]), 100)
    
    plt.figure(figsize=(4.5, 4.5)) 

    mask_highlight = (df[tool1metric] > 0.75) & (df[tool2metric] > 0.5)

    plt.scatter(df[~mask_highlight][tool1metric], df[~mask_highlight][tool2metric], marker='.', color='black') 
    plt.scatter(df[mask_highlight][tool1metric], df[mask_highlight][tool2metric], marker='.', color='peru') 

    plt.plot(y, y, label='y = x', color='lightgray', linestyle='--')

    plt.axvline(x=0.75, color='y', linestyle='--', label='x = 0.75')
    plt.axhline(y=0.5, color='sienna', linestyle='--', label='y = 0.5')
    
    plt.text(0.88, -0.05, '0.75', color='crimson', verticalalignment='bottom', horizontalalignment='right', fontsize=12)
    plt.text(0.04, 0.56, '0.5', color='crimson', verticalalignment='top', horizontalalignment='right', fontsize=12)

    pearsonr_results = f"Pearson's r = {round(pearsonr_coefficient, 3)} p-value = " + "{:.2e}".format(pearsonr_pval)
    spearman_results = f"Spearman's " + r"$\rho$" + f" = {round(spearmanr_coefficient, 3)} p-value = " + "{:.2e}".format(spearmanr_pval)
    
    plt.xlabel(f"{tool1} {tool1metric}\n{pearsonr_results}\n{spearman_results}", fontsize=12)
    plt.ylabel(f"{tool2} {tool2metric}", fontsize=12)
    plt.title(f"{tool1} {tool1prior} vs {tool2} {tool2prior}", fontsize=13)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.savefig(f"/home/e1124850/github/locuscompare2_private/revise_20241015/{tool1}_{tool1prior}_vs_{tool2}_{tool2prior}.pdf", format='pdf', bbox_inches='tight')


coloc = pd.read_csv('/home/e1124850/project/locuscompare2file/full_coloc_astle/GCST004599.full.coloc.tsv',sep='\t')
coloc.index = coloc['gene_id']
fastenloc = pd.read_csv('/home/e1124850/project/locuscompare2file/fastenlocv3.1_astle_gene_level_nosetupprior/GCST004599_withoutcolocprior.enloc.gene.out',sep='\s+')
fastenloc = fastenloc.drop_duplicates('Gene')
fastenloc.index = fastenloc['Gene']
df = pd.concat([coloc['overall_H4'], fastenloc['GRCP']], join='inner', axis=1)
scatterplot(df, tool1='COLOC', tool2='fastENLOC', tool1metric='overall_H4', tool2metric='GRCP', tool1prior='default_prior', tool2prior='nospecified_prior')


coloc = pd.read_csv('/home/e1124850/project/locuscompare2file/full_coloc_astle/GCST004599.full.coloc.tsv',sep='\t')
coloc.index = coloc['gene_id']
coloc = coloc.drop_duplicates('gene_id')
fastenloc = pd.read_csv('/home/e1124850/project/locuscompare2file/fastenlocv3.1_astle_gene_level_colocdefaultprior/GCST004599_withcolocdefaultprior.enloc.gene.out',sep='\s+')
fastenloc = fastenloc.drop_duplicates('Gene')
fastenloc.index = fastenloc['Gene']
df = pd.concat([coloc['overall_H4'], fastenloc['GRCP']], join='inner', axis=1)
scatterplot(df, tool1='COLOC', tool2='fastENLOC', tool1metric='overall_H4', tool2metric='GRCP', tool1prior='default_prior', tool2prior='coloc_default_prior')