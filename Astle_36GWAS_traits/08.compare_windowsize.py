import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,pearsonr


def scatterplot(df, tool1=None, tool2=None, tool1metric=None, tool2metric=None, metric=None, tool1level=None, tool2level=None, threshold=None):
    spearmanr_coefficient,spearmanr_pval = spearmanr(df[tool1metric], df[tool2metric]) 
    pearsonr_coefficient,pearsonr_pval = pearsonr(df[tool1metric], df[tool2metric]) 
    y = np.linspace(min(df[tool1metric]), max(df[tool1metric]), 100) 
    plt.figure(figsize=(4.5, 4.5))
    

    mask_highlight = ((df[tool1metric] > threshold) & (df[tool2metric] < threshold)) | ((df[tool1metric] < threshold) & (df[tool2metric] > threshold))


    plt.scatter(df[~mask_highlight][tool1metric], df[~mask_highlight][tool2metric], marker='.', color='black')  
    plt.scatter(df[mask_highlight][tool1metric], df[mask_highlight][tool2metric], marker='.', color='royalblue')  
    
    if tool1 == tool2 == 'COLOC':
        plt.text(0.88, -0.05, '0.75', color='crimson', verticalalignment='bottom', horizontalalignment='right', fontsize=12)
        plt.text(0.08, 0.81, '0.75', color='crimson', verticalalignment='top', horizontalalignment='right', fontsize=12)
    elif tool1 == tool2 == 'fastENLOC':
        plt.text(0.60, -0.04, '0.5', color='crimson', verticalalignment='bottom', horizontalalignment='right', fontsize=12)
        plt.text(0.045, 0.56, '0.5', color='crimson', verticalalignment='top', horizontalalignment='right', fontsize=12)
    else:
        plt.text(0.19, -0.04, '0.1', color='crimson', verticalalignment='bottom', horizontalalignment='right', fontsize=12)
        plt.text(0.045, 0.15, '0.1', color='crimson', verticalalignment='top', horizontalalignment='right', fontsize=12)
        

    #plt.plot(y, y, label='y = x', color='#00FF00', linestyle='--')
    plt.plot(y, y, label='y = x', color='lightgray', linestyle='--')

    plt.axvline(x=threshold, color='y', linestyle='--')
    plt.axhline(y=threshold, color='sienna', linestyle='--')

    pearsonr_results = f"Pearson's r" + f"={round(pearsonr_coefficient,3)} p-value = " + "{:.2e}".format(pearsonr_pval)
    spearman_results = f"Spearmanr's " + r"$\rho$" + f"={round(spearmanr_coefficient,3)} p-value = " + "{:.2e}".format(spearmanr_pval)
    
    #plt.xlabel(f"{tool1} {tool1metric}\n{pearsonr_results}\n{spearman_results}", fontsize=14)
    plt.xlabel(f"{tool1level}\n{spearman_results}", fontsize=12)
    plt.ylabel(f"{tool2level}", fontsize=12)

    plt.title(f"{tool1} {metric}", fontsize=13)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(f"/home/e1124850/github/locuscompare2_private/revise_20241015/{tool1}_{tool1level}_vs_{tool2}_{tool2level}.pdf",format='pdf', bbox_inches='tight')


# coloc gene vs coloc loci
coloc_gene = pd.read_csv('/home/e1124850/project/locuscompare2file/full_coloc_astle/GCST004599.full.coloc.tsv',sep='\t')
coloc_gene.index = coloc_gene['gene_id']
coloc_loci = pd.read_csv('/home/e1124850/project/locuscompare2file/Astle_blood_trait/Astle_blood_trait_results/GCST004599/coloc_output_20240518223806.tsv.gz',sep='\t')
coloc_loci.index = coloc_loci['gene_id']
coloc_loci = coloc_loci.drop_duplicates('gene_id')
df = pd.concat([coloc_gene['overall_H4'], coloc_loci['overall_H4']], join='inner', axis=1)
df.columns = ["gene-level overall_H4", "loci-level overall_H4"]
scatterplot(df, tool1='COLOC', tool2='COLOC', tool1metric='gene-level overall_H4', tool2metric='loci-level overall_H4', metric='overall_H4', tool1level='eQTL-level', tool2level='GWAS-level', threshold=0.75)

# fastenloc gene vs fastenloc loci

fastenloc_gene = pd.read_csv('~/project/locuscompare2file/fastenlocv3.1_astle_gene_level_nosetupprior/GCST004599_withoutcolocprior.enloc.gene.out',sep='\s+')
fastenloc_gene.index = fastenloc_gene['Gene']
fastenloc_loci = pd.read_csv('/home/e1124850/project/locuscompare2file/fastenlocv3gwasloci/GCST004599/processed/default/GCST004599/Whole_Blood/EUR/fastenloc/analyzed/fastenloc_output_20240923213145.tsv.gz',sep='\t')
fastenloc_loci = fastenloc_loci.drop_duplicates('gene_id')
fastenloc_loci.index = fastenloc_loci['gene_id']
df = pd.concat([fastenloc_gene['GRCP'], fastenloc_loci['GRCP']], join='inner', axis=1)
df.columns = ["gene-level GRCP", "loci-level GRCP"]
scatterplot(df, tool1='fastENLOC', tool2='fastENLOC', tool1metric='gene-level GRCP', tool2metric='loci-level GRCP', metric='GRCP', tool1level='eQTL-level', tool2level='GWAS-level', threshold = 0.5)


# ecaviar gene vs ecaviar loci
ecaviar_gene = pd.read_csv('/home/e1124850/project/locuscompare2file/Astle_blood_trait/GCST004599_ecaviar_genelevel_20241025/processed/default/GCST004599/Whole_Blood/EUR/ecaviar/analyzed/ecaviar_output_20241026211608.tsv.gz',sep='\t')
ecaviar_gene.index = ecaviar_gene['gene_id']
ecaviar_loci = pd.read_csv('/home/e1124850/project/locuscompare2file/Astle_blood_trait/Astle_blood_trait_results/GCST004599/ecaviar_output_20240521213302.tsv.gz',sep='\t')
ecaviar_loci.index = ecaviar_loci['gene_id']
ecaviar_loci = ecaviar_loci.drop_duplicates('gene_id')
df = pd.concat([ecaviar_gene['clpp'], ecaviar_loci['clpp']], join='inner', axis=1)
df.columns = ["gene-level CLPP", "loci-level CLPP"]
scatterplot(df, tool1='eCAVIAR', tool2='eCAVIAR', tool1metric='gene-level CLPP', tool2metric='loci-level CLPP', metric='CLPP', tool1level='eQTL-level', tool2level='GWAS-level', threshold=0.1)