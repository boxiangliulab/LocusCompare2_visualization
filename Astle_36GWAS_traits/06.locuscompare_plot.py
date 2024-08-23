import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats

gwasfiledict = {
    'GCST004599':'preprocessed_mpv_build37_164454_20161212.tsv.gz',
    'GCST004601':'preprocessed_rbc_build37_172952_20161212.tsv.gz',
    'GCST004600':'preprocessed_eo_p_build37_172378_20161212.tsv.gz',
    'GCST004602':'preprocessed_mcv_build37_172433_20161212.tsv.gz',
    'GCST004603':'preprocessed_plt_build37_166066_20161212.tsv.gz',
    'GCST004604':'preprocessed_hct_build37_173039_20161212.tsv.gz',
    'GCST004605':'preprocessed_mchc_build37_172851_20161212.tsv.gz',
    'GCST004606':'preprocessed_eo_build37_172275_20161212.tsv.gz',
    'GCST004607':'preprocessed_pct_build37_164339_20161212.tsv.gz',
    'GCST004608':'preprocessed_gran_p_myeloid_wbc_build37_169545_20161212.tsv.gz',
    'GCST004609':'preprocessed_mono_p_build37_170494_20161212.tsv.gz',
    'GCST004610':'preprocessed_wbc_build37_172435_20161212.tsv.gz',
    'GCST004611':'preprocessed_hlr_build37_170761_20161212.tsv.gz',
    'GCST004612':'preprocessed_hlr_p_build37_170763_20161212.tsv.gz',
    'GCST004613':'preprocessed_neut_eo_sum_build37_170384_20161212.tsv.gz',
    'GCST004614':'preprocessed_gran_build37_169822_20161212.tsv.gz',
    'GCST004615':'preprocessed_hgb_build37_172925_20161212.tsv.gz',
    'GCST004616':'preprocessed_pdw_build37_164433_20161212.tsv.gz',
    'GCST004617':'preprocessed_eo_p_gran_build37_170536_20161212.tsv.gz',
    'GCST004618':'preprocessed_baso_build37_171846_20161212.tsv.gz',
    'GCST004619':'preprocessed_ret_p_build37_170690_20161212.tsv.gz',
    'GCST004620':'preprocessed_baso_neut_sum_build37_170143_20161212.tsv.gz',
    'GCST004621':'preprocessed_rdw_build37_171529_20161212.tsv.gz',
    'GCST004622':'preprocessed_ret_build37_170641_20161212.tsv.gz',
    'GCST004623':'preprocessed_neut_p_gran_build37_170672_20161212.tsv.gz',
    'GCST004624':'preprocessed_eo_baso_sum_build37_171771_20161212.tsv.gz',
    'GCST004625':'preprocessed_mono_build37_170721_20161212.tsv.gz',
    'GCST004626':'preprocessed_myeloid_wbc_build37_169219_20161212.tsv.gz',
    'GCST004627':'preprocessed_lymph_build37_171643_20161212.tsv.gz',
    'GCST004628':'preprocessed_irf_build37_170548_20161212.tsv.gz',
    'GCST004629':'preprocessed_neut_build37_170702_20161212.tsv.gz',
    'GCST004630':'preprocessed_mch_build37_172332_20161212.tsv.gz',
    'GCST004631':'preprocessed_baso_p_build37_171996_20161212.tsv.gz',
    'GCST004632':'preprocessed_lymph_p_build37_171748_20161212.tsv.gz',
    'GCST004633':'preprocessed_neut_p_build37_171542_20161212.tsv.gz',
    'GCST004634':'preprocessed_baso_p_gran_build37_170223_20161212.tsv.gz',
}

gwastraitname = {
    "GCST004599": "Mean platelet volume",
    "GCST004601": "Red blood cell count",
    "GCST004600": "Eosinophil percentage of white cells",
    "GCST004602": "Mean corpuscular volume",
    "GCST004603": "Platelet count",
    "GCST004604": "Hematocrit",
    "GCST004605": "Mean corpuscular hemoglobin concentration",
    "GCST004606": "Eosinophil counts",
    "GCST004607": "Plateletcrit",
    "GCST004608": "Granulocyte percentage of myeloid white cells",
    "GCST004609": "Monocyte percentage of white cells",
    "GCST004610": "White blood cell count",
    "GCST004611": "High light scatter reticulocyte count",
    "GCST004612": "High light scatter reticulocyte percentage of red cells",
    "GCST004613": "Sum neutrophil eosinophil counts",
    "GCST004614": "Granulocyte count",
    "GCST004615": "Hemoglobin concentration",
    "GCST004616": "Platelet distribution width",
    "GCST004617": "Eosinophil percentage of granulocytes",
    "GCST004618": "White blood cell count (basophil)",
    "GCST004619": "Reticulocyte fraction of red cells",
    "GCST004620": "Sum basophil neutrophil counts",
    "GCST004621": "Red cell distribution width",
    "GCST004622": "Reticulocyte count",
    "GCST004623": "Neutrophil percentage of granulocytes",
    "GCST004624": "Sum eosinophil basophil counts",
    "GCST004625": "Monocyte count",
    "GCST004626": "Myeloid white cell count",
    "GCST004627": "Lymphocyte count",
    "GCST004628": "Immature fraction of reticulocytes",
    "GCST004629": "Neutrophil count",
    "GCST004630": "Mean corpuscular hemoglobin",
    "GCST004631": "Basophil percentage of white cells",
    "GCST004632": "Lymphocyte percentage of white cells",
    "GCST004633": "Neutrophil percentage of white cells",
    "GCST004634": "Basophil percentage of granulocytes",
}

gene_symbol = pd.read_csv('/home/liulab/datasource/gencode.v26.basic.annotation.gtf.gz',
                          sep='\t',comment='#',header=None)
gene_symbol['gene_symbol'] = gene_symbol[8].apply(\
    lambda x: x.split('gene_name ')[1].split('"')[1])
gene_symbol['gene'] = gene_symbol[8].apply(\
    lambda x: x.split('gene_id ')[1].split('"')[1].split('.')[0])
gene_symbol = gene_symbol[gene_symbol[2] == 'gene']
gene_symbol.index = gene_symbol['gene']
gene_symbol = gene_symbol.drop_duplicates('gene')

color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]

def color_map(LD):
    color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]
    if LD < 0.2:
        return color[0]
    elif LD < 0.4:
        return color[1]
    elif LD < 0.6:
        return color[2]
    elif LD < 0.8:
        return color[3]
    else:
        return color[4]

def matchsnp_preprocess(eqtl, gwas, lead_snp, ld):
    matchsnp = pd.concat([eqtl[['pvalue','beta','position']],gwas[['pval','beta']]],
                         join='inner',axis=1)
    matchsnp.columns = ['eqtl','eqtl_beta','position','gwas','gwas_beta']
    matchsnp['eqtl'] = -np.log10(matchsnp['eqtl'])
    matchsnp['gwas'] = -np.log10(matchsnp['gwas'])

    # color by LD
    
    targetsnp_ld = ld[ld[1] == lead_snp].iloc[:,1:]
    targetsnp_ld.index = targetsnp_ld[2]
    targetsnp_ld = targetsnp_ld.loc[list(set(targetsnp_ld.index) & set(matchsnp.index))]
    matchsnp = pd.concat([matchsnp,targetsnp_ld[3]],axis=1)
    matchsnp.columns = ['eqtl','eqtl_beta','position','gwas','gwas_beta','LD']
    matchsnp = matchsnp.fillna(0)
    matchsnp['LD'][lead_snp] = 1
    matchsnp['LD_color'] = matchsnp['LD'].apply(lambda x: color_map(x))


    matchsnp = matchsnp.sort_values('position')
    pos = matchsnp.loc[lead_snp]['position']
    beta_true = []
    for ix, row in matchsnp.iterrows():
        if row['position'] > pos-1e5 and row['position'] < pos+1e5:
            beta_true.append('1')
        else:
            beta_true.append('0')

    matchsnp['beta_true'] = beta_true
    matchsnp['position'] = matchsnp['position']/1e6
    
    return matchsnp


def gwas_matchsnp(gwas_only, gwas_both, bete_df):
    rs_ls = list(bete_df.index)
    matchsnp = pd.concat([gwas_only.loc[rs_ls][['rsid','pos','beta']],
                          gwas_both.loc[rs_ls]['beta']],join='inner',axis=1)
    matchsnp.columns = ['rsid','position','gwasonly_beta','gwasboth_beta']

    matchsnp = pd.concat([matchsnp, bete_df[['LD','LD_color']]],join='inner',axis=1)
    return matchsnp

def find_coloc_top_causal_SNP(gwas_session, gene):
    result_folder = '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results'
    result_path = os.path.join(result_folder, gwas_session)
    files = os.listdir(result_path)
    for file in files:
        if file[:5] == 'coloc':
            coloc = pd.read_csv(os.path.join(result_path, file), sep='\t')
            break


    coloc = coloc[coloc['gene_id'] == gene]
    coloc = coloc.sort_values('SNP.PP.H4',ascending=False)
    snp_chr = coloc.iloc[0,:]['snp'].split('_')[0][3:]
    snp_pos =  coloc.iloc[0,:]['snp'].split('_')[1]
    return snp_chr, snp_pos
    

def locuscompareplot(gene, chrom, ld, outfolder, only_tmp, both_tmp):
    try:
        gwas_only_session = list(set(only_tmp['gwas_session']))[0]
        gwas_only = pd.read_csv(\
            f'/home/liulab/datasource/Astle_blood_trait/{gwasfiledict[gwas_only_session]}',
            sep='\t')
        gwas_only.index = gwas_only['rsid']
        
        gwas_both_session = list(set(both_tmp['gwas_session']))[0]
        gwas_both = pd.read_csv(\
            f'/home/liulab/datasource/Astle_blood_trait/{gwasfiledict[gwas_both_session]}',
            sep='\t')
        gwas_both.index = gwas_both['rsid']
        
        category = list(set(only_tmp['category']))[0]
    
        color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]
        classes = ['<0.2','0.2-0.4','0.4-0.6','0.6-0.8','>0.8']
        colors = ListedColormap(color)
        eqtl = pd.read_csv(\
            f'/home/liulab/datasource/preprocessed/eqtl/Whole_Blood/grouped/{chrom}/{gene}.tsv.gz',
            sep='\t')
        eqtl = eqtl.sort_values('pvalue')
        eqtl.index = eqtl['rsid']
        lead_snp = eqtl.iloc[0,:]['rsid']
        matchsnp_only = matchsnp_preprocess(eqtl, gwas_only, lead_snp, ld)
        matchsnp_both = matchsnp_preprocess(eqtl, gwas_both, lead_snp, ld)

        if category == 'colocalization_only':
            only_annotation_text = "colocalization signal only"
        elif category == 'twas_only':
            only_annotation_text = "TWAS/MR signal only"

        both_annotation_text = "All signals"		
        
        snp_chr, snp_pos  = find_coloc_top_causal_SNP(gwas_only_session, gene)
        topcoloc_snp_identify_by_gwas_only_session = \
        list(gwas_only[(gwas_only['pos'] == int(snp_pos)) & \
            (gwas_only['chr'] == int(snp_chr))]['rsid'])[0]

        snp_chr, snp_pos  = find_coloc_top_causal_SNP(gwas_both_session, gene)
        topcoloc_identify_by_gwas_both_session = \
        list(gwas_both[(gwas_both['pos'] == int(snp_pos)) & \
            (gwas_both['chr'] == int(snp_chr))]['rsid'])[0]

        
        if gene == 'ENSG00000245937':
            print('1')
            topcoloc_identify_by_gwas_both_session = 'rs2250127'
            topcoloc_snp_identify_by_gwas_only_session = 'rs2250127'
            
        
        ################################################
        #          1. locus coloc or TWAS only         #
        ################################################
        fig = plt.figure(figsize = (8, 4.5))
        ax1 = plt.subplot(211)
        ax1.scatter(matchsnp_only['position'],matchsnp_only['gwas'], 
                    c=matchsnp_only['LD_color'], cmap=colors, alpha=0.8, s=50, 
                    edgecolors='black', linewidths=1)
        if topcoloc_snp_identify_by_gwas_only_session == lead_snp:
            color_ = 'darkviolet'
            ax1.scatter(matchsnp_only.loc[lead_snp]['position'], 
                        matchsnp_only.loc[lead_snp]['gwas'], c='darkviolet', 
                        marker='D', s=60, edgecolors='black', linewidths=1)
            ax1.text(matchsnp_only.loc[lead_snp]['position'], 
                     matchsnp_only.loc[lead_snp]['gwas'], f"{lead_snp}")
        else:
            color_ = 'pink'
            ax1.scatter(matchsnp_only.loc[lead_snp]['position'], 
                        matchsnp_only.loc[lead_snp]['gwas'],c='darkviolet', 
                        marker='D',s=60, edgecolors='black', linewidths=1)
            ax1.text(matchsnp_only.loc[lead_snp]['position'], 
                     matchsnp_only.loc[lead_snp]['gwas'], f"{lead_snp}")
            ax1.scatter(matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['position'], 
                        matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['gwas'], 
                        c=color_,marker='D',s=60, edgecolors='black', linewidths=1)
            ax1.text(matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['position'], 
                     matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['gwas'], 
                     f"{topcoloc_snp_identify_by_gwas_only_session}")
        #plt.text(min(matchsnp_both['position']), max(matchsnp_only['gwas']), f"{only_annotation_text}")
        # plt.xlabel('Chromosome 7 (BP)')
        ax1.set_ylabel(f'{gwastraitname[gwas_only_session]}\n' + 'GWAS -log$_{10}$(${P}$)',fontsize=13)
        ax1.spines.right.set_visible(False)
        ax1.spines.top.set_visible(False)
    
        # locus eQTL
        ax2 = plt.subplot(212)
        ax2.scatter(matchsnp_only['position'],matchsnp_only['eqtl'], 
                    c=matchsnp_only['LD_color'], cmap=colors, alpha=0.8, s=50, 
                    edgecolors='black', linewidths=1)
        ax2.scatter(matchsnp_only.loc[lead_snp]['position'], 
                    matchsnp_only.loc[lead_snp]['eqtl'], c='darkviolet', 
                    marker='D',s=60, edgecolors='black', linewidths=1)
        ax2.text(matchsnp_only.loc[lead_snp]['position'], 
                 matchsnp_only.loc[lead_snp]['eqtl'], f"{lead_snp}")
        ax2.set_xlabel(f'Chromosome {chrom} (MB)',fontsize=13)
        # plt.legend(*sc2.legend_elements("sizes", num=6))
        ax2.set_ylabel('Whole blood\neQTL -log$_{10}$(${P}$)',fontsize=13)
        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_locus_{category}_{gwas_only_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()
        
        ################################################
        #         2. locus coloc and TWAS both         #
        ################################################
        fig = plt.figure(figsize = (8, 4.5))
        ax1 = plt.subplot(211)
        ax1.scatter(matchsnp_both['position'],matchsnp_both['gwas'], 
                    c=matchsnp_both['LD_color'],cmap=colors, alpha=0.8, s=50, 
                    edgecolors='black', linewidths=1)
        if topcoloc_identify_by_gwas_both_session == lead_snp:
            color_ = 'darkviolet'
            ax1.scatter(matchsnp_both.loc[lead_snp]['position'], 
                        matchsnp_both.loc[lead_snp]['gwas'], c='darkviolet', 
                        marker='D',s=60, edgecolors='black', linewidths=1)
            ax1.text(matchsnp_both.loc[lead_snp]['position'], 
                     matchsnp_both.loc[lead_snp]['gwas'], f"{lead_snp}")
        else:
            color_ = 'pink'
            ax1.scatter(matchsnp_both.loc[lead_snp]['position'], 
                        matchsnp_both.loc[lead_snp]['gwas'], c='darkviolet', 
                        marker='D',s=60, edgecolors='black', linewidths=1)
            ax1.text(matchsnp_both.loc[lead_snp]['position'], 
                     matchsnp_both.loc[lead_snp]['gwas'], f"{lead_snp}")
            ax1.scatter(matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['position'], 
                        matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['gwas'], 
                        c=color_,marker='D',s=60, edgecolors='black', linewidths=1)
            ax1.text(matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['position'], 
                     matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['gwas'], 
                     f"{topcoloc_identify_by_gwas_both_session}")

        # plt.xlabel('Chromosome 7 (BP)')
        #plt.text(min(matchsnp_both['position']), max(matchsnp_both['gwas']), f"{both_annotation_text}")
        ax1.set_ylabel(f'{gwastraitname[gwas_both_session]}\n' + 'GWAS -log$_{10}$(${P}$)',
                       fontsize=13)
        ax1.spines.right.set_visible(False)
        ax1.spines.top.set_visible(False)

        # locus eQTL
        ax2 = plt.subplot(212)
        ax2.scatter(matchsnp_only['position'],matchsnp_only['eqtl'], 
                    c=matchsnp_only['LD_color'], cmap=colors, alpha=0.8, s=50, 
                    edgecolors='black', linewidths=1)
        ax2.scatter(matchsnp_only.loc[lead_snp]['position'], 
                    matchsnp_only.loc[lead_snp]['eqtl'], c='darkviolet', marker='D', 
                    s=60, edgecolors='black', linewidths=1)
        ax2.set_xlabel(f'Chromosome {chrom} (MB)',fontsize=13)
        ax2.text(matchsnp_only.loc[lead_snp]['position'], 
                 matchsnp_only.loc[lead_snp]['eqtl'], f"{lead_snp}")
        # plt.legend(*sc2.legend_elements("sizes", num=6))
        ax2.set_ylabel('Whole blood\neQTL -log$_{10}$(${P}$)',fontsize=13)
        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_locus_both_{gwas_both_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()
        
        ################################################
        #         3. pvalue coloc or TWAS only         #
        ################################################
        fig = plt.figure(figsize = (4.5, 4.5))
        ax = plt.subplot(111)
        #ax.set_title(f"{gene_symbol.loc[gene]['gene_symbol']}\n", fontsize=13)
        plt.scatter(matchsnp_only['eqtl'],matchsnp_only['gwas'], 
                    c=matchsnp_only['LD_color'],cmap=colors, alpha=0.8, s = 50, 
                    edgecolors='black', linewidths=1)
        if topcoloc_snp_identify_by_gwas_only_session == lead_snp:
            color_ = 'darkviolet'
            ax.scatter(matchsnp_only.loc[lead_snp]['eqtl'], 
                       matchsnp_only.loc[lead_snp]['gwas'], c='darkviolet', 
                       marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(matchsnp_only.loc[lead_snp]['eqtl'], 
                    matchsnp_only.loc[lead_snp]['gwas'], f"{lead_snp}")
        else:
            color_ = 'pink'
            ax.scatter(matchsnp_only.loc[lead_snp]['eqtl'], 
                       matchsnp_only.loc[lead_snp]['gwas'], c='darkviolet', 
                       marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(matchsnp_only.loc[lead_snp]['eqtl'], 
                    matchsnp_only.loc[lead_snp]['gwas'], f"{lead_snp}")
            ax.scatter(matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['eqtl'], 
                       matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['gwas'], 
                       c=color_,marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['eqtl'], 
                    matchsnp_only.loc[topcoloc_snp_identify_by_gwas_only_session]['gwas'], 
                    f"{topcoloc_snp_identify_by_gwas_only_session}")

        
        ax.scatter(matchsnp_only.loc[lead_snp]['eqtl'], 
                   matchsnp_only.loc[lead_snp]['gwas'], c='darkviolet', 
                   marker='D',s=60, edgecolors='black', linewidths=1)
        ax.text(matchsnp_only.loc[lead_snp]['eqtl'], 
                matchsnp_only.loc[lead_snp]['gwas'], f"{lead_snp}")
        #ax.set_xlabel('Whole blood eQTL -log$_{10}$(${P}$)',fontsize=13)
        ax.set_ylabel(f'{gwastraitname[gwas_only_session]}\n'+'GWAS -log$_{10}$(${P}$)', 
                      fontsize=12)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        #plt.text(0, max(matchsnp_only['gwas']), f"{only_annotation_text}")
        spear = stats.pearsonr(matchsnp_only['eqtl'],matchsnp_only['gwas'])
        #annot_text = f"Spearman coefficient={spear.correlation.round(3)}, P-value={'%.3g' % spear.pvalue}"
        #plt.text(0, min(matchsnp_only['gwas']), f"{annot_text}")
        ax.set_xlabel('Whole blood eQTL -log$_{10}$(${P}$)'+f"\nSpearman coefficient={spear.correlation.round(3)}, P-value={'%.3g' % spear.pvalue}",
                      fontsize=12)
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_pvalue_{category}_{gwas_only_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()
        
        ################################################
        #       4. effect size coloc or TWAS only      #
        ################################################
        fig = plt.figure(figsize = (4.5, 4.5))
        ax = plt.subplot(111)
        bete_df = matchsnp_only[matchsnp_only['beta_true'] == '1']
        ax.scatter(bete_df['eqtl_beta'], bete_df['gwas_beta'], c=bete_df['LD_color'],
                   cmap=colors, alpha=0.8, s = 50, edgecolors='black', linewidths=1)
        if topcoloc_snp_identify_by_gwas_only_session == lead_snp:
            color_ = 'darkviolet'
            ax.scatter(bete_df.loc[lead_snp]['eqtl_beta'], 
                       bete_df.loc[lead_snp]['gwas_beta'], c='darkviolet', 
                       marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(bete_df.loc[lead_snp]['eqtl_beta'], 
                    bete_df.loc[lead_snp]['gwas_beta'], f"{lead_snp}")
        else:
            color_ = 'pink'
            ax.scatter(bete_df.loc[lead_snp]['eqtl_beta'], 
                    bete_df.loc[lead_snp]['gwas_beta'], c='darkviolet', 
                    marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(bete_df.loc[lead_snp]['eqtl_beta'], 
                    bete_df.loc[lead_snp]['gwas_beta'], f"{lead_snp}")
            if topcoloc_snp_identify_by_gwas_only_session in set(bete_df.index):
                ax.scatter(bete_df.loc[topcoloc_snp_identify_by_gwas_only_session]['eqtl_beta'], 
                        bete_df.loc[topcoloc_snp_identify_by_gwas_only_session]['gwas_beta'], 
                        c=color_,marker='D',s=60, edgecolors='black', linewidths=1)
                ax.text(bete_df.loc[topcoloc_snp_identify_by_gwas_only_session]['eqtl_beta'], 
                        bete_df.loc[topcoloc_snp_identify_by_gwas_only_session]['gwas_beta'], 
                        f"{topcoloc_snp_identify_by_gwas_only_session}")
        
        ax.set_ylabel(f'{gwastraitname[gwas_only_session]}\nGWAS effect size',fontsize=12)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        pears = stats.pearsonr(bete_df['eqtl_beta'],bete_df['gwas_beta'])
        #annot_text = f"{only_annotation_text}"
        #plt.text(min(bete_df['eqtl_beta']), max(bete_df['gwas_beta']), f"{annot_text}")
        annot_text = f"Pearson coefficient={pears.correlation.round(3)}, P-value={'%.3g' % pears.pvalue}"
        #plt.text(min(bete_df['eqtl_beta']), min(bete_df['gwas_beta']), f"{annot_text}")
        ax.set_xlabel('Whole blood eQTL effect size'+f"\n{annot_text}",fontsize=12)
        x = np.array(bete_df['eqtl_beta']).reshape((-1, 1))
        y = np.array(bete_df['gwas_beta'])
        model = LinearRegression().fit(x, y)
        y_pred = model.predict(x)
        ax.plot(x,y_pred,color='silver',linestyle='dashed')
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_effectsize_{category}_{gwas_only_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()
        
        ################################################
        #        5. pvalue coloc and TWAS both         #
        ################################################
        fig = plt.figure(figsize = (4.5, 4.5))
        ax = plt.subplot(111)
        ax.scatter(matchsnp_both['eqtl'],matchsnp_both['gwas'], 
                   c=matchsnp_both['LD_color'],cmap=colors, alpha=0.8, s = 50, 
                   edgecolors='black', linewidths=1)
        if topcoloc_identify_by_gwas_both_session == lead_snp:
            color_ = 'darkviolet'
            ax.scatter(matchsnp_both.loc[lead_snp]['eqtl'], 
                       matchsnp_both.loc[lead_snp]['gwas'], c='darkviolet', 
                       marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(matchsnp_both.loc[lead_snp]['eqtl'], 
                    matchsnp_both.loc[lead_snp]['gwas'], f"{lead_snp}")
        else:
            color_ = 'pink'
            ax.scatter(matchsnp_both.loc[lead_snp]['eqtl'], 
                       matchsnp_both.loc[lead_snp]['gwas'], c='darkviolet', 
                       marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(matchsnp_both.loc[lead_snp]['eqtl'], 
                    matchsnp_both.loc[lead_snp]['gwas'], f"{lead_snp}")
            ax.scatter(matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['eqtl'], 
                       matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['gwas'],
                       c=color_,marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['eqtl'],
                    matchsnp_both.loc[topcoloc_identify_by_gwas_both_session]['gwas'], 
                    f"{topcoloc_identify_by_gwas_both_session}")

        #plt.xlabel('Whole blood eQTL -log$_{10}$(${P}$)',fontsize=13)
        ax.set_ylabel(f'{gwastraitname[gwas_both_session]}\n' + 'GWAS -log$_{10}$(${P}$)',
                      fontsize=12)
        #plt.text(0, max(matchsnp_both['gwas']), f"{both_annotation_text}")
        spear = stats.spearmanr(matchsnp_both['eqtl'],matchsnp_both['gwas'])
        annot_text = f"Spearman coefficient={spear.correlation.round(3)}, P-value={'%.3g' % spear.pvalue}"
        #plt.text(0, min(matchsnp_both['gwas']), f"{annot_text}")
        ax.set_xlabel('Whole blood eQTL -log$_{10}$(${P}$)'+f"\nSpearman coefficient={spear.correlation.round(3)}, P-value={'%.3g' % spear.pvalue}",
                      fontsize=12)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_pvalue_both_{gwas_both_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()
        
        ################################################
        #      6. effect size coloc and TWAS both      #
        ################################################
        fig = plt.figure(figsize = (4.5, 4.5))
        ax = plt.subplot(111)
        bete_df = matchsnp_both[matchsnp_both['beta_true'] == '1']
        ax.scatter(bete_df['eqtl_beta'], bete_df['gwas_beta'], c=bete_df['LD_color'], 
                   cmap=colors, alpha=0.8, s = 50, edgecolors='black', linewidths=1)
        if topcoloc_snp_identify_by_gwas_only_session == lead_snp:
            color_ = 'darkviolet'
            ax.scatter(bete_df.loc[lead_snp]['eqtl_beta'], 
                       bete_df.loc[lead_snp]['gwas_beta'], c='darkviolet', 
                       marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(bete_df.loc[lead_snp]['eqtl_beta'], 
                    bete_df.loc[lead_snp]['gwas_beta'], f"{lead_snp}")
        else:
            color_ = 'pink'
            ax.scatter(bete_df.loc[lead_snp]['eqtl_beta'], 
                    bete_df.loc[lead_snp]['gwas_beta'], c='darkviolet', marker='D',
                    s=60, edgecolors='black', linewidths=1)
            ax.text(bete_df.loc[lead_snp]['eqtl_beta'], 
                    bete_df.loc[lead_snp]['gwas_beta'], f"{lead_snp}")
            if topcoloc_identify_by_gwas_both_session in set(bete_df.index):
                ax.scatter(bete_df.loc[topcoloc_identify_by_gwas_both_session]['eqtl_beta'], 
                        bete_df.loc[topcoloc_identify_by_gwas_both_session]['gwas_beta'], 
                        c='pink',marker='D',s=60, edgecolors='black', linewidths=1)
                ax.text(bete_df.loc[topcoloc_identify_by_gwas_both_session]['eqtl_beta'], 
                        bete_df.loc[topcoloc_identify_by_gwas_both_session]['gwas_beta'], 
                        f"{topcoloc_identify_by_gwas_both_session}")
                

        #plt.xlabel('Whole blood eQTL effect size',fontsize=13)
        ax.set_ylabel(f'{gwastraitname[gwas_both_session]}\nGWAS effect size',fontsize=12)
        pears = stats.pearsonr(bete_df['eqtl_beta'],bete_df['gwas_beta'])
        #annot_text = f"{both_annotation_text}"
        #plt.text(min(bete_df['eqtl_beta']), max(bete_df['gwas_beta']), f"{annot_text}")
        annot_text = f"Pearson coefficient={pears.correlation.round(3)}, P-value={'%.3g' % pears.pvalue}"
        #plt.text(min(bete_df['eqtl_beta']), min(bete_df['gwas_beta']), f"{annot_text}")
        ax.set_xlabel('Whole blood eQTL effect size'+f"\n{annot_text}",fontsize=12)
        #plt.text(0.0, min(bete_df['gwas_beta'])+0.015, f"Pearson coefficient={pears.correlation.round(4)}")
        #plt.text(0.0, min(bete_df['gwas_beta']), f"P-value={'%.3g' % pears.pvalue}")
        x = np.array(bete_df['eqtl_beta']).reshape((-1, 1))
        y = np.array(bete_df['gwas_beta'])
        model = LinearRegression().fit(x, y)
        y_pred = model.predict(x)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.plot(x,y_pred,color='silver',linestyle='dashed')
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_effectsize_both_{gwas_both_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()


        ################################################
        #     7. gwas only & both effect size plot     #
        ################################################
        gwas_beta = gwas_matchsnp(gwas_only, gwas_both, bete_df)
        fig = plt.figure(figsize = (4.5, 4.5))
        ax = plt.subplot(111)

        ax.scatter(gwas_beta['gwasonly_beta'], gwas_beta['gwasboth_beta'], 
                   c=gwas_beta['LD_color'],cmap=colors, alpha=0.8, s = 50, 
                   edgecolors='black', linewidths=1)
        if topcoloc_snp_identify_by_gwas_only_session in set(gwas_beta.index):
            ax.scatter(gwas_beta.loc[topcoloc_snp_identify_by_gwas_only_session]['gwasonly_beta'],
                    gwas_beta.loc[topcoloc_snp_identify_by_gwas_only_session]['gwasboth_beta'],
                    c='pink',marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(gwas_beta.loc[topcoloc_snp_identify_by_gwas_only_session]['gwasonly_beta'],
                    gwas_beta.loc[topcoloc_snp_identify_by_gwas_only_session]['gwasboth_beta'], 
                    f"{topcoloc_snp_identify_by_gwas_only_session}")
        
        if topcoloc_identify_by_gwas_both_session in set(gwas_beta.index):
            ax.scatter(gwas_beta.loc[topcoloc_identify_by_gwas_both_session]['gwasonly_beta'],
                    gwas_beta.loc[topcoloc_identify_by_gwas_both_session]['gwasboth_beta'],
                    c='pink',marker='D',s=60, edgecolors='black', linewidths=1)
            ax.text(gwas_beta.loc[topcoloc_identify_by_gwas_both_session]['gwasonly_beta'],
                    gwas_beta.loc[topcoloc_identify_by_gwas_both_session]['gwasboth_beta'], 
                    f"{topcoloc_identify_by_gwas_both_session}")
            
        
        ax.scatter(gwas_beta.loc[lead_snp]['gwasonly_beta'], 
                   gwas_beta.loc[lead_snp]['gwasboth_beta'], c='darkviolet',
                   marker='D', s=60, edgecolors='black', linewidths=1)
        
        ax.text(gwas_beta.loc[lead_snp]['gwasonly_beta'],
                gwas_beta.loc[lead_snp]['gwasboth_beta'], f"{lead_snp}")
        #plt.xlabel(f'{gwastraitname[gwas_only_session]} GWAS effect size',fontsize=13)
        ax.set_ylabel(f'{gwastraitname[gwas_both_session]}\nGWAS effect size',
                      fontsize=12)
        pears = stats.pearsonr(gwas_beta['gwasonly_beta'],gwas_beta['gwasboth_beta'])
        annot_text = f"Pearson coefficient={pears.correlation.round(3)}, P-value={'%.3g' % pears.pvalue}"
        ax.set_xlabel(f'{gwastraitname[gwas_only_session]} GWAS effect size'+f"\n{annot_text}",
                      fontsize=12)
        #plt.text(0.0, min(bete_df['gwas_beta'])+0.015, f"Pearson coefficient={pears.correlation.round(4)}")
        #plt.text(0.0, min(bete_df['gwas_beta']), f"P-value={'%.3g' % pears.pvalue}")
        x = np.array(gwas_beta['gwasonly_beta']).reshape((-1, 1))
        y = np.array(gwas_beta['gwasboth_beta'])
        model = LinearRegression().fit(x, y)
        y_pred = model.predict(x)
        ax.plot(x,y_pred,color='silver',linestyle='dashed')
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        #plt.savefig(f'~/datasource/test.pdf',format='pdf',bbox_inches='tight')

        plt.savefig(f'/home/liulab/codes/locuscompare2_private/figure/supp_figure/{outfolder}/{gene}_effectsize_{gwas_both_session}_{gwas_only_session}_fdrthreshold.pdf',
                    format='pdf',bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error: {e}")
        print(f"Gene: {gene}, chrom: {chrom}")


outfolder = 'illustrate_example'

#for gene in ['ENSG00000054523']:
for gene in ['ENSG00000245937']:
    only_tmp = coloc_only[coloc_only['gene'] == gene]
    both_tmp = both[both['gene'] == gene]
    chrom = list(set(only_tmp['chr']))[0]
    ld = pd.read_csv(f'/home/liulab/datasource/ld/tkg_p3v5a_ld_chr{chrom}_EUR.csv.gz',
                     sep='\t',header=None)
    print(f"Gene: {gene}")
    locuscompareplot(gene,chrom, ld, outfolder, only_tmp, both_tmp)

order = both.iloc[:,1:7].sum(axis=1)
both['order'] = order
both = both.sort_values('order',ascending=False)

order = twas_only.iloc[:,4:7].sum(axis=1)/twas_only.iloc[:,1:4].sum(axis=1)
twas_only['order'] = order
twas_only = twas_only.sort_values('order',ascending=False)

for gene in ['ENSG00000174529','ENSG00000187122']:
    only_tmp = twas_only[twas_only['gene'] == gene]
    both_tmp = both[both['gene'] == gene]
    chrom = list(set(only_tmp['chr']))[0]
    ld = pd.read_csv(f'/home/liulab/datasource/ld/tkg_p3v5a_ld_chr{chrom}_EUR.csv.gz',
                     sep='\t',header=None)
    print(f"Gene: {gene}")
    for _, row1 in only_tmp.iterrows():
        only_tmp_ = pd.DataFrame(row1).T
        for _, row2 in both_tmp.iterrows():
            both_tmp_ = pd.DataFrame(row2).T
            locuscompareplot(gene,chrom, ld, outfolder, only_tmp, both_tmp_)
            break
        break