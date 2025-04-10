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

    return matchsnp


def gwas_matchsnp(gwas_df, gwas_both, bete_df):
    rs_ls = list(bete_df.index)
    matchsnp = pd.concat([gwas_df.loc[rs_ls][['rsid','pos','beta']],
                        gwas_both.loc[rs_ls]['beta']],join='inner',axis=1)
    matchsnp.columns = ['rsid','position','gwasonly_beta','gwasboth_beta']

    matchsnp = pd.concat([matchsnp, bete_df[['LD','LD_color']]],join='inner',axis=1)
    return matchsnp

    

def locuszoomplot(gene, chrom, ld, gwas_session, category):
    try:
        gwas_df = pd.read_csv(\
            f'/work/lc_data/datasource/Astle_blood_trait/{gwasfiledict[gwas_session]}',
            sep='\t')
        gwas_df.index = gwas_df['rsid']
        
    
        color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]
        classes = ['<0.2','0.2-0.4','0.4-0.6','0.6-0.8','>0.8']
        colors = ListedColormap(color)
        gene_ = gene.split('.')[0]
        eqtl = pd.read_csv(\
            f'/work/lc_data/datasource/preprocessed/eqtl/Whole_Blood/grouped/{chrom}/{gene_}.tsv.gz',
            sep='\t')
        eqtl = eqtl.sort_values('pvalue')
        eqtl.index = eqtl['rsid']
        lead_snp = eqtl.iloc[0,:]['rsid']
        matchsnp_df = matchsnp_preprocess(eqtl, gwas_df, lead_snp, ld)


        ################################################
        #                  1. locuszoom                #
        ################################################
        fig = plt.figure(figsize = (8, 4.5))
        ax1 = plt.subplot(211)
        ax1.scatter(matchsnp_df['position'],matchsnp_df['gwas'], 
                    c=matchsnp_df['LD_color'], cmap=colors, s=50, 
                    edgecolors='black', linewidths=1)
        color_ = 'darkviolet'
        ax1.scatter(matchsnp_df.loc[lead_snp]['position'], 
                matchsnp_df.loc[lead_snp]['gwas'], c='darkviolet', 
                marker='D', s=60, edgecolors='black', linewidths=1)
        ax1.text(matchsnp_df.loc[lead_snp]['position'], 
                matchsnp_df.loc[lead_snp]['gwas'], f"{lead_snp}", fontsize=10)

        ax1.set_ylabel(f'{gwastraitname[gwas_session]}\n' + 'GWAS -log$_{10}$(${P}$)',fontsize=13)
        ax1.spines.right.set_visible(False)
        ax1.spines.top.set_visible(False)
    
        # locus eQTL
        ax2 = plt.subplot(212)
        ax2.scatter(matchsnp_df['position'],matchsnp_df['eqtl'], 
                    c=matchsnp_df['LD_color'], cmap=colors, s=50, 
                    edgecolors='black', linewidths=1)
        ax2.scatter(matchsnp_df.loc[lead_snp]['position'], 
                    matchsnp_df.loc[lead_snp]['eqtl'], c='darkviolet', 
                    marker='D',s=60, edgecolors='black', linewidths=1)
        ax2.text(matchsnp_df.loc[lead_snp]['position'], 
                matchsnp_df.loc[lead_snp]['eqtl'], f"{lead_snp}", fontsize=10)
        ax2.set_xlabel(f'Chromosome {chrom} (MB)',fontsize=13)
        # plt.legend(*sc2.legend_elements("sizes", num=6))
        ax2.set_ylabel('Whole blood\neQTL -log$_{10}$(${P}$)',fontsize=13)
        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/suppfigure/locuszoom_{gene}_{category}_{gwas_session}.png',
                    format='png',bbox_inches='tight')
        plt.close()
        

        ################################################
        #                  2. effect size              #
        ################################################
        fig = plt.figure(figsize = (4.5, 4.5))
        ax = plt.subplot(111)
        
        ax.scatter(matchsnp_df['eqtl_beta'], matchsnp_df['gwas_beta'], c=matchsnp_df['LD_color'],
                cmap=colors, s = 50, edgecolors='black', linewidths=1)

        color_ = 'darkviolet'
        ax.scatter(matchsnp_df.loc[lead_snp]['eqtl_beta'], 
                matchsnp_df.loc[lead_snp]['gwas_beta'], c='darkviolet', 
                marker='D',s=60, edgecolors='black', linewidths=1)
        ax.text(matchsnp_df.loc[lead_snp]['eqtl_beta'], 
                matchsnp_df.loc[lead_snp]['gwas_beta'], f"{lead_snp}", fontsize=10)
        
        ax.set_ylabel(f'{gwastraitname[gwas_session]}\nGWAS effect size',fontsize=12)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        pears = stats.pearsonr(matchsnp_df['eqtl_beta'],matchsnp_df['gwas_beta'])
        #annot_text = f"{only_annotation_text}"
        #plt.text(min(bete_df['eqtl_beta']), max(bete_df['gwas_beta']), f"{annot_text}")
        annot_text = f"Pearson's $\\it{{r}}$ = {pears.correlation.round(3)}, P-value={'%.3g' % pears.pvalue}"
        #plt.text(min(bete_df['eqtl_beta']), min(bete_df['gwas_beta']), f"{annot_text}")
        ax.set_xlabel('Whole blood eQTL effect size'+f"\n{annot_text}",fontsize=12)
        x = np.array(matchsnp_df['eqtl_beta']).reshape((-1, 1))
        y = np.array(matchsnp_df['gwas_beta'])
        model = LinearRegression().fit(x, y)
        y_pred = model.predict(x)
        ax.plot(x,y_pred,color='silver',linestyle='dashed')
        plt.savefig(f'/home/liulab/codes/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/suppfigure/effectsize_{gene}_{category}_{gwas_session}_fdrthreshold.png',
                    format='png',bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error: {e}")
        print(f"Gene: {gene}, chrom: {chrom}")

        

df = pd.read_csv('/home/liulab/codes/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/coloc_twas_signal_summary_20250224.tsv',sep='\t')



coloc_only = df[df['category'] == 'colocalization_only']
both = df[df['category'] == 'both']
twas_only = df[df['category'] == 'twas_only']

twas_only['twas'] = twas_only[['smr','predixcan','fusion']].mean(axis=1)
twas_only = twas_only.sort_values('twas', ascending=False)
twas_only = twas_only.drop_duplicates('gene')


chrom_ls = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
for chrom in chrom_ls:
    twas_only_ = twas_only[twas_only['chrom'] == chrom]
    if len(twas_only_) >0:
        ld = pd.read_csv(f'/work/lc_data/datasource/ld/tkg_p3v5a_ld_chr{chrom}_EUR.csv.gz',sep='\t',header=None)
        for ix, row in twas_only_.iterrows():
            gene = row['gene']
            gwas_session = row['gwas_session']
            category = row['category']
            locuszoomplot(gene, chrom, ld, gwas_session, category)


for chrom in chrom_ls:
    coloc_only_ = coloc_only[coloc_only['chrom'] == chrom]
    if len(coloc_only_) >0:
        ld = pd.read_csv(f'/work/lc_data/datasource/ld/tkg_p3v5a_ld_chr{chrom}_EUR.csv.gz',sep='\t',header=None)
        for ix, row in coloc_only_.iterrows():
            gene = row['gene']
            gwas_session = row['gwas_session']
            category = row['category']
            locuszoomplot(gene, chrom, ld, gwas_session, category)


for chrom in chrom_ls:
    both_ = both[both['chrom'] == chrom]
    if len(both_) >0:
        ld = pd.read_csv(f'/work/lc_data/datasource/ld/tkg_p3v5a_ld_chr{chrom}_EUR.csv.gz',sep='\t',header=None)
        for ix, row in both_.iterrows():
            gene = row['gene']
            gwas_session = row['gwas_session']
            category = row['category']
            locuszoomplot(gene, chrom, ld, gwas_session, category)
            
            
            
df = pd.read_csv('/home/liulab/codes/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/notwas_4599_20250225.tsv.gz',sep='\t')


for chrom in chrom_ls:
    df_ = df[df['chrom'] == chrom]
    if len(df_) >0:
        ld = pd.read_csv(f'/work/lc_data/datasource/ld/tkg_p3v5a_ld_chr{chrom}_EUR.csv.gz',sep='\t',header=None)
        for ix, row in df_.iterrows():
            gene = row['gene']
            gwas_session = row['gwas_session']
            category = 'notwas'
            locuszoomplot(gene, chrom, ld, gwas_session, category)