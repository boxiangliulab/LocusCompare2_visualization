import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
from pathlib import Path

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

def matchsnp_preprocess(eqtl, gwas, smallwindow_min, smallwindow_max):
    matchsnp = pd.concat([eqtl[['chr','rs_id_dbSNP151_GRCh38p7', 'pval_nominal','slope','pos']],gwas[['pval','beta']]],
                        join='inner',axis=1)
    matchsnp.columns = ['chr', 'rsid', 'eqtl_p','eqtl_beta','position','gwas_p','gwas_beta']

    matchsnp_ = matchsnp[(matchsnp['position']>smallwindow_min) &(matchsnp['position']<smallwindow_max)]
    matchsnp_ = matchsnp_.sort_values('gwas_p')
    lead_snp = matchsnp_.iloc[0,:]['rsid']
    del matchsnp_


    matchsnp['eqtl_p'] = -np.log10(matchsnp['eqtl_p'])
    matchsnp['gwas_p'] = -np.log10(matchsnp['gwas_p'])

    chrom = matchsnp.iloc[0,:]['chr']
    shell_command_plink_execute = 'plink --silent --vcf {} --r2 --ld-snp {} --ld-window 99999 --ld-window-kb 500 --ld-window-r2 0 --out {}'
    vcf_output_dir = f'/home/users/nus/locuscompare2file/vcf_hg38/EUR/chr{chrom}.vcf.gz'
    output_ld_file = f'/home/users/nus/scratch/test'
    os.system(shell_command_plink_execute.format(vcf_output_dir, lead_snp, output_ld_file))
    ld = pd.read_csv(f'{output_ld_file}.ld',sep='\s+')
    ld.index = ld['SNP_B']
    matchsnp["R2"] = matchsnp.index.map(ld["R2"]).fillna(0)
    matchsnp.columns = ['chr', 'rsid', 'eqtl_p','eqtl_beta','position','gwas_p','gwas_beta','LD']
    matchsnp = matchsnp.fillna(0)
    matchsnp['LD'][lead_snp] = 1
    matchsnp['LD_color'] = matchsnp['LD'].apply(lambda x: color_map(x))

    matchsnp = matchsnp.sort_values('position')
    matchsnp['position'] = matchsnp['position']/1e6

    return lead_snp, matchsnp


def gwas_matchsnp(gwas_df, gwas_both, bete_df):
    rs_ls = list(bete_df.index)
    matchsnp = pd.concat([gwas_df.loc[rs_ls][['rsid','pos','beta']],
                        gwas_both.loc[rs_ls]['beta']],join='inner',axis=1)
    matchsnp.columns = ['rsid','position','gwasonly_beta','gwasboth_beta']

    matchsnp = pd.concat([matchsnp, bete_df[['LD','LD_color']]],join='inner',axis=1)
    return matchsnp

    

def locuszoomplot(gene, chrom, gwas_session, category, smallwindow_min, smallwindow_max, largewindow_min, largewindow_max, figpath, smallwindowcolor = 'cornflowerblue', largewindowcolor = 'tan'):
    try:
        gwas_df = pd.read_csv(\
            f'/home/users/nus/locuscompare2file/Astle_blood_trait/{gwasfiledict[gwas_session]}',
            sep='\t')
        gwas_df.index = gwas_df['rsid']
        
        gwas_df = gwas_df.sort_values('pval')
        
    
        color = ["#0c028b", "#87ceeb", "#186400", "#f8a402", "#f50703"]
        classes = ['<0.2','0.2-0.4','0.4-0.6','0.6-0.8','>0.8']
        colors = ListedColormap(color)
        eqtl = pd.read_csv(\
            f'/home/users/nus/scratch/astle_202503/preprocessed/eqtl/Whole_Blood/grouped/{chrom}/{gene}.tsv.gz',
            sep='\t')
        eqtl = eqtl.sort_values('pval_nominal')
        eqtl = eqtl.drop_duplicates('rs_id_dbSNP151_GRCh38p7')
        eqtl.index = eqtl['rs_id_dbSNP151_GRCh38p7']

        lead_snp, matchsnp_df = matchsnp_preprocess(eqtl, gwas_df, smallwindow_min, smallwindow_max)


        ################################################
        #                  1. locuszoom                #
        ################################################
        # gwas
        fig = plt.figure(figsize = (5.5, 4.5))
        ax1 = plt.subplot(211)
        ax1.scatter(matchsnp_df['position'],matchsnp_df['gwas_p'], 
                    c=matchsnp_df['LD_color'], cmap=colors, s=50, 
                    edgecolors='grey', linewidths=0.01)
        color_ = 'darkviolet'
        ax1.scatter(matchsnp_df.loc[lead_snp]['position'], 
                matchsnp_df.loc[lead_snp]['gwas_p'], c='darkviolet', 
                marker='D', s=60, edgecolors='black', linewidths=1)
        ax1.text(matchsnp_df.loc[lead_snp]['position'], 
                matchsnp_df.loc[lead_snp]['gwas_p'], f"{lead_snp}", fontsize=10)
        ## gwas ld-based loci test region
        ax1.axvline(x=smallwindow_min/1e6, color=smallwindowcolor, linestyle='--', linewidth=1.5)
        ax1.axvline(x=smallwindow_max/1e6, color=smallwindowcolor, linestyle='--', linewidth=1.5)
        ## global ld-based loci test region
        ax1.axvline(x=largewindow_min/1e6, color=largewindowcolor, linestyle='--', linewidth=1.5)
        ax1.axvline(x=largewindow_max/1e6, color=largewindowcolor, linestyle='--', linewidth=1.5)
        

        ax1.set_ylabel(f'{gwastraitname[gwas_session]}\n' + 'GWAS -log$_{10}$(${P}$)',fontsize=13)
        ax1.spines.right.set_visible(False)
        ax1.spines.top.set_visible(False)
    
        # locus eQTL
        ax2 = plt.subplot(212)
        ax2.scatter(matchsnp_df['position'],matchsnp_df['eqtl_p'], 
                    c=matchsnp_df['LD_color'], cmap=colors, s=50, 
                    edgecolors='grey', linewidths=0.01)
        ax2.scatter(matchsnp_df.loc[lead_snp]['position'], 
                    matchsnp_df.loc[lead_snp]['eqtl_p'], c='darkviolet', 
                    marker='D',s=60, edgecolors='black', linewidths=1)
        ax2.text(matchsnp_df.loc[lead_snp]['position'], 
                matchsnp_df.loc[lead_snp]['eqtl_p'], f"{lead_snp}", fontsize=10)
        ax2.set_xlabel(f'Chromosome {chrom} (MB)',fontsize=13)
        ## gwas ld-based loci test region
        ax2.axvline(x=smallwindow_min/1e6, color=smallwindowcolor, linestyle='--', linewidth=1.5)
        ax2.axvline(x=smallwindow_max/1e6, color=smallwindowcolor, linestyle='--', linewidth=1.5)
        ## global ld-based loci test region
        ax2.axvline(x=largewindow_min/1e6, color=largewindowcolor, linestyle='--', linewidth=1.5)
        ax2.axvline(x=largewindow_max/1e6, color=largewindowcolor, linestyle='--', linewidth=1.5)
        ax2.set_ylabel('Whole blood\neQTL -log$_{10}$(${P}$)',fontsize=13)

        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
        plt.savefig(os.path.join(figpath, f'locuszoom_{gene}_{gwas_session}_20250410.pdf'),
                    format='pdf',bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error: {e}")
        print(f"Gene: {gene}, chrom: {chrom}")

figpath = '/home/users/nus/locuscompare2_private/rebuttle_20250328'
Path(figpath).mkdir(exist_ok=True, parents=True)

## done
chrom = 5
gene = 'ENSG00000113504.20'
gwas_session = 'GCST004602'
category = 'both'
smallwindow_min=990592 # gwasld H4=1.0 != combinedld 112597 H4=1.0
smallwindow_max=1170208 # gwasld = combinedld
largewindow_min= 784710  # fixed
largewindow_max=1783741# fixed
locuszoomplot(gene, chrom, gwas_session, category, smallwindow_min, smallwindow_max, largewindow_min, largewindow_max, figpath, smallwindowcolor = 'cornflowerblue', largewindowcolor = 'mediumspringgreen')

## done
gene = 'ENSG00000113504.20'
gwas_session = 'GCST004616'
category = 'both'
smallwindow_min= 990592# gwasld H4 = 0.999965 != combinedld 112597 H4 = 0.999965
smallwindow_max=1170208 # gwasld = combinedld 1170208
largewindow_min=604949 # fixed
largewindow_max=1604597# fixed
locuszoomplot(gene, chrom, gwas_session, category, smallwindow_min, smallwindow_max, largewindow_min, largewindow_max, figpath, smallwindowcolor = 'cornflowerblue', largewindowcolor = 'mediumspringgreen')



## done
chrom = 3
gene = 'ENSG00000161202.17'  ## DVL3
gwas_session = 'GCST004603'
category = 'both'
smallwindow_min=183563886 # gwasld # H4 = 0.840992040807663
smallwindow_max=184101069 # gwasld # 184101069
largewindow_min= 183155714  # combinedld 183549375 H4 = 0.840728
largewindow_max= 184097063 # combinedld 184101069
locuszoomplot(gene, chrom, gwas_session, category, smallwindow_min, smallwindow_max, largewindow_min, largewindow_max, figpath, smallwindowcolor = 'cornflowerblue', largewindowcolor = 'tan')

## done
chrom = 17
gene = 'ENSG00000204345.1' # CD300LD
gwas_session = 'GCST004607'
category = 'both'
smallwindow_min=74627502 # gwasld # H4 = 0.78433399092216
smallwindow_max=74796149 # gwasld # 74796149
largewindow_min= 74557354 # combinedld H4 = 0.000005
largewindow_max= 74806149 # combinedld
locuszoomplot(gene, chrom, gwas_session, category, smallwindow_min, smallwindow_max, largewindow_min, largewindow_max, figpath, smallwindowcolor = 'cornflowerblue', largewindowcolor = 'tan')

