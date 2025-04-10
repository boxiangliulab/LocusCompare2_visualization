import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
from tqdm import tqdm
import seaborn as sns


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
total_smr = pd.DataFrame()
total_predixcan = pd.DataFrame()
total_fusion = pd.DataFrame()

gwas_session_ls = \
    ['GCST004599','GCST004600','GCST004601','GCST004602','GCST004603','GCST004604',
     'GCST004605','GCST004606','GCST004607','GCST004608','GCST004609','GCST004610',
     'GCST004611','GCST004613','GCST004614','GCST004612','GCST004615','GCST004616',
     'GCST004617','GCST004618','GCST004619','GCST004620','GCST004621','GCST004622',
     'GCST004623','GCST004624','GCST004625','GCST004626','GCST004627','GCST004628',
     'GCST004629','GCST004630','GCST004631','GCST004632','GCST004633','GCST004634']
    
for gwas_session in gwas_session_ls:
    smr_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/smr/analyzed/smr*')
    smr = pd.read_csv(smr_paths[0], sep='\t')
    smr['ix'] = smr['gene_id'] + gwas_session + '_chr_' + smr['chrom'].astype(str)
    total_smr = pd.concat([total_smr,smr],axis=0)

    predixcan_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/predixcan/analyzed/predixcan*')
    predixcan = pd.read_csv(predixcan_paths[0], sep='\t')
    predixcan = predixcan[['gene', 'gene_name', 'zscore', 'effect_size', 'pvalue','chrom', 'gene_id']]
    predixcan = predixcan.dropna()
    predixcan['ix'] = predixcan['gene_id']+gwas_session+'_chr_'+predixcan['chrom'].astype(str)
    total_predixcan = pd.concat([total_predixcan,predixcan],axis=0)
    
    fusion_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/fusion/analyzed/fusion*')
    fusion = pd.read_csv(fusion_paths[0], sep='\t')
    fusion['ix'] = fusion['gene_id']+gwas_session+'_chr_'+fusion['CHR'].astype(str)
    total_fusion = pd.concat([total_fusion,fusion],axis=0)

total_fusion.index=total_fusion['ix']
total_predixcan.index=total_predixcan['ix']
total_smr.index=total_smr['ix']


total_df = pd.concat([total_fusion['TWAS.P'], total_predixcan['pvalue'], total_smr['p_SMR']],join='inner',axis=1)
total_df.columns = ['fusion','predixcan','smr']
total_df = -np.log10(total_df)
        
condition1 = total_df['fusion'] > 2.40
condition2 = total_df['predixcan'] > 2.27
condition3 = total_df['smr'] > 1.64


total_df['score'] = condition1.astype(int) + condition2.astype(int) + condition3.astype(int)
total_df['ix'] = total_df.index
total_df['gwas_session'] = total_df['ix'].apply(lambda x: 'GCST'+x.split('GCST')[1].split('_')[0])
total_df['gene'] = total_df['ix'].apply(lambda x: x.split('GCST')[0])
total_df['chrom'] = total_df['ix'].apply(lambda x: x.split('_')[-1])

notwas_df = total_df[total_df['score'] < 2]
twas_df = total_df[total_df['score'] >= 2]

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

notwas_effectsize_ls = []
twas_effectsize_ls = []

for gwas in gwasfiledict.keys():
    gwas_df = pd.read_csv(f'/home/users/nus/locuscompare2file/Astle_blood_trait/{gwasfiledict[gwas]}',sep='\t')
    gwas_df.index = "chr" + gwas_df['chr'].astype(str) + "_" + gwas_df['pos'].astype(str) + "_" + gwas_df['ref'] + "_" + gwas_df['alt']
    notwas_df_ = notwas_df[notwas_df['gwas_session'] == gwas]
    for ix, row in tqdm(notwas_df_.iterrows()):
        eqtl = pd.read_csv(f'/home/users/nus/scratch/astle_202502_global_LD_window/preprocessed/eqtl/Whole_Blood/grouped/{row['chrom']}/{row['gene']}.tsv.gz', sep='\t')
        eqtl.index = eqtl['var_id_']
        matchsnp = pd.concat([eqtl['slope'],gwas_df['beta']], join='inner',axis=1)
        matchsnp.columns = ['eqtl_beta','gwas_beta']
        pears = stats.pearsonr(matchsnp['eqtl_beta'],matchsnp['gwas_beta'])
        notwas_effectsize_ls.append(pears.correlation)

    twas_df_ = twas_df[twas_df['gwas_session'] == gwas]
    for ix, row in tqdm(twas_df_.iterrows()):
        eqtl = pd.read_csv(f'/home/users/nus/scratch/astle_202502_global_LD_window/preprocessed/eqtl/Whole_Blood/grouped/{row['chrom']}/{row['gene']}.tsv.gz', sep='\t')
        eqtl.index = eqtl['var_id_']
        matchsnp = pd.concat([eqtl['slope'],gwas_df['beta']], join='inner',axis=1)
        matchsnp.columns = ['eqtl_beta','gwas_beta']
        pears = stats.pearsonr(matchsnp['eqtl_beta'],matchsnp['gwas_beta'])
        twas_effectsize_ls.append(pears.correlation)

        
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Create a DataFrame for easier manipulation
data = pd.DataFrame({
    'Effect Size': notwas_effectsize_ls + twas_effectsize_ls,
    'Category': ['NoTWAS'] * len(notwas_effectsize_ls) + ['TWAS'] * len(twas_effectsize_ls)
})

data.to_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/effectsize_wo_twas.tsv.gz',sep='\t',index=False,compression='gzip')



########### plot ###########
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
from tqdm import tqdm
import seaborn as sns

def box_and_whisker(data, title, ylabel, xticklabels, outname):
    """
    Create a box-and-whisker plot with significance bars without duplicate axis labels.
    """
    fig, ax = plt.subplots(figsize=(3.5,6.5))  # Create a figure with proper size

    # Create boxplot
    bp = ax.boxplot(data, widths=0.6, patch_artist=True)

    # Graph title and labels
    # ax.set_title(title, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    # Remove duplicate x-ticks and only set labels
    ax.set_xticks(range(1, len(xticklabels) + 1))  # Set correct tick positions
    ax.set_xticklabels(xticklabels, fontsize=16)

    # Remove unwanted x-axis major ticks
    ax.tick_params(axis='x', which='major', length=0)

    # Change the colour of the boxes using Seaborn's pastel palette
    colors = sns.color_palette('pastel')
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    # Colour the median lines black
    plt.setp(bp['medians'], color='k')

    # Check for statistical significance
    significant_combinations = []
    ls = list(range(1, len(data) + 1))
    combinations = [(ls[x], ls[x + y]) for y in reversed(ls) for x in range((len(ls) - y))]

    for c in combinations:
        data1 = data[c[0] - 1]
        data2 = data[c[1] - 1]
        U, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        if p < 0.05:
            significant_combinations.append([c, p])

    # Get y-axis range
    bottom, top = ax.get_ylim()
    yrange = top - bottom

    # Add significance bars
    for i, significant_combination in enumerate(significant_combinations):
        # Columns corresponding to the datasets of interest
        x1 = significant_combination[0][0]
        x2 = significant_combination[0][1]
        # What level is this bar among the bars above the plot?
        level = len(significant_combinations) - i
        # Plot the bar
        bar_height = (yrange * 0.08 * level) + top
        bar_tips = bar_height - (yrange * 0.02)
        plt.plot(
            [x1, x1, x2, x2],
            [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k')
        # Significance level
        p = significant_combination[1]
        if p < 0.001:
            sig_symbol = '***'
        elif p < 0.01:
            sig_symbol = '**'
        elif p < 0.05:
            sig_symbol = '*'
        text_height = bar_height + (yrange * 0.01)
        plt.text((x1 + x2) * 0.5, text_height, sig_symbol, ha='center', c='k')

    # Adjust y-axis
    bottom, top = ax.get_ylim()
    yrange = top - bottom
    ax.set_ylim(bottom - 0.02 * yrange, top)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    plt.savefig(f'/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/{outname}.pdf', format='pdf',bbox_inches='tight')
    plt.close()


data = pd.read_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/effectsize_wo_twas.tsv.gz',sep='\t')
data['Effect Size'] = abs(data['Effect Size'])
plotdata = []
plotdata.append(np.array(data[data['Category']=='NoTWAS']['Effect Size']))
plotdata.append(np.array(data[data['Category']=='TWAS']['Effect Size']))
title = 'Comparison of Correlation of Effect Sizes\nBetween Genes w/o TWAS Signal'
ylabel = 'Pearson Correlation of GWAS/eQTL Effect size'
xticklabels = ['TWAS/MR-', 'TWAS/MR+']
outname = 'effectsize_wo_twas'
box_and_whisker(plotdata, title, ylabel, xticklabels, outname)