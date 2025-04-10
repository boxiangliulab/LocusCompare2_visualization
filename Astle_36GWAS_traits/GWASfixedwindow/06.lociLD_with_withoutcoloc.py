import pandas as pd
import glob
import os

gwas_session_ls = \
    ['GCST004599','GCST004600','GCST004601','GCST004602','GCST004603','GCST004604',
     'GCST004605','GCST004606','GCST004607','GCST004608','GCST004609','GCST004610',
     'GCST004611','GCST004613','GCST004614','GCST004612','GCST004615','GCST004616',
     'GCST004617','GCST004618','GCST004619','GCST004620','GCST004621','GCST004622',
     'GCST004623','GCST004624','GCST004625','GCST004626','GCST004627','GCST004628',
     'GCST004629','GCST004630','GCST004631','GCST004632','GCST004633','GCST004634']

total_ecaviar_df = pd.DataFrame()
ix_ls = []
qtl_causal_ls = []
gwas_causal_ls = []
distance_gwas_qtl_ls = []

for gwas_session in gwas_session_ls:
    total_gene = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/ecaviar/candidate/*')
    for each_gene in total_gene:
        gene_name = os.path.basename(each_gene)
        total_loci_eachgene = glob.glob(f"{each_gene}/*")
        for locus in total_loci_eachgene:
            locus_name = os.path.basename(locus)
            if glob.glob(f"{locus}/gwas_*snp") == []:
                continue
            if glob.glob(f"{locus}/qtl_*snp") == []:
                continue
            gwas_file = glob.glob(f"{locus}/gwas_*snp")[0]
            qtl_file = glob.glob(f"{locus}/qtl_*snp")[0]
            gwas_df = pd.read_csv(gwas_file,sep=' ')
            qtl_df = pd.read_csv(qtl_file,sep=' ')
            gwas_causal_line = gwas_df.sort_values('prob',ascending=False).iloc[0,:]
            qtl_causal_line = qtl_df.sort_values('prob',ascending=False).iloc[0,:]
            gwas_causal = 'chr'+str(gwas_causal_line['chromosome'])+\
                '_'+str(gwas_causal_line['position'])+\
                '_'+gwas_causal_line['allele2']+\
                '_'+gwas_causal_line['allele1']
            qtl_causal = 'chr'+str(qtl_causal_line['chromosome'])+\
                '_'+str(qtl_causal_line['position'])+\
                '_'+qtl_causal_line['allele2']+\
                '_'+qtl_causal_line['allele1']
            gwas_causal_pos = gwas_df.sort_values('prob',ascending=False).iloc[0,:]['position']
            qtl_causal_pos = qtl_df.sort_values('prob',ascending=False).iloc[0,:]['position']
            distance_gwas_qtl_causals = abs(gwas_causal_pos-qtl_causal_pos)
            ix_ls.append(f"{gwas_session}_{gene_name}_{locus_name}")
            qtl_causal_ls.append(qtl_causal)
            gwas_causal_ls.append(gwas_causal)
            distance_gwas_qtl_ls.append(distance_gwas_qtl_causals)
        
        
    
total_ecaviar_df['ix'] = ix_ls
total_ecaviar_df['qtl_causal'] = qtl_causal_ls
total_ecaviar_df['gwas_causal'] = gwas_causal_ls
total_ecaviar_df['distance_gwas_qtl'] = distance_gwas_qtl_ls
total_ecaviar_df.to_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/ecaviar_gwas_qtl_causal.tsv.gz',sep='\t',compression='gzip',index=False)

reference = pd.read_csv('/home/users/nus/scratch/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz',sep='\t')
total_ecaviar_df = pd.read_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/ecaviar_gwas_qtl_causal.tsv.gz',sep='\t')
reference['variant_id'] = reference['variant_id'].apply(lambda x: '_'.join(x.split('_')[:4]))

# 合并 A 和 B，左连接，依据 qtl_causal 和 variant_id
df_merged = total_ecaviar_df.merge(reference[['variant_id','rs_id_dbSNP151_GRCh38p7']], how='left', left_on='qtl_causal', right_on='variant_id')

# 只保留你感兴趣的列
df_merged = df_merged[['ix', 'qtl_causal', 'gwas_causal', 'distance_gwas_qtl',
       'rs_id_dbSNP151_GRCh38p7']]
df_merged.columns = ['ix', 'qtl_causal', 'gwas_causal', 'distance_gwas_qtl',
       'qtl_rsid']
df_merged = df_merged.merge(reference[['variant_id','rs_id_dbSNP151_GRCh38p7']], how='left', left_on='gwas_causal', right_on='variant_id')

df_merged = df_merged[['ix', 'qtl_causal', 'gwas_causal', 'distance_gwas_qtl', 'qtl_rsid',
       'rs_id_dbSNP151_GRCh38p7']]

df_merged.columns = ['ix', 'qtl_causal', 'gwas_causal', 'distance_gwas_qtl', 'qtl_rsid',
       'gwas_rsid']
df_merged['chrom'] = df_merged['gwas_causal'].apply(lambda x: x.split('_')[0].strip('chr'))
df_merged.to_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/ecaviar_gwas_qtl_causal.tsv.gz',sep='\t',compression='gzip',index=False)

print(df_merged)
df_merged.sort_values('chrom',inplace=True)

vcf_home = '/home/users/nus/locuscompare2file/vcf_hg38/EUR/'
for chrom in list(set(df_merged['chrom'])):
    vcf_path = os.path.join(vcf_home, f"chr{chrom}.vcf.gz")
    tmp = df_merged[df_merged['chrom'] == chrom]
    chr_snp_ls = list(set(list(tmp['qtl_rsid'])+list(tmp['gwas_rsid'])))
    chr_snp = pd.DataFrame(chr_snp_ls)
    snp_list_path = f"/home/users/nus/scratch/ecaviar_gwas_qtl_causal_chr{chrom}_snplist.txt"
    chr_snp.to_csv(snp_list_path,index=False,header=False)
    out_prefix = f"/home/users/nus/scratch/ecaviar_gwas_qtl_causal_chr{chrom}_snpld"
    os.system(f"plink --vcf {vcf_path} \
          --ld-snp-list {snp_list_path} \
          --r2 \
          --ld-window 99999 \
          --ld-window-kb 10000 \
          --ld-window-r2 0 \
          --out {out_prefix}")



import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
from tqdm import tqdm
import seaborn as sns


total_ecaviar_df  = pd.read_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/ecaviar_gwas_qtl_causal.tsv.gz',sep='\t')
total_ecaviar_df['rsidpair'] = total_ecaviar_df['qtl_rsid']+total_ecaviar_df['gwas_rsid']
totalnew_ecaviar_df = pd.DataFrame()
for chrom in list(set(total_ecaviar_df['chrom'])):
    tmp = total_ecaviar_df[total_ecaviar_df['chrom'] == chrom]
    out_path = f"/home/users/nus/scratch/ecaviar_gwas_qtl_causal_chr{chrom}_snpld.ld"
    ld = pd.read_csv(out_path,sep='\s+')
    rsidpair1 = list(ld['SNP_A']+ld['SNP_B'])
    rsidpair2 = list(ld['SNP_B']+ld['SNP_A'])
    ld = pd.concat([ld,ld])
    ld['rsidpair'] = rsidpair1 + rsidpair2
    ld.index = ld['rsidpair']
    overlappair = list(set(ld['rsidpair']) & set(tmp['rsidpair']))
    ld = ld.loc[overlappair]
    ld.drop_duplicates('rsidpair',inplace=True)
    ld_dict = ld['R2'].to_dict()
    tmp.index = tmp['rsidpair']
    tmp['R2'] = tmp['rsidpair'].map(ld_dict)
    tmp = tmp[['chrom','qtl_causal', 'gwas_causal', 'R2','distance_gwas_qtl', 'qtl_rsid',
       'gwas_rsid', 'ix']]
    totalnew_ecaviar_df = pd.concat([totalnew_ecaviar_df, tmp])
    del tmp
    del ld
    del rsidpair1
    del rsidpair2
    del overlappair
    
    
totalnew_ecaviar_df.to_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/ecaviar_gwas_qtl_causal_r2.tsv.gz',sep='\t',compression='gzip',index=False)


coloc_df = pd.DataFrame()
fastenloc_df = pd.DataFrame()
ecaviar_df = pd.DataFrame()

gwas_session_ls = \
    ['GCST004599','GCST004600','GCST004601','GCST004602','GCST004603','GCST004604',
     'GCST004605','GCST004606','GCST004607','GCST004608','GCST004609','GCST004610',
     'GCST004611','GCST004613','GCST004614','GCST004612','GCST004615','GCST004616',
     'GCST004617','GCST004618','GCST004619','GCST004620','GCST004621','GCST004622',
     'GCST004623','GCST004624','GCST004625','GCST004626','GCST004627','GCST004628',
     'GCST004629','GCST004630','GCST004631','GCST004632','GCST004633','GCST004634']

    
for gwas_session in gwas_session_ls:
    coloc_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/coloc/analyzed/coloc*')
    coloc = pd.read_csv(coloc_paths[0], sep='\t')
    coloc = coloc.sort_values('overall_H4',ascending=False)
    coloc = coloc.drop_duplicates(['gene_id', 'gwas_lead_snp'],keep='first')
    coloc.index = gwas_session+'_'+ coloc['gene_id']+'_'+ coloc['gwas_lead_snp']
    
    fastenloc_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/fastenloc/analyzed/fastenloc_outp*')
    fastenloc = pd.read_csv(fastenloc_paths[0], sep='\t')
    fastenloc = fastenloc.sort_values('GRCP',ascending=False)
    fastenloc = fastenloc.drop_duplicates(['gene_id', 'lead_variant'],keep='first')
    fastenloc.index = gwas_session+'_'+fastenloc['gene_id']+'_'+fastenloc['lead_variant']
    
    ecaviar_paths = glob.glob(f'/home/users/nus/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/ecaviar/analyzed/ecaviar*')
    ecaviar = pd.read_csv(ecaviar_paths[0], sep='\t')
    ecaviar = ecaviar.sort_values('clpp',ascending=False)
    ecaviar = ecaviar.drop_duplicates(['gene_id', 'lead_variant'],keep='first')
    ecaviar.index = gwas_session+'_'+ecaviar['gene_id']+'_'+ecaviar['lead_variant']

    
    coloc_df = pd.concat([coloc_df, coloc[['overall_H4']]], axis=0)
    fastenloc_df = pd.concat([fastenloc_df, fastenloc['GRCP']], axis=0)
    ecaviar_df = pd.concat([ecaviar_df, ecaviar['clpp']], axis=0)
    
    
total_df = pd.concat([coloc_df['overall_H4'], fastenloc_df['GRCP'], ecaviar_df['clpp']],join='inner',axis=1)
        
condition1 = total_df['overall_H4'] > 0.75
condition2 = total_df['GRCP'] > 0.5
condition3 = total_df['clpp'] > 0.01

# 将布尔值转换为整数 (True=1, False=0)，然后按行相加
total_df['score'] = condition1.astype(int) + condition2.astype(int) + condition3.astype(int)
totalnew_ecaviar_df.index = totalnew_ecaviar_df['ix']
total_df = pd.concat([total_df, totalnew_ecaviar_df[['distance_gwas_qtl','R2']]],join='inner',axis=1)
nocoloc_df_sum = total_df[total_df['score'] < 2]
coloc_df_sum = total_df[total_df['score'] >= 2]


nocoloc_distance_ls = list(nocoloc_df_sum['distance_gwas_qtl'])
coloc_distance_ls = list(coloc_df_sum['distance_gwas_qtl'])
nocoloc_LD_ls = list(nocoloc_df_sum['R2'])
coloc_LD_ls = list(coloc_df_sum['R2'])

data = pd.DataFrame({
    'Distance_G/E': nocoloc_distance_ls + coloc_distance_ls,
    'LD': nocoloc_LD_ls + coloc_LD_ls,
    'Category': ['NoColocalization'] * len(nocoloc_df_sum) + ['Colocalization'] * len(coloc_df_sum)
})

data.to_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/distance_ld_wo_coloc_withecaviarresults.tsv.gz',sep='\t',index=False,compression='gzip')
data = pd.read_csv('/home/users/nus/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/distance_ld_wo_coloc_withecaviarresults.tsv.gz',sep='\t')


############ plot ############
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
    #plt.show()
    plt.savefig(f'/Users/phoebel/github/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/{outname}.pdf', format='pdf',bbox_inches='tight')
    plt.close()

        
        
data = pd.read_csv('~/github/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/distance_ld_wo_coloc_withecaviarresults.tsv.gz',sep='\t')
data = data.dropna()

plotdata = []
plotdata.append(np.array(data[data['Category']=='NoColocalization']['Distance_G/E']))
plotdata.append(np.array(data[data['Category']=='Colocalization']['Distance_G/E']))
title = ''
ylabel = 'Distance between GWAS/eQTL causal SNPs'
xticklabels = ['Coloc−', 'Coloc+']
outname = 'distance_wo_coloc_withecaviarresults'
box_and_whisker(plotdata, title, ylabel, xticklabels, outname)


plotdata = []
plotdata.append(np.array(data[data['Category']=='NoColocalization']['LD']))
plotdata.append(np.array(data[data['Category']=='Colocalization']['LD']))
title = ''
ylabel = f"LD($\\it{{r^2}}$) between GWAS/eQTL causal SNPs"
xticklabels = ['Coloc−', 'Coloc+']
outname = 'ld_wo_coloc_withecaviarresults'
box_and_whisker(plotdata, title, ylabel, xticklabels, outname)




