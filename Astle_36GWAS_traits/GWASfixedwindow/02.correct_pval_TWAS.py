import pandas as pd
import statsmodels.api as sm
import glob
import pandas as pd
import os
import scipy.stats as ss
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

gwas_session_ls = \
    ['GCST004599','GCST004600','GCST004601','GCST004602','GCST004603','GCST004604',
     'GCST004605','GCST004606','GCST004607','GCST004608','GCST004609','GCST004610',
     'GCST004611','GCST004613','GCST004614','GCST004612','GCST004615','GCST004616',
     'GCST004617','GCST004618','GCST004619','GCST004620','GCST004621','GCST004622',
     'GCST004623','GCST004624','GCST004625','GCST004626','GCST004627','GCST004628',
     'GCST004629','GCST004630','GCST004631','GCST004632','GCST004633','GCST004634']


total_fusion = pd.DataFrame()
total_predixcan = pd.DataFrame()
total_smr = pd.DataFrame()


for gwas_session in gwas_session_ls:
    fusion_paths = glob.glob(f'~/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/fusion/analyzed/fusion*')
    fusion = pd.read_csv(fusion_paths[0], sep='\t')
    fusion['ix'] = fusion['gene_id']+gwas_session
    total_fusion = pd.concat([total_fusion,fusion],axis=0)

total_fusion = total_fusion.sort_values('TWAS.P') # 567864
total_fusion = total_fusion.drop_duplicates('ix')
total_fusion.iloc[:,2:].to_csv('~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/fusion_qvalue.tsv.gz',sep='\t',compression='gzip',index=False)


for gwas_session in gwas_session_ls:
    smr_paths = glob.glob(f'~/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/smr/analyzed/smr*')
    smr = pd.read_csv(smr_paths[0], sep='\t')
    smr['ix'] = smr['gene_id']+gwas_session
    total_smr = pd.concat([total_smr,smr],axis=0)
total_smr = total_smr.sort_values('p_SMR') 

total_smr = total_smr.drop_duplicates('ix')
total_smr.to_csv('~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/smr_qvalue.tsv.gz',sep='\t',compression='gzip',index=False)



for gwas_session in gwas_session_ls:
    predixcan_paths = glob.glob(f'~/scratch/astle_202501_global_LD_window/processed/default/{gwas_session}_fixed_GWAS_Loci_window/Whole_Blood/EUR/eqtl/predixcan/analyzed/predixcan*')
    predixcan = pd.read_csv(predixcan_paths[0], sep='\t')
    predixcan = predixcan[['gene', 'gene_name', 'zscore', 'effect_size', 'pvalue','chrom', 'gene_id']]
    predixcan = predixcan.dropna()
    predixcan['ix'] = predixcan['gene_id']+gwas_session
    total_predixcan = pd.concat([total_predixcan,predixcan],axis=0)

total_predixcan = total_predixcan.sort_values('pvalue') 
total_predixcan = total_predixcan.drop_duplicates('ix')
total_predixcan.to_csv('~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/predixcan_qvalue.tsv.gz',sep='\t',compression='gzip',index=False)

