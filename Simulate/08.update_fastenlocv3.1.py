import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from scipy import stats
import scipy.stats as ss
import seaborn as sns
from sklearn.metrics import average_precision_score, precision_recall_curve
import os

gwas = pd.read_csv('sim_gwas_20240430/out.tsv.gz',sep='\t')
eqtl = pd.read_csv('sim_eqtl_20240430/out.tsv.gz',sep='\t')
gwas.index = gwas['phe_id']+gwas['var_id']
eqtl.index = eqtl['phe_id']+eqtl['var_id']
df = pd.concat([gwas[['phe_id','var_id','beta','se']], eqtl[['beta','se']]], join='inner',axis=1)
df.to_csv('~/scratch/simulate_pre_fastenloc.tsv.gz',sep='\t',header=False, index=False,compression='gzip')

os.system("fastenloc -sum ~/scratch/simulate_pre_fastenloc.tsv.gz -prefix ~/scratch/simulate_fastenlocv3.1")


df = pd.read_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/new_sim_gwaseqtlresult_20240511_qvalue_intact.tsv',sep='\t')
df = df.iloc[:,:-11]
df.index = df['gene_id']
df = df.drop('fastenloc',axis=1)

newfastenloc = pd.read_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/simulate_fastenlocv3.1.enloc.gene.out',sep='\t')
newfastenloc.index = newfastenloc['Gene']
df = pd.concat([df,newfastenloc['GRCP']],join='inner',axis=1)
df.columns = ['gene_id', 'coloc', 'SMR.b', 'SMR.se', 'SMR.p', 'ecaviar',
       'TWAS.Z', 'TWAS.P', 'PrediXcan.Z', 'PrediXcan.P', 'SMR.z',
       'IS_POSITIVE', 'smrqvalue', 'twasqvalue', 'predixcanqvalue', 'fastenloc']
df.to_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/new_sim_gwaseqtlresult_20240922_grcp_qvalue.tsv',sep='\t',index=False)

coloc = pd.read_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newcoloc_merged.tsv',sep='\t')
coloc.index = coloc['gene_id']
newfastenloc = pd.concat([newfastenloc[['Gene','GRCP','GLCP']], coloc['IS_POSITIVE']],join='inner',axis=1)
newfastenloc.columns = ['gene_id', 'GRCP', 'GLCP', 'IS_POSITIVE']
newfastenloc.to_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/newfastenlocv3.1_GRCP_merged.tsv',sep='\t',index=False)
