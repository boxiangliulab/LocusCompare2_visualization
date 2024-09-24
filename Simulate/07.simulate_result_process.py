


# 他们的predixcan model和新的twas model
# 设置threshold pvalue
#   eqtl: 1
#   gwas: 1
import os
import pandas as pd
from multiprocessing import Pool
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from tqdm import tqdm
import random
import scipy.stats as ss
from sklearn.metrics import average_precision_score, precision_recall_curve
import seaborn as sns
from scipy.stats import pearsonr


gene_name_ls = pd.read_csv('/home/liulab/datasource/gene_name.tsv')

tools = ['coloc','ecaviar','fastenloc','predixcan','smr','twas']
tools_metrics = {'coloc': 'overall_H4', 
                 'ecaviar': 'clpp',
                 'predixcan': 'pvalue',
                 'twas': 'TWAS.P',
                 'smr': 'p_SMR',
                 'fastenloc': 'LCP'}
coloc = pd.read_csv('/home/liulab/datasource/coloc_result.tsv',sep='\t')
coloc.index = coloc['gene_id']



twas_ls = pd.read_csv('/home/liulab/datasource/locuscompare2_twas20240507/twas_remain_ls',header=None)
twas_ls = list(twas_ls[0])

## new
## twas 
newtwas = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_twas20240510/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/twas/analyzed/')

    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 0:
            tmp_df = pd.read_csv(os.path.join(tmp_path, file_name[0]),sep='\t')
            if gene_name in list(tmp_df['gene_id']):
                tmp_df = tmp_df[tmp_df['gene_id'] == gene_name]
                newtwas = pd.concat([newtwas, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)

newtwas.index = newtwas['gene_id']
newtwas = pd.concat([newtwas,coloc['IS_POSITIVE']],axis=1)

## predixcan 
newpredixcan = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_20240503_predixcan/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/predixcan/analyzed/')
    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 0:
            tmp_df = pd.read_csv(os.path.join(tmp_path, file_name[0]),sep='\t')
            if gene_name in list(tmp_df['gene_id']):
                tmp_df = tmp_df[tmp_df['gene_id'] == gene_name]
                newpredixcan = pd.concat([newpredixcan, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)

newpredixcan.index = newpredixcan['gene_id'] # 716 genes
newpredixcan = pd.concat([newpredixcan,coloc['IS_POSITIVE']],axis=1)

## ecaviar 
newecaviar = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_fastenlocdebug20240504/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/ecaviar/analyzed/')
    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 0:
            tmp_df = pd.read_csv(os.path.join(tmp_path, file_name[0]),sep='\t')
            if gene_name in list(tmp_df['gene_id']):
                tmp_df = tmp_df[tmp_df['gene_id'] == gene_name]
                newecaviar = pd.concat([newecaviar, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)

newecaviar.index = newecaviar['gene_id']
newecaviar = pd.concat([newecaviar,coloc['IS_POSITIVE']],axis=1)

## coloc 
newcoloc = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_fastenlocdebug20240504/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/coloc/analyzed/')
    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 0:
            tmp_df = pd.read_csv(os.path.join(tmp_path, file_name[0]),sep='\t')
            if gene_name in list(tmp_df['gene_id']):
                tmp_df = tmp_df[tmp_df['gene_id'] == gene_name]
                newcoloc = pd.concat([newcoloc, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)

newcoloc.index = newcoloc['gene_id']
newcoloc = pd.concat([newcoloc,coloc['IS_POSITIVE']],axis=1)

## smr
newsmr = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_fastenlocdebug20240428/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/smr/analyzed/')
    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 0:
            tmp_df = pd.read_csv(os.path.join(tmp_path, file_name[0]),sep='\t')
            if gene_name in list(tmp_df['gene_id']):
                tmp_df = tmp_df[tmp_df['gene_id'] == gene_name]
                newsmr = pd.concat([newsmr, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)

newsmr.index = newsmr['gene_id']
newsmr = pd.concat([newsmr,coloc['IS_POSITIVE']],axis=1)

## fastenloc 
newfastenloc = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_fastenlocdebug20240504/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/fastenloc/analyzed/')
    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 3:
            for results in file_name:
                if results[:9] == 'fastenloc':
                    break
            tmp_df = pd.read_csv(os.path.join(tmp_path, results),sep='\t')
            if gene_name in list(tmp_df['gene_id']):
                tmp_df = tmp_df[tmp_df['gene_id'] == gene_name]
                newfastenloc = pd.concat([newfastenloc, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)


newfastenloc_GLCP = pd.DataFrame()
for gene_name in gene_name_ls['gene']:
    tmp_path = os.path.join('/home/liulab/datasource/locuscompare2_fastenlocdebug20240504/', gene_name,'processed/default/', gene_name, 'test_tissue/EUR/fastenloc/analyzed/')
    if not os.path.exists(tmp_path):
        continue
    else:
        file_name = os.listdir(tmp_path)
        if len(file_name) > 3:
            for results in file_name:
                if results[-8:] == 'gene.out':
                    break
            tmp_df = pd.read_csv(os.path.join(tmp_path, results),sep='\s+')
            if gene_name in list(tmp_df['Gene']):
                tmp_df = tmp_df[tmp_df['Gene'] == gene_name]
                newfastenloc_GLCP = pd.concat([newfastenloc_GLCP, pd.DataFrame(tmp_df.iloc[0,:]).T], axis=0)

newfastenloc_GLCP.index = newfastenloc_GLCP['Gene'] # 1137
newfastenloc_GLCP = pd.concat([newfastenloc_GLCP,coloc['IS_POSITIVE']],axis=1)
newfastenloc_GLCP['Gene'] = newfastenloc_GLCP.index
newfastenloc_GLCP= newfastenloc_GLCP.fillna(0)

newfastenloc.index = newfastenloc['gene_id'] # 1137
newfastenloc = pd.concat([newfastenloc,coloc['IS_POSITIVE']],axis=1)

new_df = pd.concat([newcoloc['overall_H4'], newsmr[['b_SMR','se_SMR','p_SMR']], newfastenloc['LCP'], newecaviar['clpp'], newtwas[['TWAS.Z','TWAS.P']],newpredixcan[['zscore','pvalue']]],join='inner', axis=1)

new_df.columns = ['coloc','SMR.b','SMR.se','SMR.p','fastenloc','ecaviar','TWAS.Z', 'TWAS.P', 'PrediXcan.Z', 'PrediXcan.P']
new_df = pd.concat([new_df,coloc['IS_POSITIVE']],join='inner',axis=1)
new_df_dropna = new_df.dropna() # 637 gene
new_df_dropna.to_csv('LocusCompare2_visualization/Simulate//new_sim_gwaseqtlresult_dropna_20240511.tsv',sep='\t')

new_df['PrediXcan.P'] = new_df['PrediXcan.P'].fillna(1)
new_df['SMR.p'] = new_df['SMR.p'].fillna(1)
new_df['TWAS.P'] = new_df['TWAS.P'].fillna(1)
new_df['coloc'] = new_df['coloc'].fillna(0)
new_df['fastenloc'] = new_df['fastenloc'].fillna(0)
new_df['ecaviar'] = new_df['ecaviar'].fillna(0)
new_df['PrediXcan.Z'] = new_df['PrediXcan.Z'].fillna(0)
new_df['TWAS.Z'] = new_df['TWAS.Z'].fillna(0)
new_df['SMR.z'] = new_df['SMR.b']/new_df['SMR.se']
new_df.to_csv('LocusCompare2_visualization/Simulate//new_sim_gwaseqtlresult_20240511.tsv',sep='\t')

## compare twas
compare_twas_old_new = pd.concat([old_df[['TWAS.Z','TWAS.P']],newtwas[['TWAS.Z','TWAS.P']]],join='inner', axis=1)
compare_twas_old_new.columns = ['oldz','oldp','newz','newp']
compare_twas_old_new = pd.concat([compare_twas_old_new,coloc['IS_POSITIVE']],axis=1)
compare_twas_old_new = compare_twas_old_new.dropna()

pearsonr(compare_twas_old_new['oldz'],compare_twas_old_new['newz'])
# Out[34]: PearsonRResult(statistic=0.9878833025928685, pvalue=0.0)
pearsonr(compare_twas_old_new['oldp'],compare_twas_old_new['newp'])
# Out[35]: PearsonRResult(statistic=0.9258481353688646, pvalue=0.0)
group = list(compare_twas_old_new['IS_POSITIVE'])


newtwas['TWAS.P'] = newtwas['TWAS.P'].fillna(1)
newtwas['TWAS.Z'] = newtwas['TWAS.Z'].fillna(0)
newtwas.to_csv('LocusCompare2_visualization/Simulate/newtwas2_merged.tsv',sep='\t',index=False)

