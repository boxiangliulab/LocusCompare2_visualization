import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from scipy import stats
import scipy.stats as ss
import seaborn as sns
from sklearn.metrics import average_precision_score, precision_recall_curve

df = pd.read_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/new_sim_gwaseqtlresult_20240922_grcp_qvalue.tsv',sep='\t')
df.index = df['gene_id']


#coloc = pd.read_csv('/Users/phoebel/tmp_doc/coloc_result.tsv',sep='\t')
#coloc.index = coloc['gene_id']

#df = pd.concat([df,coloc['IS_POSITIVE']],join='inner',axis=1)

def roc_plot(df, column_name, color = None, linewidth=2.0):
    df.sort_values(column_name, inplace=True, ascending=False)
    fpr, tpr, thresholds = roc_curve(df['IS_POSITIVE'], df[column_name])
    AUC = auc(fpr, tpr)
    print(f"{column_name} AUC = {round(AUC, 4)}")
    plt.plot(fpr, tpr, linewidth=linewidth, color=color, label=f"{column_name} AUC = {round(AUC, 4)}")
    


df.index = df['gene_id']

df['smrqvalue'] = 1-df['smrqvalue']
df['predixcanqvalue'] = 1-df['predixcanqvalue']
df['twasqvalue'] = 1-df['twasqvalue']

def ranking(array):
    array = ss.rankdata(array)
    array = [int(x) for x in array]
    array = [x/1197 for x in array]
    return np.array(array)

def arithemtic_mean_func(x):
    x=list(x)
    while 99 in x:
        x.remove(99)
    x=np.array(x)
    return np.average(x)

  
clean_df = df[['coloc','fastenloc','ecaviar','smrqvalue','predixcanqvalue','twasqvalue','IS_POSITIVE']]
clean_df.columns = ['coloc','fastenloc','ecaviar','smr','predixcan','fusion','IS_POSITIVE']

## prob value
mean_top2 = (clean_df.iloc[:,:6].apply(lambda x: arithemtic_mean_func(sorted(x)[-2:]),axis=1))
ensemble = (clean_df.iloc[:,:6].apply(lambda x: np.average(sorted(x), weights=[0,0,0,0,0,1]),axis=1))
meanofmaxprobin2class = (clean_df.iloc[:,:6].apply(lambda x: (np.max([x['fusion'],x['smr'],x['predixcan']]) + np.max([x['ecaviar'], x['coloc'], x['fastenloc']]))/2, axis=1))



clean_df['mean_top2'] = mean_top2
clean_df['ensemble'] = ensemble
clean_df['meanofmaxprobin2class'] = meanofmaxprobin2class

## all metrics
plt.figure(figsize=(6,6))
plt.title('ROC Curve',fontsize=15)

roc_plot(clean_df, 'coloc') # coloc AUC = 0.7804
roc_plot(clean_df, 'predixcan') # predixcan AUC = 0.7106
roc_plot(clean_df, 'smr') # smr AUC = 0.8304
roc_plot(clean_df, 'fastenloc') # fastenloc AUC = 0.7658
roc_plot(clean_df, 'ecaviar') # ecaviar AUC = 0.8264
roc_plot(clean_df, 'fusion') # twas AUC = 0.7716

plt.legend(loc='lower right', prop = { "size": 10})
plt.plot([(0, 0), (1, 1)], 'k--', alpha=0.1)
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 01.01])
plt.yticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.xticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.ylabel('True Positive Rate', fontsize=15)
plt.xlabel('False Positive Rate', fontsize=15)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_roc_6tools_grcp.pdf',format='pdf',dpi=300)
plt.close()

plt.figure(figsize=(6,6))
plt.title('ROC Curve',fontsize=15)
#sns.set_style("whitegrid")
roc_plot(clean_df, 'coloc', color = '#1f77b4') # coloc AUC = 0.7804
roc_plot(clean_df, 'predixcan', color = '#ff7f0e') # predixcan AUC = 0.7106
roc_plot(clean_df, 'smr', color = '#2ca02c') # smr AUC = 0.8304
roc_plot(clean_df, 'fastenloc', color = '#d62728') # fastenloc AUC = 0.7658
roc_plot(clean_df, 'ecaviar', color = '#9467bd') # ecaviar AUC = 0.8264
roc_plot(clean_df, 'fusion', color = '#8c564b') # twas AUC = 0.7716
roc_plot(clean_df, 'meanofmaxprobin2class', color = '#bcff94', linewidth=3) # meanofmaxprobin2class AUC = 0.837
roc_plot(clean_df, 'ensemble', color = '#1eb2c4', linewidth=3) # ensemble AUC = 0.8388
plt.legend(loc='lower right', prop = { "size": 10})
plt.plot([(0, 0), (1, 1)], 'k--', alpha=0.1)
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 01.01])
plt.yticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.xticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.ylabel('True Positive Rate', fontsize=15)
plt.xlabel('False Positive Rate', fontsize=15)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_roc_allmetrics_grcp.pdf',format='pdf',dpi=300)
plt.close()



## ensemble
plt.figure(figsize=(6,6))
plt.title('ROC Curve',fontsize=15)
sns.set_style("whitegrid")
roc_plot(clean_df, 'coloc', color = '#1f77b4', linewidth=2.0)
roc_plot(clean_df, 'predixcan', color = '#ff7f0e', linewidth=2.0)
roc_plot(clean_df, 'smr', color = '#2ca02c', linewidth=2.0)
roc_plot(clean_df, 'fastenloc', color = '#d62728', linewidth=2.0)
roc_plot(clean_df, 'ecaviar',color = '#9467bd', linewidth=2.0)
roc_plot(clean_df, 'fusion', color = '#8c564b', linewidth=2.0)
roc_plot(clean_df, 'ensemble', color = '#1eb2c4', linewidth=3.0)
plt.legend(loc='lower right', prop = { "size": 10})
plt.plot([(0, 0), (1, 1)], 'k--', alpha=0.1)
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 01.01])
plt.yticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.xticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.ylabel('True Positive Rate', fontsize=15)
plt.xlabel('False Positive Rate', fontsize=15)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_roc_ensemble_grcp.png',format='png',dpi=300)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_roc_ensemble_grcp.pdf',format='pdf',dpi=300)
plt.close()



plt.figure(figsize=(4,9))
roc_plot(clean_df, 'coloc', color = '#1f77b4', linewidth=3.0)
roc_plot(clean_df, 'predixcan', color = '#ff7f0e', linewidth=3.0)
roc_plot(clean_df, 'smr', color = '#2ca02c', linewidth=3.0)
roc_plot(clean_df, 'fastenloc', color = '#d62728', linewidth=3.0)
roc_plot(clean_df, 'ecaviar',color = '#9467bd', linewidth=3.0)
roc_plot(clean_df, 'fusion', color = '#8c564b', linewidth=3.0)
roc_plot(clean_df, 'ensemble', color = '#1eb2c4', linewidth=4.0)
# plt.legend(loc='lower right', prop = { "size": 7})
# plt.plot([(0, 0), (1, 1)], 'k--', alpha=0.1)
plt.xlim([-0.01, 0.06])
plt.ylim([-0.01, 0.8])
# plt.ylabel('True Positive Rate', fontsize=15)
# plt.xlabel('False Positive Rate', fontsize=15)
plt.yticks(np.array([0,0.2,0.4,0.6,0.8]), fontsize=19)
plt.xticks(np.array([0,0.01,0.02,0.03,0.04,0.05]), fontsize=19)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_roc_ensemble_x_0.05_grcp.png',format='png',dpi=300)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_roc_ensemble_x_0.05_grcp.pdf',format='pdf',dpi=300)
plt.close()


## Heatmap

plt.figure()
plt.title("prob qvalue Heatmap with na")
plt.rcParams['font.size'] = '20'
sns.set(font_scale=1.5)
sns.clustermap(clean_df.iloc[:,:6].corr(), annot=True, cmap="crest",annot_kws={"size": 20})
# plt.yticks(fontsize=25)
# plt.xticks(fontsize=25)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_heatmap_withna_grcp.png',format='png',dpi=300)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_heatmap_withna_grcp.pdf',format='pdf',dpi=300)
plt.close()

  
clean_df2 = df[['coloc','fastenloc','ecaviar','SMR.p','PrediXcan.P','TWAS.P','IS_POSITIVE']]
clean_df2.columns = ['coloc','fastenloc','ecaviar','smr','predixcan','fusion','IS_POSITIVE']
clean_df2['smr'] = [-np.log10(i) for i in list(clean_df2['smr'])]
clean_df2['predixcan'] = [-np.log10(i) for i in list(clean_df2['predixcan'])]
clean_df2['fusion'] = [-np.log10(i) for i in list(clean_df2['fusion'])]

# clean_df2['predixcan'] = -np.log10(clean_df2['predixcan'])
# clean_df2['twas'] = -np.log10(clean_df2['twas'])

plt.figure()
plt.title("prob -log10(pval) Heatmap")
plt.rcParams['font.size'] = '20'
sns.set(font_scale=1.7)
sns.clustermap(clean_df2.iloc[:,:6].corr(), annot=True, cmap="crest",annot_kws={"size": 20})
# plt.yticks(fontsize=25)
# plt.xticks(fontsize=25)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_heatmap_withna_log10p_grcp.png',format='png',dpi=300)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/prob_qvalue_heatmap_withna_log10p_grcp.pdf',format='pdf',dpi=300)
plt.close()
