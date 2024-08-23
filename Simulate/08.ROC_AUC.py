import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from scipy import stats
import scipy.stats as ss
import seaborn as sns
from sklearn.metrics import average_precision_score, precision_recall_curve
df = pd.read_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/new_sim_gwaseqtlresult_20240511_qvalue_intact.tsv',sep='\t')
df = df.iloc[:,:-11]
df.index = df['gene_id']

newfastenloc_GLCP = pd.read_csv('github/locuscompare2_private/new_sim_result_20240511/newfastenloc_GLCP_merged.tsv',sep='\t')
newfastenloc_GLCP.index = newfastenloc_GLCP['Gene']
df['fastenloc'] = newfastenloc_GLCP['GLCP']
df.to_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/new_sim_gwaseqtlresult_20240511_glcp_qvalue.tsv',sep='\t',index=False)

df = pd.read_csv('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/new_sim_gwaseqtlresult_20240511_glcp_qvalue_intact.tsv',sep='\t')
df.index = df['gene_id']



def roc_plot(df, column_name, linewidth=2.0):
    df.sort_values(column_name, inplace=True, ascending=False)
    fpr, tpr, thresholds = roc_curve(df['IS_POSITIVE'], df[column_name])
    AUC = auc(fpr, tpr)
    # m = pd.DataFrame()
    # m['fpr'] = fpr
    # m['tpr'] = tpr
    # AUC = auc(m[m['fpr']<=0.05]['fpr'],m[m['fpr']<=0.05]['tpr'])
    print(f"{column_name} AUC = {round(AUC, 4)}")
    plt.plot(fpr, tpr, linewidth=linewidth, label=f"{column_name} AUC = {round(AUC, 4)}")
    

df.index = df['gene_id']

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
roc_plot(clean_df, 'fastenloc') # fastenloc AUC = 0.0.7709
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
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/prob_qvalue_roc_6tools_glcp.png',format='png',dpi=300)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/prob_qvalue_roc_6tools_glcp.pdf',format='pdf',dpi=300)
plt.close()

plt.figure(figsize=(6,6))
plt.title('ROC Curve',fontsize=15)

## ensemble
plt.figure(figsize=(6,6))
plt.title('ROC Curve',fontsize=15)
roc_plot(clean_df, 'coloc') # coloc AUC = 0.7804
roc_plot(clean_df, 'predixcan') # predixcan AUC = 0.7106
roc_plot(clean_df, 'smr') # smr AUC = 0.8304
roc_plot(clean_df, 'fastenloc') # fastenloc AUC = 0.7709
roc_plot(clean_df, 'ecaviar') # ecaviar AUC = 0.8264
roc_plot(clean_df, 'fusion') # twas AUC = 0.7716
roc_plot(clean_df, 'ensemble',linewidth=3) # ensemble AUC = 0.8407
plt.legend(loc='lower right', prop = { "size": 10})
plt.plot([(0, 0), (1, 1)], 'k--', alpha=0.1)
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 01.01])
plt.yticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.xticks(np.array([0,0.2,0.4,0.6,0.8,1.0]), fontsize=15)
plt.ylabel('True Positive Rate', fontsize=15)
plt.xlabel('False Positive Rate', fontsize=15)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/prob_qvalue_roc_ensemble_glcp.png',format='png',dpi=300)
plt.savefig('/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/prob_qvalue_roc_ensemble_glcp.pdf',format='pdf',dpi=300)
plt.close()

