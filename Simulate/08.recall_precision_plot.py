import os
from itertools import combinations
from pathlib import Path
import uuid
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve
import numpy as np

GENE_ID_COL_NAME = 'gene_id'
RESULT_TYPE_PVAL = 'pval'
RESULT_TYPE_PROB = 'prob'
# list of (tool, key_col_in_tool_result, tool_result_type)
TOOL_SIG_COL_INFO = [('coloc', 'overall_H4', RESULT_TYPE_PROB),
                     ('fastenloc', 'GRCP', RESULT_TYPE_PROB),
                     ('smr', 'PROB', RESULT_TYPE_PROB),
                     ('predixcan', 'PROB', RESULT_TYPE_PROB),
                     ('ecaviar', 'clpp', RESULT_TYPE_PROB),
                     ('twas', 'PROB', RESULT_TYPE_PROB)
                     ]


AVG_RANKING_COL_NAME = 'avg_ranking'

prob_col_name = 'PROB'
is_positive_col_name = 'IS_POSITIVE'

tools_threshold = {'coloc': 0.75, 
                 'ecaviar': 0.01,
                 'predixcan': 0.8,
                 'twas': 0.8,
                 'smr': 0.8,
                 'fastenloc': 0.5}


def retrieve_std_df(generated_list):
    std_df = pd.read_table(generated_list, header=None, sep=r'\s+', usecols=[0, 3])
    std_df.columns = [GENE_ID_COL_NAME, is_positive_col_name]
    # std_df[GENE_ID_COL_NAME] = std_df[GENE_ID_COL_NAME] + '_' + std_df[GENE_ID_COL_NAME]
    return std_df

def prepare_plot_data(generated_list, h1_report,
                      rpt_prob_col_name=None, rpt_pval_col_name=None, tool=None):
    if rpt_prob_col_name is None and rpt_pval_col_name is None:
        raise ValueError('p-value column name and probability column name can not be both null')
    std_df = retrieve_std_df(generated_list)

    if rpt_prob_col_name is None:
        reading_cols = [rpt_pval_col_name, GENE_ID_COL_NAME]
        if tool == 'predixcan':
            reading_cols.append('zscore')
        elif tool == 'twas':
            reading_cols.append('TWAS.Z')
        elif tool == 'smr':
            reading_cols.append('b_SMR')
            reading_cols.append('se_SMR')
            reading_cols.append('p_HEIDI')
        h1_report_df = pd.read_table(h1_report, usecols=reading_cols)
        h1_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        h1_report_df[prob_col_name] = 1 - h1_report_df[rpt_pval_col_name]
    else:
        h1_report_df = pd.read_table(h1_report, usecols=[rpt_prob_col_name, GENE_ID_COL_NAME])
        h1_report_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
        h1_report_df[prob_col_name] = h1_report_df[rpt_prob_col_name]
    report_df = h1_report_df
    # TODO MHY is it right to use how as outer?
    result_df = pd.merge(left=std_df, right=report_df,
                         left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                         how='outer')
    tp_na_bool_series = result_df[is_positive_col_name].isna()
    result_df.loc[result_df[tp_na_bool_series].index, is_positive_col_name] = 0
    prob_na_bool_series = result_df[prob_col_name].isna()
    result_df.loc[result_df[prob_na_bool_series].index, prob_col_name] = 0
    if rpt_prob_col_name is None:
        # TODO MHY here 这里还在补值, plot value不能为NA
        result_df.loc[result_df[prob_na_bool_series].index, rpt_pval_col_name] = 1
    else:
        result_df.loc[result_df[prob_na_bool_series].index, rpt_prob_col_name] = 0
    # result_df has 4 columns: [gene_id, prob_col_name, is_positive_col_name, (rpt_prob_col_name or rpt_pval_col_name)]
    return result_df


def plot_precision_recall_f1_comb_bar(
        generated_file_path=None,
        h1_rpt_obj=None,
        output_figure_prefix=None,
        typ='UNION'):
    rpt_obj = {}
    thresh_obj = tools_threshold
    # total = {}
    for tool, sig_column, sig_type in TOOL_SIG_COL_INFO:
        print(tool)
        h1_rpt = h1_rpt_obj.get(tool)
        tool_df = prepare_plot_data(generated_file_path, h1_rpt,
                                    rpt_prob_col_name=sig_column if sig_type == RESULT_TYPE_PROB else None,
                                    rpt_pval_col_name=sig_column if sig_type == RESULT_TYPE_PVAL else None,
                                    tool=tool)
        
        # total[tool]=tool_df
        tool_df['mark'] = pd.NA
        tool_df['mark'].mask(tool_df[sig_column] > thresh_obj[tool], 1, inplace=True)
        tool_df['mark'].mask(tool_df[sig_column] <= thresh_obj[tool], 0, inplace=True)
        rpt_obj[tool] = tool_df[[GENE_ID_COL_NAME, 'mark', is_positive_col_name]]
    # return total
    tools = [tool for tool, _, _ in TOOL_SIG_COL_INFO]
    count = []
    precision_means = []
    precision_stds = []
    recall_means = []
    recall_stds = []
    f1_means = []
    f1_stds = []
    coms = []
    precisions = []
    recalls = []
    f1s = []
    for n in range(1, len(tools) + 1):
        count.append(n)
        ntool_precisions = []
        ntool_recalls = []
        ntool_f1s = []
        n_coms = []
        for com in combinations(tools, n):
            # union n tool results, for the same gene across diff tools,
            # keep the one with the largest mark if typ=UNION,
            # else keep the one with the smallest mark
            df = pd.concat([rpt_obj[tool] for tool in com])
            df.sort_values(by='mark', ascending=typ.upper() != 'UNION', inplace=True)
            # 某基因有一个工具是1就排在前面，则该基因只要在任意工具是1就会被选择，以确保union
            # 是inter则0排在前面，则该基因必须在六种工具中全是1才会被选择
            df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
            if df.empty:
                ntool_precisions.append(0)
                ntool_recalls.append(0)
                ntool_f1s.append(0)
                continue
            tp = sum((df[f'mark'] == 1) & (df[is_positive_col_name] == 1))
            fp = sum((df[is_positive_col_name] == 0) & (df[f'mark'] == 1))
            fn = sum((df[is_positive_col_name] == 1) & (df[f'mark'] == 0))
            if tp + fp == 0:
                precision = 0
            else:
                precision = tp / (tp + fp)
            ntool_precisions.append(precision)
            if tp + fn == 0:
                recall = 0
            else:
                recall = tp / (tp + fn)
            ntool_recalls.append(recall)
            if 2 * tp + fp + fn == 0:
                f1 = 0
            else:
                f1 = 2 * tp / (2 * tp + fp + fn)
            ntool_f1s.append(f1)
            n_coms.append(com)
        n_union_df = pd.DataFrame({'precision': ntool_precisions, 'recall': ntool_recalls, 'f1': ntool_f1s})
        # return n_union_df
        precision_means.append(n_union_df['precision'].mean())
        precision_stds.append(0 if n_union_df.shape[0] == 1 else n_union_df['precision'].std(ddof=0))
        recall_means.append(n_union_df['recall'].mean())
        recall_stds.append(0 if n_union_df.shape[0] == 1 else n_union_df['recall'].std(ddof=0))
        f1_means.append(n_union_df['f1'].mean())
        f1_stds.append(0 if n_union_df.shape[0] == 1 else n_union_df['f1'].std(ddof=0))
        coms.append(n_coms)
        precisions.append(ntool_precisions)
        recalls.append(ntool_recalls)
        f1s.append(ntool_f1s)

    # merge all genes
    gene_list = []
    for tool in tools:
        gene_list.append(rpt_obj[tool][[GENE_ID_COL_NAME, is_positive_col_name]])
    gene_list_df = pd.concat(gene_list)
    del gene_list
    gene_list_df.drop_duplicates(subset=GENE_ID_COL_NAME, inplace=True)
    # merge all results in rpt_obj
    mark_cols = []
    merged_df = gene_list_df
    for tool in tools:
        merged_df = pd.merge(left=merged_df, right=rpt_obj[tool][[GENE_ID_COL_NAME, 'mark']],
                             left_on=GENE_ID_COL_NAME, right_on=GENE_ID_COL_NAME,
                             how='outer')
        merged_df.rename(columns={'mark': f'{tool}_mark'}, inplace=True)
        mark_cols.append(f'{tool}_mark')

    xticklbl = [str(c) for c in count]
    # precision mean
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(2.5, 5.5))
    plt.rcParams['font.size'] = '16'
    ax.grid(False)
    ax.errorbar(count, precision_means, yerr=precision_stds, fmt='o', color='black',
             ecolor='lightgray', elinewidth=1.5, capsize=4)
    ax.set(xticks=[*count], xticklabels=xticklbl)
    ax.set_ylim(0.75,1.04)
    ax.set_xticklabels(np.array([1,2,3,4,5,6]), fontsize=13)
    ax.set_yticklabels(np.array([0.7,0.75,0.80,0.85,0.9,0.95,1.00]), fontsize=13)
    ax.set_ylabel('Precision', fontsize=13)
    ax.set_xlabel('Num of GIMs', fontsize=13)
    for n,v in zip(count,precision_means):
        ax.text(n-0.35, v+0.01, round(v,2),size=10)
    plt.savefig(f'{output_figure_prefix}_{typ}_precision.png', dpi=300, format='png', bbox_inches='tight')
    plt.savefig(f'{output_figure_prefix}_{typ}_precision.pdf', dpi=300, format='pdf', bbox_inches='tight')
    plt.close()


    # recall mean
    plt.figure().clear()
    fig, ax = plt.subplots(figsize=(2.5, 5.5))
    plt.rcParams['font.size'] = '16'
    ax.grid(False)
    ax.errorbar(count, recall_means, yerr=recall_stds, fmt='o', color='black',
             ecolor='lightgray', elinewidth=1.5, capsize=4)
    ax.set(xticks=[*count], xticklabels=xticklbl)
    ax.set_ylim(0.30,0.6)
    ax.set_xticklabels(np.array([1,2,3,4,5,6]), fontsize=13)
    ax.set_yticklabels(np.array([0.35, 0.4,0.45,0.5,0.55,0.6]), fontsize=13)
    ax.set_ylabel('Recall', fontsize=13)
    ax.set_xlabel('Num of GIMs', fontsize=13)
    for n,v in zip(count,recall_means):
        ax.text(n-0.39, v+0.005, round(v,2),size=10)
    #ax.set_title(f'Recall Mean of different combinations of tools', fontsize=15)
    plt.savefig(f'{output_figure_prefix}_{typ}_recall.png', dpi=300, format='png', bbox_inches='tight')
    plt.savefig(f'{output_figure_prefix}_{typ}_recall.pdf', dpi=300, format='pdf', bbox_inches='tight')
    plt.close()

    ## f1 mean
    #plt.figure().clear()
    #fig, ax = plt.subplots(figsize=(5.5, 5.5))
    #plt.rcParams['font.size'] = '16'
    ## ax.set_facecolor("white")
    ## ax.grid(False)
    #rects = ax.bar(count, f1_means, 0.5, label='F1 Mean', color='tan')
    #ax.bar_label(rects, fmt='{:0.2f}', fontsize=15)
    #ax.errorbar(count, f1_means, yerr=f1_stds, fmt=',', ecolor='grey', capsize=2, elinewidth=1)
    #ax.set(xticks=[*count], xticklabels=xticklbl)
    ## ax.set(yticks=np.array([0,0.2,0.4,0.6,0.8,1.0]))
    ## ax.set_xticklabels(np.array([1,2,3,4,5,6]), fontsize=13)
    ## ax.set_yticklabels(np.array([0,0.1,0.2,0.3,0.4,0.5,0.6]), fontsize=13)
    #ax.set_ylabel('F1', fontsize=15)
    #ax.set_xlabel('Num of tools', fontsize=15)
    #ax.set_title(f'F1 Mean of different combinations of tools', fontsize=15)
    #plt.savefig(f'{output_figure_prefix}_{typ}_f1.png', dpi=300, format='png', bbox_inches='tight')
    #plt.savefig(f'{output_figure_prefix}_{typ}_f1.pdf', dpi=300, format='pdf', bbox_inches='tight')
    #plt.close()
    return n_union_df

n_union_df = plot_precision_recall_f1_comb_bar(
    '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240504/gwas.truth',
    {
        'twas': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newtwas2_merged_qvalue.tsv',
        'smr': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newsmr_merged_qvalue.tsv',
        'coloc': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newcoloc_merged.tsv',
        'ecaviar': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newecaviar_merged.tsv', 
        'fastenloc': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/newfastenlocv3.1_GRCP_merged.tsv',
        'predixcan': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newpredixcan_merged_qvalue.tsv'
    },
    output_figure_prefix='/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/new_sim_without_mv_grcp',
    typ='UNION')
    
    
plot_precision_recall_f1_comb_bar(
    '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240504/gwas.truth',
    {
        'twas': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newtwas2_merged_qvalue.tsv',
        'smr': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newsmr_merged_qvalue.tsv',
        'coloc': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newcoloc_merged.tsv',
        'ecaviar': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newecaviar_merged.tsv', 
        'fastenloc': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/newfastenlocv3.1_GRCP_merged.tsv',
        'predixcan': '/Users/phoebel/github/locuscompare2_private/new_sim_result_20240511/newpredixcan_merged_qvalue.tsv'
    },
    output_figure_prefix='/Users/phoebel/github/locuscompare2_private/new_sim_result_20240922/new_sim_wihtout_mv_grcp',
    typ='INTER')


