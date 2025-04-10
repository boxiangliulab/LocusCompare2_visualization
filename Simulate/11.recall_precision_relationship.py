
def plot_precision_recall_f1_comb_bar(
        generated_file_path=None,
        h1_rpt_obj=None,
        output_figure_prefix=None,
        type='UNION'):
    rpt_obj = {}
    thresh_obj = calc_threshold(h1_rpt_obj)
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
        if sig_type == RESULT_TYPE_PVAL:
            if tool == 'smr':
                positive_series = (tool_df[sig_column] < thresh_obj[tool]) & (tool_df['p_HEIDI'] > 0.05)
                tool_df['mark'].mask(positive_series, 1, inplace=True)
                tool_df['mark'].mask(~positive_series, 0, inplace=True)
            else:
                tool_df['mark'].mask(tool_df[sig_column] < thresh_obj[tool], 1, inplace=True)
                tool_df['mark'].mask(tool_df[sig_column] >= thresh_obj[tool], 0, inplace=True)
        else:
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
            # keep the one with the largest mark if type=UNION,
            # else keep the one with the smallest mark
            df = pd.concat([rpt_obj[tool] for tool in com])
            df.sort_values(by='mark', ascending=type.upper() != 'UNION', inplace=True)
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

    ################################### plot ###################################
    fig, ax = plt.subplots(figsize=(5, 5))  # 增加图的宽度以便显示更多内容
    annos = ['1 GIM', '2 GIMs', '3 GIMs', '4 GIMs', '5 GIMs', '6 GIMs']

    ax.errorbar(recall_means, precision_means, xerr=recall_stds, yerr=precision_stds, color='tan',
            fmt='o', ecolor='lightgray', elinewidth=1, capsize=3)
    ax.scatter(recall_means, precision_means, marker='X', color='tan', edgecolors='grey', s=200)
    ax.set_ylabel(f'Precision', fontsize=16)
    ax.set_xlabel('Recall', fontsize=16)
    ax.tick_params(axis='both', labelsize=14)  
    ax.set_title(f"{type}", fontsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for i, txt in enumerate(annos):
        ax.annotate(txt, (recall_means[i]+0.01, precision_means[i]+0.01), fontsize=9, ha='right')


    plt.savefig(f'{output_figure_prefix}_{type}_combined_precision_recall.pdf', dpi=300, format='pdf', bbox_inches='tight')


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
    output_figure_prefix='/Users/phoebel/github/locuscompare2_private/rebuttle_20250121/relationship',
    type='UNION')
    
    
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
    output_figure_prefix='/Users/phoebel/github/locuscompare2_private/rebuttle_20250121/relationship',
    type='INTER')
    