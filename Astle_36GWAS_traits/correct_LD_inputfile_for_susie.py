import pandas as pd
input_df = pd.read_csv('~/datasource/4621/input.tsv',sep='\t')
ref = pd.read_csv('/home/liulab/datasource/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz',sep='\t')

ref.index = ref['variant_id'].apply(lambda x: '_'.join(x.split('_')[:2]))
input_df.index = input_df['var_id_']
ref = ref.loc[input_df.index]

ref.index = ref.index + ref['ref']+ref['alt']
input_df.index = input_df['var_id_'] + input_df['other_allele'] + input_df['effect_allele']
ref = ref.loc[input_df.index]
input_df = pd.concat([input_df,ref['rs_id_dbSNP151_GRCh38p7']],join='inner',axis=1)

LD = pd.read_csv('~/datasource/4621/LD.matrix',sep='\t')
input_df.index = input_df['rs_id_dbSNP151_GRCh38p7']
input_df = input_df.loc[LD.index]
input_df.to_csv('~/datasource/4621/input.tsv',sep='\t')