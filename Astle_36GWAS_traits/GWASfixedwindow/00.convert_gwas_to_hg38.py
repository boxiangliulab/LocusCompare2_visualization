
import pandas as pd
import numpy as np
import os

gwasfiledict = {
    'GCST004626':'myeloid_wbc_build37_169219_20161212.tsv.gz',
    'GCST004627':'lymph_build37_171643_20161212.tsv.gz',
    'GCST004628':'irf_build37_170548_20161212.tsv.gz',
    'GCST004629':'neut_build37_170702_20161212.tsv.gz',
    'GCST004630':'mch_build37_172332_20161212.tsv.gz',
    'GCST004631':'baso_p_build37_171996_20161212.tsv.gz',
    'GCST004632':'lymph_p_build37_171748_20161212.tsv.gz',
    'GCST004633':'neut_p_build37_171542_20161212.tsv.gz',
    'GCST004634':'baso_p_gran_build37_170223_20161212.tsv.gz',
    'GCST004599':'mpv_build37_164454_20161212.tsv.gz',
    'GCST004601':'rbc_build37_172952_20161212.tsv.gz',
    'GCST004600':'eo_p_build37_172378_20161212.tsv.gz',
    'GCST004602':'mcv_build37_172433_20161212.tsv.gz',
    'GCST004603':'plt_build37_166066_20161212.tsv.gz',
    'GCST004604':'hct_build37_173039_20161212.tsv.gz',
    'GCST004605':'mchc_build37_172851_20161212.tsv.gz',
    'GCST004606':'eo_build37_172275_20161212.tsv.gz',
    'GCST004607':'pct_build37_164339_20161212.tsv.gz',
    'GCST004608':'gran_p_myeloid_wbc_build37_169545_20161212.tsv.gz',
    'GCST004609':'mono_p_build37_170494_20161212.tsv.gz',
    'GCST004610':'wbc_build37_172435_20161212.tsv.gz',
    'GCST004611':'hlr_build37_170761_20161212.tsv.gz',
    'GCST004613':'neut_eo_sum_build37_170384_20161212.tsv.gz',
    'GCST004614':'gran_build37_169822_20161212.tsv.gz',
    'GCST004612':'hlr_p_build37_170763_20161212.tsv.gz',
    'GCST004615':'hgb_build37_172925_20161212.tsv.gz',
    'GCST004616':'pdw_build37_164433_20161212.tsv.gz',
    'GCST004617':'eo_p_gran_build37_170536_20161212.tsv.gz',
    'GCST004618':'baso_build37_171846_20161212.tsv.gz',
    'GCST004619':'ret_p_build37_170690_20161212.tsv.gz',
    'GCST004620':'baso_neut_sum_build37_170143_20161212.tsv.gz',
    'GCST004621':'rdw_build37_171529_20161212.tsv.gz',
    'GCST004622':'ret_build37_170641_20161212.tsv.gz',
    'GCST004623':'neut_p_gran_build37_170672_20161212.tsv.gz',
    'GCST004624':'eo_baso_sum_build37_171771_20161212.tsv.gz',
    'GCST004625':'mono_build37_170721_20161212.tsv.gz'
    }


def main(gwasfile,ref):
    path = '/home/liulab/datasource/Astle_blood_trait'
    df = pd.read_csv(os.path.join(path, f"{gwasfile}"),sep='\t')
    df = df[df['ID'] != '.']
    df.index = df['ID']
    df = df.loc[list(set(df['ID']) & set(ref['rs_id_dbSNP151_GRCh38p7']))]
    ref.index = ref['rs_id_dbSNP151_GRCh38p7']
    ref_ = ref.loc[list(set(df['ID']) & set(ref['rs_id_dbSNP151_GRCh38p7']))]
    ref_ = ref_.drop_duplicates('variant_id')
    ref_ = ref_.drop_duplicates('rs_id_dbSNP151_GRCh38p7')
    df = df.sort_values('P')
    df = df.drop_duplicates('ID')
    df = pd.concat([df, ref_['variant_id']], join='inner',axis=1)
    df['ref'] = df['variant_id'].apply(lambda x: x.split('_')[2])
    df['alt'] = df['variant_id'].apply(lambda x: x.split('_')[3])
    df['chr'] = df['variant_id'].apply(lambda x: str(x.split('_')[0][3:]))
    df['pos'] = df['variant_id'].apply(lambda x: x.split('_')[1])
    df = df[['CHR','pos','ref','alt','EFFECT','P','ID','SE']]

    df.columns = ['chr','pos','ref','alt','beta','pval','rsid','se']
    df.to_csv(os.path.join(path, f"preprocessed_{gwasfile}"),sep='\t',compression='gzip',index=False)



ref = pd.read_csv('/home/liulab/datasource/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz',sep='\t')


for i in gwasfiledict.keys():
    main(gwasfiledict[i], ref)
