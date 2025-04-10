import pandas as pd
import os
from datetime import datetime
import scipy.stats




def convert_sbams_to_split(sbams_dir, dest_dir,
                           genecode='/data/projects/11003054/locuscompare2file/gtf_all_gene.tsv',
                           individual_num=500,
                           dbsnp_ref_path='/data/projects/11003054/locuscompare2file/vcf_hg38/EUR/chr5.vcf.gz'):
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    sum_dir = os.path.join(dest_dir, 'sum')
    if not os.path.exists(sum_dir):
        os.makedirs(sum_dir)
    ind_ids = [f'id{k}' for k in range(1, individual_num + 1)]
    genecode_df = pd.read_table(genecode, usecols=['start', 'end', 'strand', 'gene_id'])
    genecode_df.rename({'gene_id': 'gid'}, axis='columns', inplace=True)
    iid_mapper = {k: f'id{k - 2}' for k in range(3, individual_num + 3)}
    for sbam in os.listdir(sbams_dir):
        if not sbam.endswith('.sbams.dat'):
            continue
        sbam_path = os.path.join(sbams_dir, sbam)

        dosage_df = pd.read_table(sbam_path, sep=r'\s+', header=None, skiprows=1,
                                  dtype={k: 'str' for k in range(3, individual_num + 3)})

        split_df = dosage_df[1].str.split("_", n=4, expand=True)
        chrom = int(split_df[0].loc[0].strip('chr'))
        dosage_df.drop(columns=[0, 2], inplace=True)
        dosage_df['FORMAT'] = 'GT'
        dosage_df['INFO'] = '.'
        dosage_df['FILTER'] = '.'
        dosage_df['QUAL'] = '.'
        dosage_df['ALT'] = split_df[3]
        dosage_df['REF'] = split_df[2]
        dosage_df['POS'] = split_df[1].astype(int)
        del split_df
        dosage_df['#CHROM'] = chrom
        dosage_mapper = iid_mapper.copy()
        dosage_mapper[1] = 'ID'
        dosage_df.rename(dosage_mapper, axis='columns', inplace=True)
        dosage_df = dosage_df.reindex(
            columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + ind_ids, copy=False)
        dosage_df.replace({'0': '0/0', '1': '0/1', '2': '1/1'}, inplace=True)
        dosage_df.fillna('./.', inplace=True)

        phen_df = pd.read_table(sbam_path, sep=r'\s+', header=None, nrows=1)
        gene_id = phen_df.iat[0, 1]
        pre_vcf_out = os.path.join(dest_dir, f'pre_{gene_id}.vcf')
        vcf_out = os.path.join(dest_dir, f'{gene_id}.vcf.gz')
        with open(pre_vcf_out, mode='w') as vcf_file:
            vcf_file.write(f'##fileformat=VCFv4.2\n')
            vcf_file.write(f'##FILTER=<ID=PASS,Description="All filters passed">\n')
            vcf_file.write(f'##fileDate={datetime.now().strftime("%Y%m%d")}\n')
            vcf_file.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            vcf_file.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ind_ids) + '\n')
        dosage_df.drop_duplicates(subset=['#CHROM', 'POS'], inplace=True)
        dosage_df.sort_values('POS', inplace=True)
        min_pos = dosage_df['POS'].min()
        max_pos = dosage_df['POS'].max()
        dosage_df.to_csv(pre_vcf_out, sep='\t', header=False, index=False, mode='a')
        del dosage_df
        os.system(f'bgzip -f {pre_vcf_out} && tabix -f -p vcf {pre_vcf_out}.gz && '
                  f'bcftools annotate -r {chrom}:{min_pos}-{max_pos} -a {dbsnp_ref_path} -c ID {pre_vcf_out}.gz -Oz -o {vcf_out} && tabix -f -p vcf {vcf_out}')
        # os.remove(f'{pre_vcf_out}.gz')
        # os.remove(f'{pre_vcf_out}.gz.tbi')
        phen_df.drop(columns=[0, 2], inplace=True)
        phen_mapper = iid_mapper.copy()
        phen_mapper[1] = 'gid'
        phen_df.rename(phen_mapper, axis='columns', inplace=True)
        phen_df = pd.merge(left=phen_df, right=genecode_df, left_on='gid', right_on='gid', how='left')
        phen_df['#chr'] = chrom
        phen_df['pid'] = phen_df['gid']
        phen_df = phen_df.reindex(columns=['#chr', 'start', 'end', 'pid', 'gid', 'strand'] + ind_ids, copy=False)
        bed_out = os.path.join(dest_dir, f'{gene_id}.bed')
        with open(bed_out, mode='w') as bed_file:
            bed_file.write(f'#chr\tstart\tend\tpid\tgid\tstrand\t' + '\t'.join(ind_ids) + '\n')
        phen_df.drop_duplicates(subset='gid', inplace=True)
        phen_df.sort_values('start', inplace=True)
        phen_df.to_csv(bed_out, sep='\t', header=False, index=False, mode='a')
        del phen_df
        os.system(f'bgzip -f {bed_out} && tabix -f -p bed {bed_out}.gz')
        pre_sum_out = os.path.join(dest_dir, f'{gene_id}_sum')
        # os.system(
        #     f'QTLtools_1.2_Ubuntu16.04_x86_64 cis --vcf {vcf_out} --bed {bed_out}.gz --nominal 1 --normal --std-err --out {pre_sum_out}')
        os.system(
            f'QTLtools_1.2_Ubuntu16.04_x86_64 cis --vcf {vcf_out} --bed {bed_out}.gz --nominal 1 --normal --out {pre_sum_out}')        

        if not os.path.exists(pre_sum_out):
            print(f'Warning: QTL analysis result of {gene_id} does not exist')
            continue
        ana_result_df = pd.read_table(pre_sum_out, header=None, sep=r'\s+')
        os.remove(pre_sum_out)
        if ana_result_df.shape[0] == 0:
            print(f'Warning: QTL analysis result of {gene_id} is empty')
            continue
        # ana_result_df.columns = ['phe_id', 'phe_chrom', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis',
        #                          'dist_phe_var',
        #                          'var_id', 'var_chrom', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'beta', 'se',
        #                          'best_hit']
        ana_result_df.columns = ['phe_id', 'phe_chrom', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis',
                                 'dist_phe_var',
                                 'var_id', 'var_chrom', 'var_from', 'var_to', 'nom_pval', 'beta',
                                 'best_hit']
        ana_result_df['se'] = abs(ana_result_df['beta']/scipy.stats.norm.ppf(ana_result_df['nom_pval']/2))
        os.system(f'plink2 --vcf {vcf_out} --freq --double-id --out  {os.path.join(dest_dir,gene_id)}_maf')
        # eqtl_maf_df = pd.read_table(f'{os.path.join(dest_dir,gene_id)}_maf.afreq', sep=r'\s+', usecols=['SNP', 'MAF'])
        eqtl_maf_df = pd.read_table(f'{os.path.join(dest_dir,gene_id)}_maf.afreq', sep=r'\s+', usecols=['ID', 'ALT_FREQS'])
        os.remove(f'{os.path.join(dest_dir,gene_id)}_maf.afreq')
        # os.remove(f'{os.path.join(dest_dir,gene_id)}_maf.nosex')
        os.remove(f'{os.path.join(dest_dir,gene_id)}_maf.log')
        eqtl_maf_df.columns = ['var_id', 'maf']
        ana_result_df = pd.merge(ana_result_df, eqtl_maf_df,
                                 left_on='var_id',
                                 right_on='var_id',
                                 how='left')
        # merge REF, ALT from genotype vcf
        genotype_df = pd.read_table(vcf_out, header=None, comment='#', usecols=[2, 3, 4])
        genotype_df.columns = ['var_id', 'ref', 'alt']
        ana_result_df = pd.merge(ana_result_df, genotype_df,
                                 left_on='var_id',
                                 right_on='var_id',
                                 how='left')
        sum_out = os.path.join(sum_dir, f'{gene_id}.tsv.gz')
        ana_result_df.to_csv(sum_out, sep='\t', header=True, index=False)

convert_sbams_to_split('/data/projects/11003054/locuscompare2file/sim_scheme_1/sim_data/eqtl_sbams', '/data/projects/11003054/locuscompare2file/sim_eqtl_20240422')
convert_sbams_to_split('/data/projects/11003054/locuscompare2file/sim_scheme_1/sim_data/gwas_sbams', '/data/projects/11003054/locuscompare2file/sim_gwas_20240422')


def convert_sbams_to_merged(sbams_dir, dest_dir,prefix='out',
                           genecode='/data/projects/11003054/locuscompare2file/gtf_all_gene.tsv',
                           individual_num=500,
                           dbsnp_ref_path='/data/projects/11003054/locuscompare2file/vcf_hg38/EUR/chr5.vcf.gz'):
    pre_vcf_out = os.path.join(dest_dir, f'pre_{prefix}.vcf')
    bed_out = os.path.join(dest_dir, f'{prefix}.bed')
    ind_ids = [f'id{k}' for k in range(1, individual_num + 1)]
    genecode_df = pd.read_table(genecode, usecols=['start', 'end', 'strand', 'gene_id'])
    genecode_df.rename({'gene_id': 'gid'}, axis='columns', inplace=True)
    iid_mapper = {k: f'id{k - 2}' for k in range(3, individual_num + 3)}
    phen_df_list = []
    dosage_df_list = []
    for sbam in os.listdir(sbams_dir):
        if not sbam.endswith('.sbams.dat'):
            continue
        sbam_path = os.path.join(sbams_dir, sbam)

        dosage_df = pd.read_table(sbam_path, sep=r'\s+', header=None, skiprows=1,
                                  dtype={k: 'str' for k in range(3, individual_num + 3)})

        split_df = dosage_df[1].str.split("_", n=4, expand=True)
        chrom = int(split_df[0].loc[0].strip('chr'))
        dosage_df.drop(columns=[0, 2], inplace=True)
        dosage_df['FORMAT'] = 'GT'
        dosage_df['INFO'] = '.'
        dosage_df['FILTER'] = '.'
        dosage_df['QUAL'] = '.'
        dosage_df['ALT'] = split_df[3]
        dosage_df['REF'] = split_df[2]
        dosage_df['POS'] = split_df[1].astype(int)
        del split_df
        dosage_df['#CHROM'] = chrom
        dosage_mapper = iid_mapper.copy()
        dosage_mapper[1] = 'ID'
        dosage_df.rename(dosage_mapper, axis='columns', inplace=True)
        dosage_df = dosage_df.reindex(
            columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + ind_ids, copy=False)
        dosage_df.replace({'0': '0/0', '1': '0/1', '2': '1/1'}, inplace=True)
        dosage_df.fillna('./.', inplace=True)
        dosage_df_list.append(dosage_df)

        phen_df = pd.read_table(sbam_path, sep=r'\s+', header=None, nrows=1)
        # gene_id = phen_df.iat[0, 1]
        phen_df.drop(columns=[0, 2], inplace=True)
        phen_mapper = iid_mapper.copy()
        phen_mapper[1] = 'gid'
        phen_df.rename(phen_mapper, axis='columns', inplace=True)
        phen_df = pd.merge(left=phen_df, right=genecode_df, left_on='gid', right_on='gid', how='left')
        phen_df['#chr'] = chrom
        phen_df['pid'] = phen_df['gid']
        phen_df = phen_df.reindex(columns=['#chr', 'start', 'end', 'pid', 'gid', 'strand'] + ind_ids, copy=False)
        phen_df_list.append(phen_df)
    with open(pre_vcf_out, mode='w') as vcf_file:
        vcf_file.write(f'##fileformat=VCFv4.2\n')
        vcf_file.write(f'##FILTER=<ID=PASS,Description="All filters passed">\n')
        vcf_file.write(f'##fileDate={datetime.now().strftime("%Y%m%d")}\n')
        vcf_file.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf_file.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ind_ids) + '\n')
    with open(bed_out, mode='w') as bed_file:
        bed_file.write(f'#chr\tstart\tend\tpid\tgid\tstrand\t' + '\t'.join(ind_ids) + '\n')
    phen = pd.concat(phen_df_list)
    phen.drop_duplicates(subset='gid', inplace=True)
    phen.sort_values('start', inplace=True)
    phen.to_csv(bed_out, sep='\t', header=False, index=False, mode='a')
    del phen
    dosage = pd.concat(dosage_df_list)
    dosage.drop_duplicates(subset=['#CHROM', 'POS'], inplace=True)
    dosage.sort_values('POS', inplace=True)
    dosage.to_csv(pre_vcf_out, sep='\t', header=False, index=False, mode='a')
    del dosage
    #
    vcf_out = os.path.join(dest_dir, f'{prefix}.vcf.gz')
    os.system(f'bgzip -f {pre_vcf_out} && tabix -f -p vcf {pre_vcf_out}.gz && '
              f'bcftools annotate -a {dbsnp_ref_path} -c ID {pre_vcf_out}.gz -Oz -o {vcf_out} && tabix -f -p vcf {vcf_out}')
    os.remove(f'{pre_vcf_out}.gz')
    os.remove(f'{pre_vcf_out}.gz.tbi')

    os.system(f'bgzip -f {bed_out} && tabix -f -p bed {bed_out}.gz')
    pca_out = os.path.join(dest_dir, f'{prefix}.pca')
    os.system(
        f'QTLtools_1.2_Ubuntu16.04_x86_64 pca --bed {bed_out}.gz --out {os.path.join(dest_dir, prefix)}_raw_pca1 --center --scale && head -4 {os.path.join(dest_dir, prefix)}_raw_pca1.pca > {pca_out} && gzip -k {pca_out}')
    # os.remove(f'{os.path.join(dest_dir,prefix)}_raw_pca1')
    if not os.path.exists(pca_out):
        raise ValueError('pca not generated')
    pre_sum_out = os.path.join(dest_dir,f'{prefix}_sum')
    os.system(
        f'QTLtools_1.2_Ubuntu16.04_x86_64 cis --vcf {vcf_out} --bed {bed_out}.gz --cov {pca_out} --nominal 1 --normal --out {pre_sum_out}')
    if not os.path.exists(pre_sum_out):
        raise ValueError('sumstats not generated')
    sumstats_df = pd.read_table(pre_sum_out, header=None, sep=r'\s+')
    if sumstats_df.shape[0] == 0:
        raise ValueError('sumstats is empty')
    os.remove(pre_sum_out)
    sumstats_df.columns = ['phe_id', 'phe_chrom', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis',
                            'dist_phe_var',
                            'var_id', 'var_chrom', 'var_from', 'var_to', 'nom_pval', 'beta',
                            'best_hit']
    sumstats_df['se'] = abs(sumstats_df['beta']/scipy.stats.norm.ppf(sumstats_df['nom_pval']/2))
    os.system(f'plink2 --vcf {vcf_out} --freq --double-id --out  {os.path.join(dest_dir,prefix)}_maf')
    maf_df = pd.read_table(f'{os.path.join(dest_dir,prefix)}_maf.afreq', sep=r'\s+', usecols=['ID', 'ALT_FREQS'])
    os.remove(f'{os.path.join(dest_dir,prefix)}_maf.afreq')
    # os.remove(f'{os.path.join(dest_dir,gene_id)}_maf.nosex')
    os.remove(f'{os.path.join(dest_dir,prefix)}_maf.log')
    maf_df.columns = ['var_id', 'maf']
    sumstats_df = pd.merge(sumstats_df, maf_df,
                           left_on='var_id',
                           right_on='var_id',
                           how='left')
    # merge REF, ALT from genotype vcf
    genotype_df = pd.read_table(vcf_out, header=None, comment='#', usecols=[2, 3, 4])
    genotype_df.columns = ['var_id', 'ref', 'alt']
    sumstats_df = pd.merge(sumstats_df, genotype_df,
                           left_on='var_id',
                           right_on='var_id',
                           how='left')
    sum_out = os.path.join(dest_dir, f'{prefix}.tsv.gz')
    sumstats_df.to_csv(sum_out, sep='\t', header=True, index=False)
    
    
convert_sbams_to_merged('/data/projects/11003054/locuscompare2file/sim_scheme_1/sim_data/eqtl_sbams', '/data/projects/11003054/locuscompare2file/sim_eqtl_20240430')
convert_sbams_to_merged('/data/projects/11003054/locuscompare2file/sim_scheme_1/sim_data/gwas_sbams', '/data/projects/11003054/locuscompare2file/sim_gwas_20240430')



