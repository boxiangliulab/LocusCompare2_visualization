import os
import pandas as pd

tissue_dict = {
    'Whole_Blood':'670',
    }


gwasfiledict = {
    'GCST004626':'preprocessed_myeloid_wbc_build37_169219_20161212.tsv.gz',
    'GCST004627':'preprocessed_lymph_build37_171643_20161212.tsv.gz',
    'GCST004628':'preprocessed_irf_build37_170548_20161212.tsv.gz',
    'GCST004629':'preprocessed_neut_build37_170702_20161212.tsv.gz',
    'GCST004630':'preprocessed_mch_build37_172332_20161212.tsv.gz',
    'GCST004631':'preprocessed_baso_p_build37_171996_20161212.tsv.gz',
    'GCST004632':'preprocessed_lymph_p_build37_171748_20161212.tsv.gz',
    'GCST004633':'preprocessed_neut_p_build37_171542_20161212.tsv.gz',
    'GCST004634':'preprocessed_baso_p_gran_build37_170223_20161212.tsv.gz',
    'GCST004599':'preprocessed_mpv_build37_164454_20161212.tsv.gz',
    'GCST004601':'preprocessed_rbc_build37_172952_20161212.tsv.gz',
    'GCST004600':'preprocessed_eo_p_build37_172378_20161212.tsv.gz',
    'GCST004602':'preprocessed_mcv_build37_172433_20161212.tsv.gz',
    'GCST004603':'preprocessed_plt_build37_166066_20161212.tsv.gz',
    'GCST004604':'preprocessed_hct_build37_173039_20161212.tsv.gz',
    'GCST004605':'preprocessed_mchc_build37_172851_20161212.tsv.gz',
    'GCST004606':'preprocessed_eo_build37_172275_20161212.tsv.gz',
    'GCST004607':'preprocessed_pct_build37_164339_20161212.tsv.gz',
    'GCST004608':'preprocessed_gran_p_myeloid_wbc_build37_169545_20161212.tsv.gz',
    'GCST004609':'preprocessed_mono_p_build37_170494_20161212.tsv.gz',
    'GCST004610':'preprocessed_wbc_build37_172435_20161212.tsv.gz',
    'GCST004611':'preprocessed_hlr_build37_170761_20161212.tsv.gz',
    'GCST004613':'preprocessed_neut_eo_sum_build37_170384_20161212.tsv.gz',
    'GCST004614':'preprocessed_gran_build37_169822_20161212.tsv.gz',
    'GCST004612':'preprocessed_hlr_p_build37_170763_20161212.tsv.gz',
    'GCST004615':'preprocessed_hgb_build37_172925_20161212.tsv.gz',
    'GCST004616':'preprocessed_pdw_build37_164433_20161212.tsv.gz',
    'GCST004617':'preprocessed_eo_p_gran_build37_170536_20161212.tsv.gz',
    'GCST004618':'preprocessed_baso_build37_171846_20161212.tsv.gz',
    'GCST004619':'preprocessed_ret_p_build37_170690_20161212.tsv.gz',
    'GCST004620':'preprocessed_baso_neut_sum_build37_170143_20161212.tsv.gz',
    'GCST004621':'preprocessed_rdw_build37_171529_20161212.tsv.gz',
    'GCST004622':'preprocessed_ret_build37_170641_20161212.tsv.gz',
    'GCST004623':'preprocessed_neut_p_gran_build37_170672_20161212.tsv.gz',
    'GCST004624':'preprocessed_eo_baso_sum_build37_171771_20161212.tsv.gz',
    'GCST004625':'preprocessed_mono_build37_170721_20161212.tsv.gz'
    }

gwassamplesize = {

    'GCST004599':'164454',
    'GCST004601':'172952',
    'GCST004600':'172378',
    'GCST004602':'172433',
    'GCST004603':'166066',
    'GCST004604':'173039',
    'GCST004605':'172851',
    'GCST004606':'172275',
    'GCST004607':'164339',
    'GCST004608':'169545',
    'GCST004609':'170494',
    'GCST004610':'172435',
    'GCST004611':'170761',
    'GCST004612':'170763',
    'GCST004613':'170384',
    'GCST004614':'169822',
    'GCST004615':'172925',
    'GCST004616':'164433',
    'GCST004617':'170536',
    'GCST004618':'171846',
    'GCST004619':'170690',
    'GCST004620':'170143',
    'GCST004621':'171529',
    'GCST004622':'170641',
    'GCST004623':'170672',
    'GCST004624':'171771',
    'GCST004625':'170721',
    'GCST004626':'169219',
    'GCST004627':'171643',
    'GCST004628':'170548',
    'GCST004629':'170702',
    'GCST004630':'172332',
    'GCST004631':'171996',
    'GCST004632':'171748',
    'GCST004633':'171542',
    'GCST004634':'170223',
    }
gwastraittype = {
    'GCST004626':'quant',
    'GCST004627':'quant',
    'GCST004628':'quant',
    'GCST004629':'quant',
    'GCST004630':'quant',
    'GCST004631':'quant',
    'GCST004632':'quant',
    'GCST004633':'quant',
    'GCST004634':'quant',
    'GCST004599':'quant',
    'GCST004601':'quant',
    'GCST004600':'quant',
    'GCST004602':'quant',
    'GCST004603':'quant',
    'GCST004604':'quant',
    'GCST004605':'quant',
    'GCST004606':'quant',
    'GCST004607':'quant',
    'GCST004608':'quant',
    'GCST004609':'quant',
    'GCST004610':'quant',
    'GCST004611':'quant',
    'GCST004613':'quant',
    'GCST004614':'quant',
    'GCST004612':'quant',
    'GCST004615':'quant',
    'GCST004616':'quant',
    'GCST004617':'quant',
    'GCST004618':'quant',
    'GCST004619':'quant',
    'GCST004620':'quant',
    'GCST004621':'quant',
    'GCST004622':'quant',
    'GCST004623':'quant',
    'GCST004624':'quant',
    'GCST004625':'quant'
    }







def main(dest_dir):

    template = '''input:
  eqtl:
    col_name_mapping:
      alt: alt
      beta: beta
      chrom: chromosome
      gene_id: molecular_trait_id
      maf: maf
      position: position
      pvalue: pvalue
      ref: ref
      se: se
      snp: rsid
    file: /home/users/nus/e1124850/scratch/tmp/TISSUE.v8.all_association.processed.txt.gz
    eqtl_preprocessed_dir: /home/users/nus/e1124850/scratch/tmp/locuscompare2file/preprocessed/
    sample_size: EQTLSAMPLESIZE
    sep: '\\t'
    tissue: TISSUE
    type: quant
  eqtl_finemapping_file: /data/projects/11003054/e1124850/locuscompare2file/v8.vcf.gz
  genecode: /data/projects/11003054/e1124850/locuscompare2file/gencode.v45.basic.annotation.gtf.gz
  gwas:
    col_name_mapping:
      beta: beta
      chrom: chr
      effect_allele: alt
      other_allele: ref
      position: pos
      pvalue: pval
      se: se
      snp: rsid
    file: /data/projects/11003054/e1124850/locuscompare2file/Astle_blood_trait/GWASFILEPATH
    sample_size: GWASSAMPLESIZE
    trait: GWAS_TRAIT
    type: GWASTRAITTYPE
    sep: '\\t'
  ld_block_loci_file: /data/projects/11003054/e1124850/locuscompare2file/eur_hg38_ld_block.bed
  prediction_dir: /data/projects/11003054/e1124850/locuscompare2file/predixcan_model/mashr
  twas_model_dir: /data/projects/11003054/e1124850/locuscompare2file/GTEx_twas_model
  vcf: /data/projects/11003054/e1124850/locuscompare2file/vcf_hg38
p-value_threshold:
  eqtl: 1.0
  gwas: 5.0e-8
population: EUR
working_dir: /data/projects/11003054/e1124850/locuscompare2file/Astle_blood_trait/OUTPUTPATH
tools: 
- fusion
- predixcan
- smr
- coloc
- fastenloc
- ecaviar

'''

    pbstemplat='''#PBS -q normal
#PBS -l select=1:ncpus=5:mem=300G
#PBS -l walltime=24:00:00
#PBS -P Personal
#PBS -N CONFIGURE
#PBS -o /home/users/nus/e1124850/scratch/qsub_dir/CONFIGURE.o
#PBS -e /home/users/nus/e1124850/scratch/qsub_dir/CONFIGURE.e



module load miniforge3/23.10
module load cuda/12.2.2
conda activate colotools


python /home/users/nus/e1124850/locuscompare2/locuscompare2-standalone/colotools.py \
--config /data/projects/11003054/e1124850/locuscompare2file/Astle_blood_trait/config/CONFIGURE \
--log /home/users/nus/e1124850/scratch/tmp/CONFIGURE.log


'''
    # gene_name_ls = pd.read_csv('/data/projects/11003054/e1124850/locuscompare2file/matrixeqtl_file/gene_name.tsv')
    #gene_name_ls = pd.read_csv('/Users/phoebel/gene_name.tsv')
    #for gene_name in gene_name_ls['gene']:
    #    cfg = template.replace('@', gene_name)
    #    with open(os.path.join(dest_dir, f'{gene_name}.configure'), 'w') as out:
    #        out.write(cfg)

    
    for gwas_ in gwasfiledict.keys():  
        cfg = template.replace('GWASFILEPATH', gwasfiledict[gwas_])
        cfg = cfg.replace('GWASSAMPLESIZE', gwassamplesize[gwas_])
        cfg = cfg.replace('GWASTRAITTYPE', gwastraittype[gwas_])
        cfg = cfg.replace('GWAS_TRAIT', gwas_)
        cfg = cfg.replace('OUTPUTPATH', gwas_)
        for tissue in tissue_dict.keys():
            cfg_ = cfg.replace('TISSUE', tissue)
            cfg_ = cfg_.replace('EQTLSAMPLESIZE', tissue_dict[tissue])

            with open(os.path.join(dest_dir, f'{gwas_}.configure'), 'w') as out:
                out.write(cfg_)
                
            pbs = pbstemplat.replace('CONFIGURE', f'{gwas_}.configure')
            with open(os.path.join('/Users/phoebel/pbsjob', f'{gwas_}.{tissue}.pbs'), 'w') as out:
                out.write(pbs)
            
            

#if __name__ == '__main__':
main('/Users/phoebel/configurefiles')
