import os
import pandas as pd

def main(dest_dir):
    template = '''input:
  eqtl:
    col_name_mapping:
      alt: alt
      beta: beta
      chrom: var_chrom
      gene_id: phe_id
      maf: maf
      position: var_from
      pvalue: nom_pval
      ref: ref
      se: se
      snp: var_id
    file: /home/liulab/datasource/eqtl_input/@.tsv.gz
    sample_size: 500
    sep: '\\t'
    tissue: test_tissue
    type: quant
  eqtl_finemapping_file: /home/liulab/datasource/fastenloc/@.vcf.gz
  genecode: /home/liulab/datasource/gencode.v26.basic.annotation.gtf.gz
  gwas:
    col_name_mapping:
      beta: beta
      chrom: var_chrom
      effect_allele: alt
      other_allele: ref
      position: var_from
      pvalue: nom_pval
      se: se
      snp: var_id
    file: /home/liulab/datasource/gwas_input/@.tsv.gz
    sample_size: 500
    trait: @
    type: quant
    sep: '\\t'
  ld_block_loci_file: /home/liulab/datasource/eur_hg38_ld_block.bed
  twas_model_dir: /home/liulab/datasource/sim_twas_model_6/
  prediction_dir: /home/liulab/datasource/sim_predixcan_4/@
  vcf: /home/liulab/datasource/phased_hg38
neighbour_snp_range: 25
p-value_threshold:
  eqtl: 1.0
  gwas: 1.0
population: EUR
working_dir: /home/liulab/datasource/locuscompare2_twas20240507/@
tools: 
- twas
- predixcan
- fastenloc
- coloc
- smr
- ecaviar
'''
    gene_name_ls = pd.read_csv('LocusCompare2_visualization/Simulate/04.generate_config.py')
    for gene_name in gene_name_ls['gene']:
       cfg = template.replace('@', gene_name)
       with open(os.path.join(dest_dir, f'{gene_name}.configure'), 'w') as out:
           out.write(cfg)
 

#if __name__ == '__main__':
main('/Users/phoebel/configurefiles')



import os
import pandas as pd

def main(dest_dir):
    template = '''input:
  eqtl:
    col_name_mapping:
      alt: alt
      beta: beta
      chrom: var_chrom
      gene_id: phe_id
      maf: maf
      position: var_from
      pvalue: nom_pval
      ref: ref
      se: se
      snp: var_id
    file: /home/liulab/datasource/eqtl_input/@.tsv.gz
    sample_size: 500
    sep: '\\t'
    tissue: test_tissue
    type: quant
  eqtl_finemapping_file: /home/liulab/datasource/fastenloc/@.vcf.gz
  genecode: /home/liulab/datasource/gencode.v26.basic.annotation.gtf.gz
  gwas:
    col_name_mapping:
      beta: beta
      chrom: var_chrom
      effect_allele: alt
      other_allele: ref
      position: var_from
      pvalue: nom_pval
      se: se
      snp: var_id
    file: /home/liulab/datasource/gwas_input/@.tsv.gz
    sample_size: 500
    trait: @
    type: quant
    sep: '\\t'
  ld_block_loci_file: /home/liulab/datasource/eur_hg38_ld_block.bed
  twas_model_dir: /home/liulab/datasource/sim_twas_model_6/
  prediction_dir: /home/liulab/datasource/sim_predixcan_4/@
  vcf: /home/liulab/datasource/phased_hg38
neighbour_snp_range: 25
p-value_threshold:
  eqtl: 1.0
  gwas: 1.0
population: EUR
working_dir: /home/liulab/datasource/locuscompare2_twas20240507/@
tools: 
- predixcan
'''
    gene_name_ls = pd.read_csv('LocusCompare2_visualization/Simulate/predixcan_gene_ls.tsv')
    for gene_name in gene_name_ls[0]:
       cfg = template.replace('@', gene_name)
       with open(os.path.join(dest_dir, f'{gene_name}.predixcan.configure'), 'w') as out:
           out.write(cfg)
 

#if __name__ == '__main__':
main('/Users/phoebel/configurefiles')