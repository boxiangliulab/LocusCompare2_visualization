library(LDlinkR)
library(coloc)
library(susieR)


gwas_sample_size = 171529
eqtl_sample_size = 670
gwas_df = read.table(file = '/home/liulab/datasource/4621/gwas_ENSG00000245937_chr5_128083158.tsv.gz', header = T)
eqtl_df = read.table(file = '/home/liulab/datasource/4621/eqtl_ENSG00000245937_chr5_128083158.tsv.gz', header = T)
gwas_type = 'quant'
eqtl_type = 'quant'

eqtl_df = eqtl_df[match(gwas_df$var_id_, eqtl_df$var_id_),]

input = merge(gwas_df, eqtl_df, by = "var_id_", all = FALSE, suffixes = c("_gwas", "_eqtl"))
input = subset(input, position_gwas>128000000, position_gwas<128200000)


rownames(input) = gsub("_",":",input$var_id_)


LD = LDmatrix(snps = rownames(input), 
         pop = "EUR", 
         r2d = "r2", 
         token = 'token',
         genome_build = "grch38",
        )
LD[is.na(LD)] <- 0
rownames(LD) = LD$RS_number
LD = (LD[,2:1217])
input = input[rownames(LD),]
write.table(LD, '/home/liulab/datasource/4621/LD.matrix',sep='\t',quote=F)
write.table(input,'~/datasource/4621/input.tsv',sep='\t',quote=F)

# use python to correct rsid: LocusCompare2_visualization/Astle_36GWAS_traits/correct_LD_inputfile_for_susie.py

input = read.table(file = '~/datasource/4621/input.tsv', header = T)
LD = read.table(file = '~/datasource/4621/LD.matrix', header = T)
d1 = list(snp = input$rs_id_dbSNP151_GRCh38p7, beta = input$beta_gwas, varbeta = input$varbeta_gwas, position = input$position_gwas, type = gwas_type, N = gwas_sample_size, MAF = as.numeric(input$maf), LD=as.matrix(LD))
d2 = list(snp = input$rs_id_dbSNP151_GRCh38p7, beta = input$beta_eqtl, varbeta = input$varbeta_eqtl, position = input$position_eqtl, type = eqtl_type, N = eqtl_sample_size, MAF = as.numeric(input$maf), LD=as.matrix(LD))


my.res <- coloc.abf(dataset1=d1, dataset2=d2)


result = susie_rss(bhat = input$beta_gwas,shat = input$se_gwas,R = as.matrix(LD),n=171529,L = 10,estimate_residual_variance = F,coverage = 0.95)
