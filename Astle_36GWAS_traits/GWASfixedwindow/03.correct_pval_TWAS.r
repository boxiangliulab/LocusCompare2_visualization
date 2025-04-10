library(qvalue)

df = read.csv('~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/fusion_qvalue.tsv.gz',sep='\t')
df$fusionqvalue = qvalue(df$TWAS.P, fdr.level=0.05)$qvalues
write.table(df, '~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/fusion_qvalue.tsv', sep='\t', quote=F, row.names=F)

df = read.csv('~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/smr_qvalue.tsv.gz',sep='\t')
df$smrqvalue = qvalue(df$p_SMR, fdr.level=0.05)$qvalues
write.table(df, '~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/smr_qvalue.tsv', sep='\t', quote=F, row.names=F)

df = read.csv('~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/predixcan_qvalue.tsv.gz',sep='\t')
df$predixcanqvalue = qvalue(df$pvalue, fdr.level=0.05)$qvalues
write.table(df, '~/locuscompare2_private/rebuttle_20250121/fix_gwasloci_window/predixcan_qvalue.tsv', sep='\t', quote=F, row.names=F)

