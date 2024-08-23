library(qvalue)

df = read.csv('LocusCompare2_visualization/Simulate//new_sim_gwaseqtlresult_20240511.tsv',sep='\t')
df$smrqvalue = qvalue(df$SMR.p, fdr.level=0.05, pfdr=F)$qvalues
df$twasqvalue = qvalue(df$TWAS.P, fdr.level=0.05, pfdr=F)$qvalues
df$predixcanqvalue = qvalue(df$PrediXcan.P, fdr.level=0.05, pfdr=F)$qvalues
write.table(df, 'LocusCompare2_visualization/Simulate//new_sim_gwaseqtlresult_20240511_qvalue.tsv', sep='\t', quote=F, row.names=F)


df = read.csv('LocusCompare2_visualization/Simulate//new_sim_gwaseqtlresult_dropna_20240511.tsv',sep='\t')
df$smrqvalue = qvalue(df$SMR.p, fdr.level=0.05, pfdr=F)$qvalues
df$twasqvalue = qvalue(df$TWAS.P, fdr.level=0.05, pfdr=F)$qvalues
df$predixcanqvalue = qvalue(df$PrediXcan.P, fdr.level=0.05, pfdr=F)$qvalues
write.table(df, 'LocusCompare2_visualization/Simulate//new_sim_gwaseqtlresult_dropna_20240511_qvalue.tsv', sep='\t', quote=F, row.names=F)

df = read.csv('LocusCompare2_visualization/Simulate//newtwas2_merged.tsv',sep='\t')
df$qvalue = qvalue(df$TWAS.P, fdr.level=0.05, pfdr=F)$qvalues
df$PROB = 1-df$qvalue
write.table(df, 'LocusCompare2_visualization/Simulate//newtwas2_merged_qvalue.tsv', sep='\t', quote=F, row.names=F)

df = read.csv('/Users/phoebel/github/locuscompare2/new_sim_result_20240504/simulate_result_v1/newsmr_merged.tsv',sep='\t')
df$qvalue = qvalue(df$p_SMR, fdr.level=0.05, pfdr=F)$qvalues
df$PROB = 1-df$qvalue
write.table(df, 'LocusCompare2_visualization/Simulate//newsmr_merged_qvalue.tsv', sep='\t', quote=F, row.names=F)


df = read.csv('LocusCompare2_visualization/Simulate//newpredixcan_merged.tsv',sep='\t')
df$qvalue = qvalue(df$pvalue, fdr.level=0.05, pfdr=F)$qvalues
df$PROB = 1-df$qvalue
write.table(df, 'LocusCompare2_visualization/Simulate//newpredixcan_merged_qvalue.tsv', sep='\t', quote=F, row.names=F)
