import pandas as pd
import statsmodels.api as sm
fusion = [
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004599/fusion_output_20240518223753.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004600/fusion_output_20240518230230.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004601/fusion_output_20240518225253.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004602/fusion_output_20240518224158.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004603/fusion_output_20240518224537.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004604/fusion_output_20240518224109.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004605/fusion_output_20240518224349.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004606/fusion_output_20240518225046.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004607/fusion_output_20240518224418.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004608/fusion_output_20240518224947.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004609/fusion_output_20240518225005.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004610/fusion_output_20240518225145.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004611/fusion_output_20240518230947.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004612/fusion_output_20240518230028.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004613/fusion_output_20240518230348.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004614/fusion_output_20240519060804.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004615/fusion_output_20240519074436.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004616/fusion_output_20240519115054.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004617/fusion_output_20240519121959.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004618/fusion_output_20240519141043.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004619/fusion_output_20240519150540.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004620/fusion_output_20240519151238.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004621/fusion_output_20240519152608.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004622/fusion_output_20240519154413.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004623/fusion_output_20240519155421.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004624/fusion_output_20240519155850.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004625/fusion_output_20240519173640.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004626/fusion_output_20240519183310.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004627/fusion_output_20240519220714.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004628/fusion_output_20240519221213.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004629/fusion_output_20240519221036.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004630/fusion_output_20240519221320.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004631/fusion_output_20240520025622.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004632/fusion_output_20240520072336.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004633/fusion_output_20240520073034.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004634/fusion_output_20240520080813.tsv.gz'
    ]

smr = [
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004599/smr_output_20240519201545.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004600/smr_output_20240519151958.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004601/smr_output_20240519145557.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004602/smr_output_20240519190131.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004603/smr_output_20240519195954.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004604/smr_output_20240519072233.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004605/smr_output_20240519054330.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004606/smr_output_20240519152701.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004607/smr_output_20240519170314.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004608/smr_output_20240519134422.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004609/smr_output_20240519144800.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004610/smr_output_20240519115546.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004611/smr_output_20240519152001.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004612/smr_output_20240519142724.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004613/smr_output_20240519112701.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004614/smr_output_20240519180709.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004615/smr_output_20240519185815.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004616/smr_output_20240520093207.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004617/smr_output_20240520070136.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004618/smr_output_20240520022712.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004619/smr_output_20240520095717.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004620/smr_output_20240520111218.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004621/smr_output_20240520134937.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004622/smr_output_20240520103717.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004623/smr_output_20240520073839.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004624/smr_output_20240520092253.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004625/smr_output_20240520102010.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004626/smr_output_20240520115341.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004627/smr_output_20240520100541.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004628/smr_output_20240520074824.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004629/smr_output_20240520094420.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004630/smr_output_20240520115627.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004631/smr_output_20240520065851.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004632/smr_output_20240520232414.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004633/smr_output_20240520210616.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004634/smr_output_20240520112359.tsv.gz'
]

predixcan = [
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004599/predixcan_output_20240518224308.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004600/predixcan_output_20240518231141.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004601/predixcan_output_20240518225917.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004602/predixcan_output_20240518224727.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004603/predixcan_output_20240518225045.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004604/predixcan_output_20240518224636.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004605/predixcan_output_20240518224915.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004606/predixcan_output_20240518225703.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004607/predixcan_output_20240518224934.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004608/predixcan_output_20240518225503.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004609/predixcan_output_20240518225512.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004610/predixcan_output_20240518225628.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004611/predixcan_output_20240518231524.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004612/predixcan_output_20240518230529.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004613/predixcan_output_20240518230855.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004614/predixcan_output_20240519061357.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004615/predixcan_output_20240519075010.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004616/predixcan_output_20240519115617.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004617/predixcan_output_20240519122512.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004618/predixcan_output_20240519141733.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004619/predixcan_output_20240519151149.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004620/predixcan_output_20240519151727.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004621/predixcan_output_20240519153215.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004622/predixcan_output_20240519154948.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004623/predixcan_output_20240519160124.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004624/predixcan_output_20240519160500.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004625/predixcan_output_20240519174226.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004626/predixcan_output_20240519183828.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004627/predixcan_output_20240519193534.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004628/predixcan_output_20240519194501.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004629/predixcan_output_20240519221541.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004630/predixcan_output_20240519221818.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004631/predixcan_output_20240520030124.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004632/predixcan_output_20240520072846.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004633/predixcan_output_20240520073550.tsv.gz',
    '/home/liulab/datasource/Astle_blood_trait/Astle_blood_trait_results/GCST004634/predixcan_output_20240520081311.tsv.gz'
]


total_fusion = pd.DataFrame()
total_predixcan = pd.DataFrame()
total_smr = pd.DataFrame()


for i in fusion:
    df = pd.read_csv(i,sep='\t')
    total_fusion = pd.concat([total_fusion,df],axis=0)

total_fusion = total_fusion.sort_values('TWAS.P') # 576360
total_fusion['fdr'] = sm.stats.multipletests(total_fusion['TWAS.P'],alpha=0.05,method='fdr_bh')[1] # 2.880000e-03
total_fusion['bonferroni'] = sm.stats.multipletests(total_fusion['TWAS.P'],alpha=0.05,method='bonferroni')[1] # 8.670000e-08  -log10=7.062
total_fusion[total_fusion['fdr']<0.05] # 33239
total_fusion[total_fusion['bonferroni']<0.05] # 6366


for i in smr:
    df = pd.read_csv(i,sep='\t')
    total_smr = pd.concat([total_smr,df],axis=0)

total_smr = total_smr.sort_values('p_SMR') # 164023
total_smr['p_SMR'] = total_smr['p_SMR'].fillna(1)
total_smr['fdr'] = sm.stats.multipletests(total_smr['p_SMR'],alpha=0.05,method='fdr_bh')[1] # 4.711000e-03
total_smr['bonferroni'] = sm.stats.multipletests(total_smr['p_SMR'],alpha=0.05,method='bonferroni')[1] # 3.047000e-07  -log10=6.516
total_smr[total_smr['fdr']<0.05] # 15455
total_smr[total_smr['bonferroni']<0.05] # 2536


for i in predixcan:
    df = pd.read_csv(i,sep='\t')
    total_predixcan = pd.concat([total_predixcan,df],axis=0)

total_predixcan = total_predixcan.sort_values('pvalue') # 409788
total_predixcan['fdr'] = sm.stats.multipletests(total_predixcan['pvalue'],alpha=0.05,method='fdr_bh')[1] # 2.831000e-03
total_predixcan['bonferroni'] = sm.stats.multipletests(total_predixcan['pvalue'],alpha=0.05,method='bonferroni')[1] # 1.215000e-07  -log10=6.915
total_predixcan[total_predixcan['fdr']<0.05] # 23208
total_predixcan[total_predixcan['bonferroni']<0.05] # 4873
