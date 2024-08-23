nohup python3 -u /home/users/nus/e1124850/locuscompare2/locuscompare2-standalone/predictive-model-pipeline/pipeline.py  \
--rscript_path /home/users/nus/e1124850/.conda/envs/colotools/bin/FUSION.compute_weights.R \   # --hsq_p 1 ensure more weights can be generated
--genotype /home/users/nus/e1124850/e1124850/locuscompare2file/sim_eqtl_20240430/out.vcf.gz \
--expression /home/users/nus/e1124850/e1124850/locuscompare2file/sim_eqtl_20240430/out.bed.gz \
--tissue_name test_tissue \
--output_dir /data/projects/11003054/e1124850/locuscompare2file/sim_coloc/sim_twas_model_6  \
--additional_compute_weight_params "--PATH_plink ~/.conda/envs/colotools/bin/plink  --PATH_gcta ~/.conda/envs/colotools/bin/gcta --PATH_gemma ~/download/GEMMA/bin/gemma"