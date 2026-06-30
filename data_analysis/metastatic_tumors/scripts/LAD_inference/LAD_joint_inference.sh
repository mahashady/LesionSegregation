python LAD_joint_inference.py \
  --tumor_tsv /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/data/multi_stat_by_chr.enriched.platinum.SBS17b_prop_lt_0.1.tsv \
  --likelihood_tsv /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/likelihood_matrices/likelihood_matrix_for_LAD_inference.no_recombination.M2-15.LAD1-10.n1000.tsv \
  --output_dir /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/inference_results/binary/no_recombination.M2-15.LAD1-10.n1000/LOW_1-2-3-4 \
  --prefix platinum.SBS17b_lt_0.1 \
  --mode binary \
  --low_lads 1,2,3,4 \
  --filter_M_not_equal 1 \
  --bootstrap \
  --n_boot 1000 \
  --seed 42


