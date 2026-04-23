python run_calibration_binary_method.py \
  --tumor_file ../../results/LAD_inference/data/multi_stat_by_chr.enriched.platinum.SBS17b_prop_lt_0.1.tsv \
  --likelihood_file ../../results/LAD_inference/likelihood_matrices/likelihood_matrix_for_LAD_inference.tsv \
  --outdir ../../results/LAD_inference/calibration/no_recombination/ \
  --q_values 0.1,0.25,0.5,0.75,0.9 \
  --n_reps 20