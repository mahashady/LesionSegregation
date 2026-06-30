#no recombination
python run_calibration_binary_method.py \
  --tumor_file /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/data/multi_stat_by_chr.enriched.platinum.SBS17b_prop_lt_0.1.tsv \
  --likelihood_file /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/likelihood_matrices/likelihood_matrix_for_LAD_inference.no_recombination.M2-15.LAD1-10.n1000.tsv \
  --outdir ../../results/LAD_inference/calibration/no_recombination/ \
  --q_values 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
  --n_reps 20
#one recombination
#python run_calibration_binary_method.py \
#  --tumor_file /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/data/multi_stat_by_chr.enriched.platinum.SBS17b_prop_lt_0.1.tsv \
#  --likelihood_file /home/bbg/mandrianova/Burst_kinetics/data_analysis/data_analysis/metastatic_tumors/results/LAD_inference/likelihood_matrices/likelihood_matrix_for_LAD_inference.recombination.M2-15.LAD1-10.lambda1.0.win10000000.n1000.tsv \
#  --outdir ../../results/LAD_inference/calibration/1_recombination/ \
#  --q_values 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
#  --n_reps 20

  