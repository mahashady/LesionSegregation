#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH -c 8
#SBATCH --mem 16gb


source activate cancer_architecture

Rscript fit_2betas_by_multiHMM_state.R
