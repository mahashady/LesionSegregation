#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH -c 8
#SBATCH --mem 16gb


source activate cancer_architecture

Rscript --vanilla fit_2betas_1div_4cells_only.R 
