#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH -c 8
#SBATCH --mem 16gb


source activate cancer_architecture

#Rscript --vanilla fit_2betas_by_driver_1div_4cells_only.R 7:145859242_T/A
#Rscript --vanilla fit_2betas_by_driver_1div_4cells_only.R 7:145859242_T/C
#Rscript --vanilla fit_2betas_by_driver_1div_4cells_only.R 6:37548568_A/T
Rscript --vanilla fit_2betas_by_driver_1div_4cells_only.R ${SLURM_ARRAY_TASK_ID}
