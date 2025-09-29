#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH -c 8
#SBATCH --mem 64gb


source activate cancer_architecture

Rscript reassign_mutations.R
