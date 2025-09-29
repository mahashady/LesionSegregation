#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH -c 8
#SBATCH --mem 16gb


source activate cancer_architecture

python annotate_gene_expression2mut.py
