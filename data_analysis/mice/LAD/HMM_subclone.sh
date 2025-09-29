#!/bin/bash
#SBATCH --job-name=run_hmm_subclone     # Job name
#SBATCH --time=02:00:00
#SBATCH -c 8
#SBATCH --mem 16gb
#SBATCH --array=0-21%20   # Job array range, limit to 5 simultaneous jobs


source activate cancer_architecture
# Path to the directory containing .hmm files
input_dir="/workspace/projects/lesion_segregation/mice/results/HMM/input_HMM"
# Get the list of .hmm files into an array
files=("$input_dir"/*.hmm)

# Get the file corresponding to the current job index ($SLURM_ARRAY_TASK_ID)
file=${files[$SLURM_ARRAY_TASK_ID]}

# Extract filename without path and extension
filename=$(basename "$file")

echo "Processing $filename"

# Run the R script and pass the filename as an argument
Rscript HMM_subclone.R "$filename"