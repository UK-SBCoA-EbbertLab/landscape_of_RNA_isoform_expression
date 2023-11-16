#!/bin/bash
#SBATCH --time=00:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=M1_annotation_compare      # Job name
#SBATCH --ntasks=1                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr         # Project allocation account name (REQUIRED)

singularity exec /project/mteb223_uksr/singularity_files/pathway_analysis_2023_10_30.sif Rscript plot_M1_annotation_comparison.R
