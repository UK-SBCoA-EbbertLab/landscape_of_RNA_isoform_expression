#!/bin/bash
#SBATCH --time=05:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=my_test_job      # Job name
#SBATCH --ntasks=1                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr         # Project allocation account name (REQUIRED)

singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif Rscript 31a_ISAR_part2.R z_output/ISAR_part1_object.RData z_output/cpat3_output.formatted.tsv z_output/pfam_output.txt z_output/combined_iupred2a_output.txt z_output/prediction_results.txt z_output/deeploc2_output_formatted.tsv z_output/biolib_results/TMRs.gff3 
