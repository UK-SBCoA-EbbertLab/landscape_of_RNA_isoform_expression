#!/bin/bash
#SBATCH --time=05:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=my_test_job      # Job name
#SBATCH --ntasks=1                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr         # Project allocation account name (REQUIRED)

singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif Rscript 03_format_gtf_for_ISAR.R z_references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109_FILTERED.gtf z_references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109_FILTERED_UPDATED_FORMAT.gtf z_references/annotation_key_info_FILTERED.tsv
