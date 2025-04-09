#!/bin/bash
#SBATCH --time=05:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=filter_for_genes      # Job name
#SBATCH --ntasks=1                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr         # Project allocation account name (REQUIRED)

singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif Rscript 02a_filter_for_gene_list.R z_data/counts_transcript.txt z_data/CPM_transcript.txt z_data/GTExDesignMatrix.tsv z_references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109.gtf z_data/genes_for_domain_analysis.txt

singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif bash 02b_filter_fasta_for_gene_list.sh z_data/CPM_transcript_FILTERED.txt z_references/transcriptome.fa z_references/transcriptome_FILTERED.fa
