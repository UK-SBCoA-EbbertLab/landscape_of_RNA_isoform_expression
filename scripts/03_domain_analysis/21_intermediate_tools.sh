#!/bin/bash
#SBATCH --time=05:15:00             # Time limit for the job (REQUIRED).
#SBATCH --job-name=my_test_job      # Job name
#SBATCH --ntasks=1                  # Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=100G                  # Total memory requested
#SBATCH --partition=normal          # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             # Error file for this job.
#SBATCH -o slurm-%j.out             # Output file for this job.
#SBATCH -A coa_mteb223_uksr         # Project allocation account name (REQUIRED)

# run cpat3
#singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif cpat -g z_references/isoformSwitchAnalyzeR_isoform_nt.fasta -o z_output/cpat3_output -d /usr/local/cpat_files/Human_logitModel.RData -x /usr/local/cpat_files/Human_Hexamer.tsv --width 1000
#
#singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif Rscript format_cpat3.R z_output/cpat3_output.ORF_prob.best.tsv 
#
#sed -i '1s/id\t//' z_output/cpat3_output.formatted.tsv
#
## pfam_scan
#rm -f z_output/pfam_output.txt
#singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif /usr/local/PfamScan/pfam_scan.pl -fasta z_references/isoformSwitchAnalyzeR_isoform_AA.fasta -dir /usr/local/PfamScan/ -outfile z_output/pfam_output.txt
#
## signalP6
#singularity exec singlarity_containers/ISAR_software_do_not_distribute_MAR_05_2025.sif signalp6 --fastafile z_references/isoformSwitchAnalyzeR_isoform_AA.fasta --output_dir z_output/ --format none --organism eukarya --mode fast
#
## split fastas
#mkdir -p z_output/split_fastas
#awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {print seq}' z_references/isoformSwitchAnalyzeR_isoform_AA.fasta > z_references/isoformSwitchAnalyzeR_isoform_AA_single_line.fasta
#awk '/^>/ {close(out); split($0, arr, " "); seqname=substr($1,2); out="z_output/split_fastas/" seqname ".fasta"} {print > out}' z_references/isoformSwitchAnalyzeR_isoform_AA_single_line.fasta
#
#rm -f z_output/combined_iupred2a_output.txt
## iupred2a
#for file in z_output/split_fastas/*.fasta; do
#	filename=$(basename "$file")
#	id=${filename%.*}
#	echo -e ">${id}\tNA\tNA\tNA" >> z_output/combined_iupred2a_output.txt
#	singularity exec singlarity_containers/ISAR_software_do_not_distribute_MAR_05_2025.sif python /usr/local/iupred2a/iupred2a.py -a $file long >> z_output/combined_iupred2a_output.txt
#done
#
#sed -i '/^#/d' z_output/combined_iupred2a_output.txt
#
#singularity exec singlarity_containers/ISAR_software_do_not_distribute_MAR_05_2025_deeploc2.sif deeploc2 -f z_references/isoformSwitchAnalyzeR_isoform_AA.fasta -o z_output/ -m Fast -d cpu
#
#singularity exec --env BIOLIB_TOKEN=mYIqWdVdO4fiwcFwfvGI2wuA5TUrnEeFdyfUQcIGTKU singlarity_containers/ISAR_software_MAR_11_2025.sif biolib run 'DTU/DeepTMHMM:1.0.24' --fasta z_references/isoformSwitchAnalyzeR_isoform_AA.fasta

singularity exec singlarity_containers/ISAR_software_MAR_11_2025.sif Rscript format_deeploc2.R z_output/results_20250414-084016.csv 
