#!/usr/bin/env bash


nextflow ../main.nf --bambu_rds "./results/GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109.gtf" \
    --out_dir "./GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/" \
    --is_discovery "False" \
    --bambu_track_reads "False" \
    --multiqc_input "./results/GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --fai "./results/GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/fai/*.fai" \
    --step 3 \
    --is_chm13 "False" -resume
