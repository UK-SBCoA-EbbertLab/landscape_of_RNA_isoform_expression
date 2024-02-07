## First pipeline step was submmitted by running:
`. ./submit_pipeline_step_2.sh`

In short, this step starts with basecalled fastq files and:
  1) Adapter trimming and read strand orientation with `pychopper`.
  2) Alignment to the human GRCh38 reference genome using `minimap2`.
  3) Only keeps reads with a mapq score > 10 after alignment using `samtools`.
  4) Prepares Bambu RDS files for quantification using ENSEMBL annotation version 109 (No new transcript discovery) using `bambu`.
  5) Performs QC steps using `pycoqc` and `rseqc`.

## Second pipeline step was submitted by running: 
`. ./submit_pipeline_step_3.sh`

In short, this step does:
  1) Creates a QC report for all files using `multiqc`.
  2) Quantifies transcripts for all files using `bambu`.
  3) Creates a transcriptome fasta file using `gffread`.

## More information:

We ran 12 aged postmortem human dorsolateral frontal cortex (Brodmann area 9/46) brain samples (50% female) through this pipeline. Samples were sequenced using
one Oxford Nanopore PromethION R9.4.1 flow cell per sample. We used kit PCS111 (PCR amplified cDNA sequencing) for library preparation. More detailed information about 
the samples, sequencing protocol, and analysis pipeline can be found in the article publication associated with this project.
