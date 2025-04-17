# Singularity Containers

#### Note: It will be better to pull the files instead of trying to rebuild them, as rebuilding them might result in slightly different software versions



2023-07-13_bambu.def - definition file for singularity container used to run the Bambu R package on the Nextflow Pipeline.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bambu:sha256.152511f143e497212d2b711dbb94eba775270ac4ee196fbce12dc245fcac404f`


guppy.def - definition file for singularity container used to run guppy basecaller... Was not used due to data being basecalled on the PromethION.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/guppy:sha256.80d73ed421e17199390658662453c2c37c8e64435d8b402822a107247882915f`
 


2023-07-13_nanopore.def - definition file for singularity container used to run the software nanopore data analysis.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/nanopore:sha256.1e78cf38407f51f1742a996adcb6bed9019e103ffd830fa210fa0194014b5064`



2023-07-13_quality_control.def - definition file for singularity container used to run the quality control software.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/quality_control:sha256.ee1d1119b5e0dda74da7d4896c2577302bd8066f4cccefc067a3253f4127c92a`



2023-07-13_bernardo_article_analysis.def - definition file for singularity container used to make some tables

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bernardo_article_analysis:sha256.bdd920de92d66d33b8a51398a6f6a5dbdb81c10accc8b841e110edd0a1730486`



pathway_analysis_2023_10_30_add_gini.def - definition file for singularity container used for tables and figures

`pull command: singularity pull --arch amd64 library://ebbertlab/rna_landscape/pathway_analysis_2023_10_30:sha256.b4fd50793609467e776b760549a8fb8f403aa8d4963c7606be7d158fd67d8330`



ISAR_software_MAR_11_2025.def - definition file for singularity container used to do analyses with IsoformSwitchAnalyzeR

`pull command: singularity pull --arch amd64 library://ebbertlab/isoform_switch_analyze_r/isar:sha256.08e4ceaa0d28cd0e96382da3ad3f7873050fc3ea08d56307c99f717ba98b2496`

