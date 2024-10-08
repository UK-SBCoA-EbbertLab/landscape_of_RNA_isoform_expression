# landscape_of_RNA_isoform_expression
## new_annotated_transcripts_per_year


### ENSEMBL annotations

Instructions on how to retrieve the ENSEMBL annotations and related documentation can be found under: `./references/`

We took ten representative ENSEMBL annotations from 2014 till 2023 (one per year) and calculate the number of genes and isoforms annotated in each, as well as conducted some other analyses. See the code and our paper for more information. The code for this is under `scripts/01_create_tables/` and `scripts/02_create_plots/`, the related scripts begin with `01` or `plot_M1`. 

### Long-read sequencing expression data

Long-read sequencing data from 15 tissues were obtained from GTEx. Fastq files were downloaded and re-analyzed using our NextFlow pipeline. For more details see the original publication with this data: https://pubmed.ncbi.nlm.nih.gov/35922509/

Deep long-read brain sequencing data were obtained from a previous experiment in the lab. We took the 700 Details about the data can be found at: https://github.com/UK-SBCoA-EbbertLab/brain_cDNA_discovery

### Nextflow Pipeline 

Fastq data were analyzed using our in-house Nanopore cDNA/RNA sequencing Nextflow pipeline. Documentation can be found under: `./nextflow_pipeline/`

### Singularity containers:

Information on how to retrieve the singularity containers with the software used in this project, as well as definition files and related documentation can be found under `./singularity/`


### Code used to generate figures and metrics

Scripts containing code used to generate plots and metrics in this project as well as related documentation can be found under `./scripts/`


### Figures

Figures generated for this project can be found under `./figures/`


### Tables 

Tables generated for this project can be found under `./tables/`

