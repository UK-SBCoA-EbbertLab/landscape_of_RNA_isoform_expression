## Year corresponding to each annotation

`Homo_sapiens.GRCh38.76.gtf` = 2014 

`Homo_sapiens.GRCh38.79.gtf` = 2015

`Homo_sapiens.GRCh38.85.gtf` = 2016

`Homo_sapiens.GRCh38.88.gtf` = 2017

`Homo_sapiens.GRCh38.92.gtf` = 2018

`Homo_sapiens.GRCh38.96.gtf` = 2019

`Homo_sapiens.GRCh38.99.gtf` = 2020

`Homo_sapiens.GRCh38.103.gtf` = 2021

`Homo_sapiens.GRCh38.106.gtf` = 2022

`Homo_sapiens.GRCh38.109.gtf` = 2023


## Commands to retrieve annotations:

Annotations were obtained from the ENSEMBL portal: https://ftp.ensembl.org

Commands used to retrieve annotations:

2014 annotation (ENSEMBL 76): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.76.gtf.gz`

2015 annotation (ENSEMBL 79): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.79.gtf.gz`

2016 annotation (ENSEMBL 85): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.85.gtf.gz`

2017 annotation (ENSEMBL 88): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz`

2018 annotation (ENSEMBL 92): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz`

2019 annotation (ENSEMBL 96): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz`

2020 annotation (ENSEMBL 99): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz`

2021 annotation (ENSEMBL 103): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz`

2022 annotation (ENSEMBL 106: `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz`

2023 annotation (ENSEMBL 109): `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz`


## Unzipping annoatations

After retrieving all annotations with the commands above, the files were unzipped with the following command: `gzip -d *.gz`


## Important to note:

Annotations from 2014 and 2015 (Version 76 and 79) include alternative assemblies and haplotypes. To make the comparison fair we only included the chromosomes that were also present in
annotations at or after 2016 (Version 85 and above).



## Other reference files:

`./brain_disease_genes_only_IDs.tsv` - Brain disease associated gene ids and gene symbols. See methods section from the publication associated with this project for a detailed description on how this gene list was generated.
`./brain_disease_genes_with_disease.tsv` - Brain disease associated gene ids and gene symbols, also includes name of disease associated with gene. See methods section from the publication associated with this project for a detailed description on how this gene list was generated.
`./medically_relevant_genes.tsv` - Gene ids and gene symbols for genes associated to disease by GWAS and/or causal variants. See methods section from the publication associated with this project for a detailed description on how this gene list was generated.


