library(readr)
library(tidyr)
library(dplyr)
library(stringr)

GTEx_samples_to_keep <- read_csv("../data/GTExSamplesToKeep.csv")

GTEx_metadata <- read_tsv("../../../data/GTEx_v9_ONT_metadata.txt") %>%
	rename(condition = tissue_site_detail, sampleID = bam_file) %>%
	mutate(sampleID = str_remove(sampleID, ".bam$")) %>%
	select(c(sampleID, condition)) %>%
	#select(c(sampleID, condition, mrna_rin)) %>%
	#select(c(sampleID, condition, date_of_sequencing, mrna_rin)) %>%
	right_join(GTEx_samples_to_keep)

write_tsv(GTEx_metadata, "../data/GTExDesignMatrix.tsv")



