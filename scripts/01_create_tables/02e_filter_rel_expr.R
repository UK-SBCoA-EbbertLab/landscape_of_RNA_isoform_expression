library(tidyverse)

# load the isoforms passing our thresholds, pivot longer, and drop any rows where the cpm is NA
cpm_isoform <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv') %>%
	select(-c(gene_biotype, threshold)) %>%
	pivot_longer(!c(transcript_id, gene_id, gene_name), names_to='tissue', values_to = 'med_isoform_cpm') %>%
	drop_na(med_isoform_cpm)

rel_abund <- read_tsv('../../tables/GTEx_expression/GTEx_tx_relative_abundance.tsv') %>%
	pivot_longer(!c(transcript_id, gene_id, gene_name), names_to='tissue', values_to = 'rel_abund') %>%
	right_join(cpm_isoform  %>%
		   select(transcript_id, gene_id, gene_name, tissue) %>%
		   distinct()) %>%
	arrange(gene_name, tissue, rel_abund)
write_tsv(rel_abund, '../../tables/GTEx_expression/GTEx_tx_relative_abundance_iso_passing_CPM_gt_1.tsv')

