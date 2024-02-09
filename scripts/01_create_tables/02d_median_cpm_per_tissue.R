library(tidyverse)

# load the isoforms passing our thresholds, pivot longer, and drop any rows where the cpm is NA
cpm_isoform <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv') %>%
	select(-c(gene_biotype, threshold)) %>%
	pivot_longer(!c(transcript_id, gene_id, gene_name), names_to='tissue', values_to = 'med_isoform_cpm') %>%
	drop_na(med_isoform_cpm)

# load in the gene cpm (note this has genes that did not have isoforms meeting our thresholds
# to deal with this, we do a join with the isoform tibble to only get those genes that 
# had at least one isoform that passed our thresholds by tissue
# then get the mean and median gene CPM per tissue
cpm_gene <- read_tsv('../../tables/GTEx_expression/GTEx_median_gene_CPM.tsv') %>%
	pivot_longer(!c(gene_id, gene_name), names_to = 'tissue', values_to = 'med_gene_cpm') %>%
	right_join(cpm_isoform %>%
		   select(gene_id, gene_name, tissue) %>%
		   distinct()) %>%
	drop_na(med_gene_cpm) %>%
	group_by(tissue) %>%
	summarize(median_gene_cpm = median(med_gene_cpm), mean_gene_cpm = mean(med_gene_cpm))

# get the mean and median isoform CPM expression per tissue
cpm_isoform <- cpm_isoform %>%
	group_by(tissue) %>%
	summarize(median_isoform_cpm = median(med_isoform_cpm), mean_isoform_cpm = mean(med_isoform_cpm))

# combine and write
means_and_med_per_tiss <- cpm_isoform %>%
	inner_join(cpm_gene) 

write_tsv(means_and_med_per_tiss, '../../tables/GTEx_expression/GTEx_median_iso_and_gene_CPM_per_tissue.tsv')

