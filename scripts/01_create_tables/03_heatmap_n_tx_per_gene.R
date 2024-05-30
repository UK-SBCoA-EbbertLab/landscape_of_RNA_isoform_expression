library(tidyverse)
library(biomaRt)

# load in the number of isoforms per gene in the gtex data, then calculate the min and max n_tx across tissues per gene
n_tx_per_gene <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_tx_per_gene_passing_thresholds_2023.tsv') %>%
	replace(is.na(.), 0) %>%
	rowwise() %>%
	mutate(min_n_tx=min(c_across(!c(gene_id, gene_biotype, gene_name)))) %>%
	mutate(max_n_tx=max(c_across(!c(gene_id, gene_biotype, gene_name))))

# create a file where ALL tisses express more than 5 isoforms for a gene
n_tx_per_gene_table <- n_tx_per_gene %>%
	filter(min_n_tx > 5)
write_tsv(n_tx_per_gene_table, '../../tables/GTEx_expression/GTEx_all_n_tx_by_tissue_gt_5.tsv')

# create a file where at least one tissue expresses more than 5 isoforms for a gene
n_tx_per_gene_by_tiss <- n_tx_per_gene %>%
	filter(max_n_tx > 5)
write_tsv(n_tx_per_gene_by_tiss, '../../tables/GTEx_expression/GTEx_any_n_tx_by_tissue_gt_5.tsv')




## Load the relative abundance, get the same genes that are in the above table, then write to file
#n_tx_per_gene_rel_abund_by_tiss <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_tx_per_gene_passing_median_rel_abund_gt_10_threshold_2023.tsv') %>%
#	replace(is.na(.), 0) %>%
#write_tsv(n_tx_per_gene_rel_abund_by_tiss, '../../tables/GTEx_expression/GTEx_number_of_tx_by_tissue_GTEX_rel_abund_gt_10_match_any_gt_5.tsv')

#############################################
### Load the relative abundance, get the same genes that are in the above table, then write to file
#n_tx_per_gene_rel_abund_by_tiss <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_tx_per_gene_passing_median_rel_abund_0_threshold_2023.tsv') %>%
#	replace(is.na(.), 0)
#write_tsv(n_tx_per_gene_rel_abund_by_tiss, '../../tables/GTEx_expression/GTEx_number_of_tx_by_tissue_GTEX_rel_abund_0_match_any_gt_5.tsv')

