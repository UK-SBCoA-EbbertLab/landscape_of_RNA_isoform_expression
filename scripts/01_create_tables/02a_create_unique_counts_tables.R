library(tidyverse)

# create tables containing information about number of isoforms per gene using total_unique_counts across all samples

######################################################################################
# load medically relevant genes
mr_genes <- read_tsv('../../references/medically_relevant_genes.tsv') %>%
	mutate(mr = 'mr_gene')

# create tables of the number of isoforms per gene, across all tissue samples, based on unique counts,
# using 4 different thresholds: total unique counts > 1, 5, 10, and 20
unique_n_tx_per_gene <- tibble()

## looking at the expression across all tissues/samples where there was at least N unique counts for the isoform.
unique_counts_n_tx <- read_tsv('../../tables/GTEx_expression/GTEx_total_unique_counts_all_samples_by_tx.tsv')

combined_unique_n_tx <- tibble()

values_at_thresholds_for_paper <- tibble()
# for each threshold, calculate the number of isoforms for each gene and the number of genes that have n isoforms
for (unique_counts_threshold in c(1,5,10,20)){
	# filter by the threshold and get the number of isoforms for each gene
	unique_n_tx <- unique_counts_n_tx %>%
		filter(total_unique_counts >= unique_counts_threshold) %>%
		group_by(GENEID, gene_name, gene_biotype) %>%
		summarise(n_tx = n()) 

	# bind the rows to a table across the thresholds
	unique_n_tx_per_gene <- bind_rows(unique_n_tx_per_gene, unique_n_tx %>%
					  mutate(threshold = unique_counts_threshold))
	tmp_unique_n_tx <- unique_n_tx

	# combine with gene_symbol, gene_biotype, and whether it is a medically relevant gene and then filter
	tmp_unique_n_tx <- tmp_unique_n_tx %>%
		left_join(mr_genes %>% rename(GENEID = gene_id)) %>%
		arrange(n_tx)
	tmp_pc_unique_n_tx <- tmp_unique_n_tx %>% filter(gene_biotype == 'protein_coding')
	tmp_mr_unique_n_tx <- tmp_unique_n_tx %>% filter(mr == 'mr_gene')

	# create table with values used in the paper, describing the isoforms per gene using all samples unique counts
	values_at_thresholds_for_paper <- bind_rows(values_at_thresholds_for_paper, tribble(
		~unique_counts_threshold, ~gene_type, ~at_least_10_isoforms, ~at_least_25_isoforms, ~at_least_60_isoforms,
		unique_counts_threshold, 'all', nrow(tmp_unique_n_tx %>% filter(n_tx >= 10)), nrow(tmp_unique_n_tx %>% filter(n_tx >= 25)), nrow(tmp_unique_n_tx %>% filter(n_tx >= 60)),
		unique_counts_threshold, 'protein-coding', nrow(tmp_pc_unique_n_tx %>% filter(n_tx >= 10)), nrow(tmp_pc_unique_n_tx %>% filter(n_tx >= 25)), nrow(tmp_pc_unique_n_tx %>% filter(n_tx >= 60)),
		unique_counts_threshold, 'medically relevant', nrow(tmp_mr_unique_n_tx %>% filter(n_tx >= 10)), nrow(tmp_mr_unique_n_tx %>% filter(n_tx >= 25)), nrow(tmp_mr_unique_n_tx %>% filter(n_tx >= 60))
		))

	# calculate the number of genes that are expressing the same number of isoforms
	unique_n_tx <- unique_n_tx %>%
		group_by(n_tx) %>%
		summarise(n_genes = n()) %>%
		mutate(threshold = unique_counts_threshold)
	combined_unique_n_tx <- bind_rows(combined_unique_n_tx, unique_n_tx)
}

# file contains information about the number of isoforms at different unique counts thresholds and different gene types
write_tsv(values_at_thresholds_for_paper, '../../tables/GTEx_expression/unique_counts_all_samples_looking_at_at_least_10_25_or_60_isoforms_per_gene.tsv')
# file contains information about the number of genes expressing n isoforms at the different thresholds
write_tsv(combined_unique_n_tx, '../../tables/GTEx_expression/unique_counts_all_samples_across_thresholds_n_tx_and_n_genes.tsv')
# file contains information about the number of isoforms a gene is expressing at different unique counts thresholds -- this one includes the gene ids
write_tsv(unique_n_tx_per_gene %>%
	arrange(threshold, n_tx), '../../tables/GTEx_expression/unique_counts_all_samples_across_thresholds_n_tx.tsv')


