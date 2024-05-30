library(tidyverse)

cpm_vs_n_iso <- read_tsv('../../tables/GTEx_expression/GTEx_n_iso_diff_vs_gene_CPM.tsv') %>%
        mutate(log_med_gene_cpm = log2(median_gene_cpm + 1))

df <- read_tsv('../../tables/GTEx_expression/GTEx_any_n_tx_by_tissue_gt_5.tsv') %>%
        dplyr::select(!c(min_n_tx, max_n_tx, gene_biotype, gene_name)) %>%
        pivot_longer(-gene_id, names_to = 'tissue', values_to = 'n_tx') %>%
        left_join(cpm_vs_n_iso %>% dplyr::select(gene_id, tissue, median_gene_cpm, log_med_gene_cpm)) %>%
	mutate(n_tx = as.numeric(n_tx)) %>%
	group_by(tissue) %>%
	nest() %>%
	mutate(correlation = map(data, ~cor.test(.$log_med_gene_cpm, .$n_tx, method='spearman'))) %>%
	mutate(correlation = map(correlation, broom::tidy)) %>%
	unnest(correlation)

dplyr::last_dplyr_warnings()

write_tsv(df, '../../tables/GTEx_expression/GTEx_spearman_correlation.tsv')





