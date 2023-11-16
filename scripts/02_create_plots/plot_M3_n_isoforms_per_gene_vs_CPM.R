library(tidyverse)

# This script plots various scatter plots of median gene cpm vs number of isoforms for various conditions

h = 10
w = 10
cpm_vs_n_iso <- read_tsv('../../tables/GTEx_expression/GTEx_n_iso_diff_vs_gene_CPM.tsv') %>%
	mutate(log_med_gene_cpm = log(median_gene_cpm))

# plot the med gene CPM vs expressed isoforms
ggplot(cpm_vs_n_iso, aes(y=n_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() + 
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_vs_gene_log_CPM_all_genes.pdf', width=w, height=h)
	
ggplot(cpm_vs_n_iso %>% filter(gene_biotype == 'protein_coding'), aes(y=n_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() +
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_vs_gene_log_CPM_pc_genes.pdf', width=w, height=h)


# plot med gene CPM vs the difference between annotated and expressed isoforms
ggplot(cpm_vs_n_iso, aes(y=n_tx_diff, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() + 
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_diff_vs_gene_log_CPM_all_genes.pdf', width=w, height=h)
	
ggplot(cpm_vs_n_iso %>% filter(gene_biotype == 'protein_coding'), aes(y=n_tx_diff, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() +
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_diff_vs_gene_log_CPM_pc_genes.pdf', width=w, height=h)

# plot med gene CPM vs the n annotated isoforms
ggplot(cpm_vs_n_iso, aes(y=n_anno_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() + 
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_anno_vs_gene_log_CPM_all_genes.pdf', width=w, height=h)
	
ggplot(cpm_vs_n_iso %>% filter(gene_biotype == 'protein_coding'), aes(y=n_anno_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() +
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_anno_vs_gene_log_CPM_pc_genes.pdf', width=w, height=h)

# plot med gene CPM vs annotated isoforms gt 5
ggplot(cpm_vs_n_iso %>% filter(n_anno_tx > 5), aes(y=n_anno_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() + 
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_anno_gt_5_vs_gene_log_CPM_all_genes.pdf', width=w, height=h)

# plot med gene CPM vs observed isoforms gt 5
ggplot(cpm_vs_n_iso %>% filter(n_tx > 5), aes(y=n_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() +
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_observed_gt_5_vs_gene_log_CPM_all_genes.pdf', width=w, height=h)


## plot median gene CPM vs observed isoforms where median CPM > 1 and at least one tissue is expressing gt 5 isoforms
cpm_threshold_gene_cpm_vs_n_iso <- read_tsv('../../tables/GTEx_expression/GTEx_any_n_tx_by_tissue_gt_5.tsv') %>%
	dplyr::select(!c(min_n_tx, max_n_tx, gene_biotype, gene_name)) %>%
	pivot_longer(-gene_id, names_to = 'tissue', values_to = 'n_tx') %>%
	left_join(cpm_vs_n_iso %>% dplyr::select(gene_id, tissue, log_med_gene_cpm)) 

ggplot(cpm_threshold_gene_cpm_vs_n_iso, aes(y=n_tx, x=log_med_gene_cpm))+
	geom_jitter(alpha=0.1) +
	geom_smooth() +
	facet_wrap(vars(tissue))
ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_observed_at_CPM_threshold_gt_5_in_at_least_one_tissue_vs_gene_log_CPM.pdf', width=w, height=h)

