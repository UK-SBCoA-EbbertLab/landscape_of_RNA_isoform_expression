library(tidyverse)

# This script plots various scatter plots of median gene cpm vs number of isoforms for various conditions

h = 10
w = 10

df <- read_tsv("../../tables/GTEx_expression/GTEx_spearman_correlation_with_gene_length.tsv")


ggplot(df, aes(y=iso_resid, x=expr_resid))+
	geom_jitter(alpha=0.1, width = 0.3, height = 0.3) +
	geom_smooth() +
	#geom_smooth(method = "lm") +
	facet_wrap(~tissue) + 
	geom_text(data = df %>% select(tissue, rho, p) %>% distinct(),
            aes(x = -Inf, y = Inf, label = paste0("rho = ", round(rho, 2), ", p ", p)),
            hjust = -0.1, vjust = 1.5, size = 3.5, inherit.aes = FALSE) +
	labs(
	     title = "Partial Correlation Residuals (Expression vs Isoforms, Gene Length Removed)",
	     x = "Expression Residual (Rank | Gene Length Removed)",
	     y = "Isoform Count Residual (Rank | Gene Length Removed)"
	     ) +
	theme_minimal() +
	theme(text = element_text(size = 15))


ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/GTEx_n_iso_observed_at_CPM_threshold_gt_5_in_at_least_one_tissue_vs_gene_log_CPM_with_gene_length.pdf', width=w, height=h)


ggplot(df %>% filter(tissue == 'Heart - Left Ventricle'), aes(y=iso_resid, x=expr_resid))+
	geom_jitter(alpha=0.1, width = 0.3, height = 0.3) +
	geom_smooth() +
	#geom_smooth(method = "lm") +
	geom_text(data = df %>% filter(tissue == 'Heart - Left Ventricle') %>% select(tissue, rho, p) %>% distinct(),
            aes(x = -Inf, y = Inf, label = paste0("rho = ", round(rho, 2), ", p ", p)),
            hjust = -0.1, vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
	labs(
	     title = "Heart - Left Ventricle Partial Correlation Residuals (Expression vs Isoforms, Gene Length Removed)",
	     x = "Expression Residual (Rank | Gene Length Removed)",
	     y = "Isoform Count Residual (Rank | Gene Length Removed)"
	     ) +
	theme_minimal() +
	theme(text = element_text(size=7))

ggsave('../../figures/M3_GTEx_n_isoforms_heatmap/Heart_Left_Ventricle_scatter_plot_with_gene_length.pdf', width=58, height=58, units='mm')
