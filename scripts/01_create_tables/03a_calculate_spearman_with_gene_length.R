library(tidyverse)

cpm_vs_n_iso <- read_tsv('../../tables/GTEx_expression/GTEx_n_iso_diff_vs_gene_CPM.tsv') %>%
        mutate(log_med_gene_cpm = log2(median_gene_cpm + 1))

gene_length <- read_tsv('../../nextflow_pipeline/references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109.gtf', comment="#", col_names=FALSE) %>%
	filter(X3 == 'gene') %>%
	mutate(gene_id = str_match(X9, 'gene_id "([^"]+)"')[,2]) %>%
	mutate(gene_length = X5 - X4) %>%
	dplyr::select(gene_id, gene_length)

df <- read_tsv('../../tables/GTEx_expression/GTEx_any_n_tx_by_tissue_gt_5.tsv') %>%
        dplyr::select(!c(min_n_tx, max_n_tx, gene_biotype, gene_name)) %>%
        pivot_longer(-gene_id, names_to = 'tissue', values_to = 'n_tx') %>%
        left_join(cpm_vs_n_iso %>% dplyr::select(gene_id, tissue, median_gene_cpm, log_med_gene_cpm)) %>%
	mutate(n_tx = as.numeric(n_tx)) %>%
	left_join(gene_length) %>%
	group_by(tissue) %>%
	nest() %>%
	mutate(data = map(data, ~{

					 dat <- .

					 # Rank transform
					 dat <- dat %>%
						 filter(!is.na(log_med_gene_cpm)) %>%
						 mutate(rank_expr = rank(log_med_gene_cpm),
							rank_iso = rank(n_tx),
							rank_length = rank(gene_length))

					 dat <- dat %>% mutate(expr_resid = residuals(lm(rank_expr ~ rank_length)),
							       iso_resid = residuals(lm(rank_iso ~ rank_length)))

					 cor_res <- cor.test(dat$expr_resid, dat$iso_resid, method = "spearman")

					 dat$rho <- cor_res$estimate
					 dat$p <- format.pval(cor_res$p.value, eps = .Machine$double.eps)

					 dat
	})) %>% unnest(data)


dplyr::last_dplyr_warnings()

write_tsv(df, '../../tables/GTEx_expression/GTEx_spearman_correlation_with_gene_length.tsv')


# Plot with rho and p-value in facet
#ggplot(df, aes(x = expr_resid, y = iso_resid)) +
#  geom_jitter(alpha = 0.1, width = 0.3, height = 0.3) +
#  geom_smooth(method = "lm", color = "red") +
#  facet_wrap(~tissue, scales = "free") +
#  geom_text(data = df %>% distinct(tissue, rho, p),
#            aes(x = -Inf, y = Inf, label = paste0("rho = ", round(rho, 2), ", p ", p)),
#            hjust = -0.1, vjust = 1.5, size = 2.5, inherit.aes = FALSE) +
#  labs(
#    title = "Partial Correlation Residuals (Expression vs Isoforms, Gene Length Removed)",
#    x = "Expression Residual (Rank | Gene Length Removed)",
#    y = "Isoform Count Residual (Rank | Gene Length Removed)"
#  ) +
#  theme_minimal() +
#  theme(text = element_text(size = 5))





