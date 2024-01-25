library(tidyverse)

n_tx_per_gene <- read_tsv('../../tables/annotation_comparison/n_tx_per_gene_2023.tsv') %>%
	mutate(percentile = percent_rank(n_tx)*100)
n_tx_per_pcgene <- read_tsv('../../tables/annotation_comparison/n_tx_per_pc_gene_2023.tsv') %>%
	mutate(percentile = percent_rank(n_tx)*100)

print(quantile(n_tx_per_gene$n_tx))
print(quantile(n_tx_per_gene$n_tx, probs = c(0.85, 0.95)))
print(quantile(n_tx_per_pcgene$n_tx))

sink('../../tables/annotation_comparison/n_tx_per_gene_percentiles.txt')
cat("n tx per gene percentiles:\n")
cat("0% 25% 50% 75% 100%\n")
cat(quantile(n_tx_per_gene$n_tx))
cat("\n85% 95%\n")
cat(quantile(n_tx_per_gene$n_tx, probs = c(0.85, 0.95)))
cat("\nn tx per protein coding gene percentiles:\n")
cat("0% 25% 50% 75% 100%\n")
cat(quantile(n_tx_per_pcgene$n_tx))

