library(tidyverse)


n_tx_per_gene <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_tx_per_gene_passing_thresholds_2023.tsv') %>%
	pivot_longer(!c('gene_id', 'gene_biotype', 'gene_name'), names_to = 'tissue', values_to = 'n_tx') %>%
	drop_na(n_tx) %>%
	group_by(tissue) %>%
	summarize(gtoet_5=sum(n_tx >= 5), gtoet_6=sum(n_tx >= 6) , gtoet_10=sum(n_tx >= 10), gtoet_30=sum(n_tx >= 30))
n_tx_per_pc_gene <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_protein_coding_tx_per_gene_passing_thresholds_2023.tsv') %>%
	pivot_longer(!c('gene_id', 'gene_biotype', 'gene_name'), names_to = 'tissue', values_to = 'n_tx') %>%
	drop_na(n_tx) %>%
	group_by(tissue) %>%
	summarize(pcgtoet_5=sum(n_tx >= 5) , pcgtoet_10=sum(n_tx >= 10))

combined <- inner_join(n_tx_per_gene, n_tx_per_pc_gene)

write_tsv(combined, '../../tables/GTEx_expression/GTEx_number_of_tx_per_gene_by_tissue.tsv')


