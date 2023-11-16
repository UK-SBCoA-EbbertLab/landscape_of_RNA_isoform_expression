library(tidyverse)
library(biomaRt)

tib <- read_tsv('../../tables/n_tx_per_gene_2023_over_60.tsv')

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

n_tx_hgnc <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype'),
	      filters = 'ensembl_gene_id',
	      values = unique(tib$gene_id),
	      mart=ensembl)
tib <- tib %>% left_join(n_tx_hgnc %>% mutate(gene_id = ensembl_gene_id, hgnc_symbol = gsub("\"", "", hgnc_symbol)))

print(tib)

write_tsv(tib, '../../tables/n_tx_per_gene_2023_over_60_with_symbol.tsv')

