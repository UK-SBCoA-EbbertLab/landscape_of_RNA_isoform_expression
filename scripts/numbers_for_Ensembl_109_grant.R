library(tidyverse)

gtf <- read_tsv('../nextflow_pipeline/references/Homo_sapiens.GRCh38.109.gtf', comment='#', col_names=FALSE) %>%	
	mutate(
	       gene_id = str_extract(X9, 'gene_id "([^"]+)"', group=1),
	       gene_name = str_extract(X9, 'gene_name "([^"]+)"', group=1),
	       gene_biotype = str_extract(X9, 'gene_biotype "([^"]+)"', group=1)
	       ) %>%
	mutate(
	       transcript_id = case_when(
					 grepl('transcript_id', X9) ~ str_extract(X9, 'transcript_id "([^"]+)"', group=1), 
					 TRUE ~ NA),
	       transcript_biotype = case_when(
					      grepl('transcript_biotype', X9) ~ str_extract(X9, 'transcript_biotype "([^"]+)"', group=1),
					      TRUE ~ NA)
	       )

n_iso <- nrow(gtf %>%
	filter(X3 == 'transcript'))
print('Number of isoforms')
print(n_iso)

n_gene_body <- nrow(gtf %>%
		    filter(X3 == 'gene'))

print('Number of gene bodies')
print(n_gene_body)

n_pc_genes <- nrow(gtf %>%
		   filter(X3 == 'gene') %>%
		   filter(gene_biotype == 'protein_coding'))

print('Number of pc genes')
print(n_pc_genes)

