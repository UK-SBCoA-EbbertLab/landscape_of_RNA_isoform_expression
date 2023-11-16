library(tidyverse)
library(biomaRt)

args = commandArgs(trailingOnly=TRUE)

isoforms_in_tissue <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv')
isoforms_in_tissue_rel_abund <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_median_rel_abund_gt_10.tsv')

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

n_tx_hgnc <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype', 'ensembl_transcript_id', 'transcript_biotype'),
	      filters = 'ensembl_transcript_id',
	      values = unique(isoforms_in_tissue$transcript_id),
	      mart=ensembl)

isoforms_in_tissue <- isoforms_in_tissue %>%
	left_join(n_tx_hgnc %>% rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id))

isoforms_in_tissue_rel_abund <- isoforms_in_tissue_rel_abund %>%
	left_join(n_tx_hgnc %>% rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id))

if(length(args) != 0){
	write_tsv(isoforms_in_tissue %>% filter(hgnc_symbol == args[1]), 
		  paste0('../../tables/GTEx_expression/GTEx_', args[1], '_isoforms_CPM.tsv'))
	write_tsv(isoforms_in_tissue_rel_abund %>% filter(hgnc_symbol == args[1]), 
		  paste0('../../tables/GTEx_expression/GTEx_', args[1], '_isoforms_rel_abund.tsv'))
}
