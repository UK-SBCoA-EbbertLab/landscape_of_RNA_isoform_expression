library(tidyverse)

hk_iso <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_input_housekeeping_genes.tsv')
p_iso <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_input_preferential_genes.tsv')


for (thresh in c(0.01, 1.01, 5.01, 10.01)) {
	hk_and_p <- read_tsv(paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_threshold_', thresh, '_housekeeping_and_preferential.tsv'))
	#hk_and_p <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_threshold_1.01_housekeeping.tsv')

	hk_deseq2 <- hk_iso %>%
		select(gene_id, transcript_id, gene_name, gene_biotype) %>%
		unite('TXNAME', c(transcript_id, gene_id, gene_name)) %>%
		left_join(hk_and_p) %>%
		distinct() %>%
		filter(pvalue >= 0.05) %>%
		group_by(TXNAME) %>%
		mutate(n_not_sig = n()) %>%
		filter(n_not_sig > 30)

	write_tsv(hk_deseq2, paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_threshold_', thresh, '_housekeeping_filtered.tsv'))


	p_deseq2 <- p_iso %>%
		select(gene_id, transcript_id, gene_name, gene_biotype) %>%
		unite('TXNAME', c(transcript_id, gene_id, gene_name)) %>%
		left_join(hk_and_p) %>%
		distinct() %>%
		filter(pvalue < 0.05) %>%
		group_by(TXNAME) %>%
		mutate(n_sig = n()) %>%
		filter(n_sig == 8 | n_sig == 14 | n_sig == 15 | n_sig == 18 | n_sig == 20 | n_sig == 21)

	write_tsv(p_deseq2, paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_threshold_', thresh, '_preferential_filtered.tsv'))
}


