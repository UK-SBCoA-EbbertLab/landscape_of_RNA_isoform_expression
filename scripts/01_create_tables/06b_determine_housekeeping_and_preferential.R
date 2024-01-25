library(tidyverse)

# load in list of tissues that we are considering for 'housekeeping' and for 'preferential' 
hk_iso <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/tx_in_nine_gtex_tissues.tsv') %>%
	unite('TXNAME', c(transcript_id, gene_id))
p_iso <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/tx_in_at_least_one_gtex_tissue.tsv') %>%
	unite('TXNAME', c(transcript_id, gene_id))

# load all the deseq2 pairwise comparisons of the new isoforms
hk_and_p <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms.tsv') %>%
	mutate(fdr = p.adjust(pvalue, 'fdr'))

# for three thresholds, calculate the housekeeping isoforms
for (thresh in c(1.01, 5.01, 10.01)) {
	for (lfc in c(1, 2)) {
		hk_deseq2 <- hk_iso %>%
			filter(threshold == thresh) %>%
			select(TXNAME, gene_name, gene_biotype) %>%
			distinct() %>%
			left_join(hk_and_p) %>%
			distinct() %>%
			filter(abs(log2FoldChange) < lfc) %>%
			group_by(TXNAME) %>%
			mutate(n_not_sig = n())
	
		write_tsv(hk_deseq2, paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_thresh_', thresh, '_lfc_', lfc, '_housekeeping_all.tsv'))
		
		hk_deseq2 <- hk_deseq2 %>%
			filter(n_not_sig > 30) %>%
			select(TXNAME, gene_name, gene_biotype, n_not_sig) %>%
			distinct()
	
		write_tsv(hk_deseq2, paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_thresh_', thresh, '_lfc_', lfc, '_housekeeping_filtered.tsv'))
	}
}

# prep for preferential isoform analysis, filter the threshold, logFold change and fdr
p_deseq2 <- p_iso %>%
	filter(threshold == 1.01) %>%
	select(TXNAME, gene_name, gene_biotype) %>%
	distinct() %>%
	left_join(hk_and_p) %>%
	distinct() %>%
	filter(abs(log2FoldChange) >= 1) %>%
	filter(fdr < 0.1)

write_tsv(p_deseq2, '../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_1.01_preferential_all.tsv')

# split out isoforms into their own tibbles
list_of_tibs <- p_deseq2 %>%
	mutate(upregulated_tissue = case_when(
					      log2FoldChange > 0 ~ tis1,
					      log2FoldChange < 0 ~ tis2)) %>%
	group_by(TXNAME) %>%
	group_split()

# create a tibble to collect them
preferential_isoforms <- tibble(
				TXNAME = character(),
				n_tis_pref = numeric(),
				tissues = character())

# for each isoform tibble, determine if it has tissues where it is preferentially expressed
for (tib in list_of_tibs) {
	# collect all tissues involved in significant pairwise comparisons and count
	tmp <- c()
	tmp <- append(tmp, tib$tis1)
	tmp <- append(tmp, tib$tis2)
	n_diff <- as_tibble(as.data.frame(table(tmp))) %>%
		filter(Freq > 5)
	# determine which tissues are upregulated compared to other tissues (needs to be up regulated compared to more than 5 tissues)
	up_reg <- as_tibble(as.data.frame(table(tib$upregulated_tissue))) %>%
		filter(Freq > 5)

	# if the isoform is differentially expressed in a tissue compared to at least 6 other tissues
	if (nrow(n_diff) > 0) {
		# A tissue needs to be condistered upregulated and differentially expressed compared to other tissues
		pref_tissues <- up_reg %>%
			rename(tmp = Var1) %>%
			select(tmp) %>%
			left_join(n_diff)

		if (nrow(pref_tissues) > 0) {
			# if it is preferentially expressed in 1,2, or 3 tissues, keep it
			if (nrow(pref_tissues) == 1 && unique(pref_tissues$Freq) == 8) {
				preferential_isoforms <- preferential_isoforms %>% add_row(
											   TXNAME=unique(tib$TXNAME),
											   n_tis_pref=1,
											   tissues = paste0(pref_tissues$tmp, collapse=','))
			} else if ((nrow(pref_tissues) == 2) && max(pref_tissues$Freq) > 6) { #TODO: CHECK THIS <7 was the orig
				preferential_isoforms <- preferential_isoforms %>% add_row(
											   TXNAME=unique(tib$TXNAME),
											   n_tis_pref=2,
											   tissues = paste0(pref_tissues$tmp, collapse=','))
			} else if (nrow(pref_tissues) == 3) {
				preferential_isoforms <- preferential_isoforms %>% add_row(
											   TXNAME=unique(tib$TXNAME),
											   n_tis_pref=3,
											   tissues = paste0(pref_tissues$tmp, collapse=','))
			} else {
				print('SKIPME')
				print(unique(tib$TXNAME))
			}

		}
	}
}
	
print(preferential_isoforms)

preferential_isoforms <- preferential_isoforms %>%
	left_join(p_deseq2 %>% select(TXNAME, gene_name, gene_biotype) %>% distinct())
# write the list of isoforms
write_tsv(preferential_isoforms %>% arrange(n_tis_pref, TXNAME), '../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_1.01_preferential_filtered.tsv')














