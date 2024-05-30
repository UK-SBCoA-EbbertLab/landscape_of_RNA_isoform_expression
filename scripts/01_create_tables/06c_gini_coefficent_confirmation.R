library(tidyverse)
library(DescTools)

# load in deseq2 normalized matrix
normalized_counts <- read_tsv('deseq2_normalized_values_bambu_iso.tsv')

# load in meta_data
meta_data <- read_tsv('../../data/GTEx_v9_ONT_metadata.txt')

normalized_counts_ginic <- normalized_counts %>%
	rowwise() %>%
	mutate(giniC = Gini(c_across(where(is.numeric))))
#	mutate(giniC_nona = Gini(na.omit(c_across(where(is.numeric)))))


above_cpm_tib <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_0_1_5_10.tsv') %>%
	filter(grepl('Bambu', transcript_id)) %>%
	unite('transcript_id', c('transcript_id','gene_id'))

id_to_gene_name <- above_cpm_tib %>%
	filter(threshold == 0.01) %>%
	select(transcript_id, gene_name)

above_cpm <- above_cpm_tib %>%
	filter(threshold == 1.01) %>%
	pull(transcript_id)

normalized_counts_ginic <- normalized_counts_ginic %>%
	left_join(id_to_gene_name)

gini_gene <- read_tsv('gini_gene_lists_GTEx.txt') %>%
	select(Genes) %>% 
	rename(gene_name = Genes) %>%
	mutate(isGiniGene = TRUE)



write_tsv(normalized_counts_ginic %>% select(transcript_id, giniC), '../../tables/GTEx_expression_our_new_isoforms/deseq_normalized_counts_ginic.tsv')

##########################################################################
######################## Housekeeping ####################################
##########################################################################

housekeeping <- normalized_counts_ginic %>%
	filter(giniC < 0.3) %>%
	select(transcript_id, gene_name, giniC) %>% 
	arrange(giniC) %>%
	left_join(gini_gene)

housekeeping_lax <- normalized_counts_ginic %>%
	filter(giniC < 0.4) %>%
	select(transcript_id, gene_name, giniC) %>% 
	arrange(giniC) %>%
	left_join(gini_gene)


write_tsv(housekeeping, '../../tables/GTEx_expression_our_new_isoforms/housekeeping_ginic.tsv')
write_tsv(housekeeping_lax, '../../tables/GTEx_expression_our_new_isoforms/housekeeping_lax_ginic.tsv')

overlap <- full_join(housekeeping, 
		     read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_thresh_1.01_lfc_2_housekeeping_filtered_with_gini_gene.tsv') %>% rename(transcript_id = TXNAME))

write_tsv(overlap, '../../tables/GTEx_expression_our_new_isoforms/housekeeping_ginic_overlap_with_lfc.tsv')

overlap_lax <- full_join(housekeeping_lax, 
		     read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_thresh_1.01_lfc_2_housekeeping_filtered_with_gini_gene.tsv') %>% rename(transcript_id = TXNAME))

write_tsv(overlap_lax, '../../tables/GTEx_expression_our_new_isoforms/housekeeping_ginic_lax_overlap_with_lfc.tsv')


##########################################################################
######################## Preferential ####################################
##########################################################################
#
#pref <- normalized_counts_ginic %>%
#	filter(giniC > 0.7) %>%
#	arrange(giniC) %>%
#	filter(transcript_id %in% above_cpm)
#	#TODO: want to filter so that at least one of the tissue groups has a median above 1 CPM
#	# I think we will also want to check out the lorenz curve and see which samples are causing the inequality --> check and see that they are all from the same tissues?
#
#pref_cpm <- above_cpm_tib %>%
#	filter(threshold == 0.01) %>%
#	mutate(across(`Brain - Frontal Cortex (BA9)`:Liver, ~replace_na(.x, 0))) %>%
#	filter(transcript_id %in% above_cpm) %>%
#	rowwise() %>%
#	mutate(giniC_cpm = Gini(c_across(where(is.numeric))))
#
#preferential <- left_join(pref, pref_cpm) %>% 
#	filter(giniC_cpm > 0.7)
#
#write_tsv(preferential, 'preferential_3_27_24.tsv')
#write_tsv(preferential %>% select(transcript_id, giniC, giniC_cpm), 'preferential_summary_3_27_24.tsv')
#
#preferential_join <- full_join(
#			       preferential %>% select(transcript_id, giniC, giniC_cpm), 
#			       read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_1.01_preferential_filtered.tsv') %>% rename(transcript_id = TXNAME))
#
#write_tsv(preferential_join, 'preferential_overlap_3_27_24.tsv')
#
#pref_combos <- unique(pref %>% pull(transcript_id))
#
#testing_theory <- normalized_counts %>% filter(transcript_id %in% pref_combos)
#
##TODO: want to calculate all the LC for each row and compare
#
#
#pref_lfc_ginis <- full_join(
#                               left_join(normalized_counts_ginic %>% select(transcript_id, giniC), pref_cpm %>% select(transcript_id, giniC_cpm)),
#                               read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_1.01_preferential_filtered.tsv') %>% rename(transcript_id = TXNAME))
#write_tsv(pref_lfc_ginis %>% filter(!is.na(n_tis_pref)) %>% arrange(n_tis_pref, giniC, giniC_cpm), 'maddy_checking_gini_on_pref.tsv')
#
########################################################################
########################################################################
################## misc testing stuff ########################
#print(normalized_counts_ginic %>% select(transcript_id, giniC, giniC_nona))
#
## probably want to filter by a CPM? don't need to look at isoforms where no tissues are above 1
##TODO: I'll want to decide a threshold for the gini genes. lets look at > .3
#
#plot(Lc(normalized_counts %>% filter(transcript_id == 'BambuTx995_BambuGene154124') %>% pivot_longer(!transcript_id, names_to = 'sample', values_to = 'cpm') %>% pull(cpm)))
#plot(Lc(normalized_counts %>% filter(transcript_id == 'BambuTx982_ENSG00000168899') %>% pivot_longer(!transcript_id, names_to = 'sample', values_to = 'cpm') %>% pull(cpm)))
#plot(Lc(normalized_counts %>% filter(transcript_id == 'BambuTx983_ENSG00000115523') %>% pivot_longer(!transcript_id, names_to = 'sample', values_to = 'cpm') %>% pull(cpm)))
#plot(Lc(normalized_counts %>% filter(transcript_id == 'BambuTx1180_ENSG00000100319') %>% pivot_longer(!transcript_id, names_to = 'sample', values_to = 'cpm') %>% pull(cpm)))
#
