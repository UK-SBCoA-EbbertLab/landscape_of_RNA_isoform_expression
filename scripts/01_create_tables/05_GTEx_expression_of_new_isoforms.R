library(tidyverse)

# load in the list of medically relevant genes and grab the ids to use for filtering later
med_rel_genes <- read_tsv("../../references/medically_relevant_genes.tsv")
med_rel_ids <- med_rel_genes %>% pull(gene_id)

#get the number of isoforms from medically relevant genes that we discovered
gtf_with_bambu <- read_tsv('../../nextflow_pipeline/references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109.gtf', comment='#', col_names = c('chr', 'source', 'type', 'start', 'end', 'dot_1', 'strand', 'dot_2', 'other') ) %>%
	filter(source == 'Bambu') %>%
	filter(type == 'transcript') %>%
	extract(other, c('gene_id'), "gene_id \"([[:alnum:]]+)\"; .*", remove= FALSE) %>%
	extract(other, c('transcript_id'), ".* transcript_id \"([[:alnum:]]+)\";.*") %>%
	filter(gene_id %in% med_rel_ids)

all_med_rel <- nrow(gtf_with_bambu)
print(all_med_rel)

# load in all the isoforms across 4 thresholds and filter for only Bambu
# tell us the status of the isoforms and if they come from medically relevant genes
# count how many tissues the isoform is expressed in across the thresholds
isoforms <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_0_1_5_10.tsv') %>%
	filter(grepl('Bambu', transcript_id))  %>%
	mutate(status = ifelse(grepl('ENSG', gene_id), 'nfk', 'nfn')) %>%
	mutate(is_med_rel = ifelse(gene_id %in% med_rel_ids, TRUE, FALSE)) %>%
	pivot_longer(!c(transcript_id, gene_id, gene_name, gene_biotype, threshold, status, is_med_rel), names_to = 'tissue', values_to = 'median_CPM') %>%
	drop_na(median_CPM) %>%
	group_by(transcript_id, gene_id, threshold) %>%
	mutate(n_tissues = n())

write_tsv(isoforms, "../../tables/GTEx_expression_our_new_isoforms/GTEx_expression_of_new_isoforms.tsv")

# Look at isoforms found in all 9 tissues and just 1 tissue
nine_tissues <- isoforms %>% 
	select(!c(tissue, median_CPM)) %>%
	filter(n_tissues == 9)

write_tsv(nine_tissues, "../../tables/GTEx_expression_our_new_isoforms/tx_in_nine_gtex_tissues.tsv")

one_tissue <- isoforms %>%
	filter(n_tissues == 1)

write_tsv(one_tissue, "../../tables/GTEx_expression_our_new_isoforms/tx_in_single_gtex_tissue.tsv")

