library(tidyr)
library(readr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
gtf_output <- args[2]
annotation_key <- args[3]

# load in the gtf
# grab the isoform and gene ids
gtf_table <- read_tsv(gtf_file, col_names = FALSE, comment="#") %>%
	mutate(isoform_id = str_match(X9, 'transcript_id "([^"]+)"')[,2]) %>%
	mutate(gene_id = str_match(X9, 'gene_id "([^"]+)"')[,2]) %>%
	mutate(isoform_id_s = str_match(X9, '(transcript_id "[^"]+")')[,2]) %>%
	mutate(gene_id_s = str_match(X9, '(gene_id "[^"]+")')[,2])

# grab all the info from X9 that we need for gene entries
gene_info <- gtf_table %>%
	filter(str_starts(gene_id, "ENSG")) %>%
	filter(X3 == 'gene') %>%
	select(X9, gene_id, gene_id_s) %>%
	mutate(gene_name = str_match(X9, 'gene_name "([^"]+)"')[,2]) %>%
	mutate(gene_name_s = str_match(X9, '(gene_name "[^"]+")')[,2]) %>%
	mutate(gene_version = str_match(X9, '(gene_version "[^"]+")')[,2]) %>%
	mutate(gene_source = str_match(X9, '(gene_source "[^"]+")')[,2]) %>%
        mutate(gene_biotype = str_match(X9, 'gene_biotype "([^"]+)"')[,2]) %>%
        mutate(gene_biotype_s = str_match(X9, '(gene_biotype "[^"]+")')[,2]) %>%
	select(gene_id, gene_id_s, gene_version, gene_name, gene_name_s, gene_source, gene_biotype, gene_biotype_s) %>%
	distinct()

# grab all the new bambu entries and give them gene information
bambu_gtf <- gtf_table %>%
	filter(str_starts(isoform_id, "Bambu")) %>%
	left_join(gene_info) %>%
	mutate(gene_name = coalesce(gene_name, gene_id)) %>%
	mutate(gene_name_s = coalesce(gene_name_s, paste0("gene_name ", gene_id))) %>%
	replace_na(list(gene_version = 'gene_version "1"', 
			gene_source = 'gene_source "bambu"',
			gene_biotype = "unknown",
			gene_biotype_s = 'gene_biotype "unknown"'))

# give bambu entries isoform information
bambu_tx_gtf <- bambu_gtf %>%
	mutate(transcript_version = 'transcript_version "1"') %>%
	unite("transcript_name", c(gene_name, isoform_id), sep = '-', remove = FALSE) %>%
	mutate(transcript_name = paste0('transcript_name "', transcript_name, '"')) %>%
	mutate(transcript_source = 'transcript_source "bambu"') %>%
	mutate(transcript_biotype = 'transcript_biotype "novel_unknown"')

# grab all the bambu exon entries and give them exon information
bambu_exon_gtf <- bambu_tx_gtf %>%
	filter(X3 == 'exon') %>%
	mutate(exon_number_s = str_match(X9, '(exon_number "[^"]+")')[,2]) %>%
	mutate(exon_number = str_match(X9, 'exon_number "([^"]+)"')[,2]) %>%
	unite("exon_id", c(isoform_id, exon_number), sep = "_E", remove = FALSE) %>%
	mutate(exon_id = paste0('exon_id "', exon_id, '"')) %>%
	mutate(exon_version = '"exon_version "1"') %>%
	unite("X9", c(gene_id_s, gene_version, isoform_id_s, transcript_version, exon_number_s, gene_name_s, gene_source, gene_biotype_s, transcript_name, transcript_source, transcript_biotype, exon_id, exon_version), sep = "; ", na.rm = TRUE) %>%
	mutate(X9 = paste0(X9, ';')) %>%
	select(-c(isoform_id, gene_id, gene_biotype))

# remove the exons from the bambu isoform entries and recreate the X9 column
bambu_tx_gtf_2 <- bambu_tx_gtf %>%
	filter(X3 != 'exon') %>%
	unite("X9", c(gene_id_s, gene_version, isoform_id_s, transcript_version, gene_name_s, gene_source, gene_biotype_s, transcript_name, transcript_source, transcript_biotype), sep = "; ", na.rm = TRUE) %>%
	mutate(X9 = paste0(X9, ';')) %>%
	select(-c(isoform_id, gene_id, gene_biotype))

# filter out bambu ids from the original gtf, add back in the bambu isoforms and exons, sort, and remove the double quotes 
combined <- gtf_table %>%
	filter(!grepl('Bambu', X9)) %>%
	select(c(X1, X2, X3, X4, X5, X6, X7, X8, X9)) %>%
	bind_rows(bambu_tx_gtf_2) %>%
	bind_rows(bambu_exon_gtf) %>%
	select(-c(exon_number, gene_name)) %>%
	arrange(X1, X4) %>%
	mutate(X9 = str_replace_all(X9, '"', ''))

write_tsv(combined, gtf_output, col_names = FALSE)

bambu_isoform_key_info <- bambu_tx_gtf %>%
	filter(X3 == 'transcript') %>%
	mutate(transcript_biotype = "novel_unknown") %>%
	select(isoform_id, gene_id, gene_name, gene_biotype, transcript_biotype)

ensembl_isoform_key_info <- gtf_table %>%
	filter(str_starts(isoform_id, "ENST")) %>%
	filter(X3 == 'transcript') %>%
        left_join(gene_info) %>%
        mutate(gene_name = coalesce(gene_name, gene_id)) %>%
	mutate(transcript_biotype = str_match(X9, 'transcript_biotype "([^"]+)"')[,2]) %>%
	select(isoform_id, gene_id, gene_name, gene_biotype, transcript_biotype)

combined_iso_key_info <- ensembl_isoform_key_info %>%
	bind_rows(bambu_isoform_key_info)

write_tsv(combined_iso_key_info, annotation_key)
	


