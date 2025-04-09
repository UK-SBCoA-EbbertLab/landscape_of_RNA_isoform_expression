#!/usr/bin/Rscript

library(readr)
library(tidyr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

counts <- args[1]
rel_expr <- args[2]
design_matrix <- args[3]
annotation <- args[4]
gene_list <- args[5]

print(counts)
print(gene_list)
print(rel_expr)
print(annotation)
print(design_matrix)

# Load list of genes to keep
genes_to_keep <- read_tsv(gene_list, col_names = c("gene_name"))

print("loaded genes_to_keep")

# load annotation 
annotation <- read_tsv(annotation, col_names = FALSE, comment="#") %>%
	mutate(gene_id = str_match(X9, 'gene_id "([^"]+)"')[,2]) %>%
	mutate(gene_name = str_match(X9, 'gene_name "([^"]+)"')[,2])

genes_to_keep <- genes_to_keep %>%
	left_join(annotation %>% select(gene_name, gene_id) %>% distinct()) %>%
	distinct() %>%
	mutate(gene_id = if_else(is.na(gene_id), gene_name, gene_id))

annotation <- annotation %>%
	filter(gene_id %in% genes_to_keep$gene_id) %>%
	select(-c(gene_id, gene_name))

print("annotation loaded and formatted")

write_tsv(annotation, "z_references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109_FILTERED.gtf", col_names = FALSE, escape = "none")

#write_tsv(genes_to_keep, "../data/gene_and_iso_IDs_for_domain_analysis.tsv", col_names = FALSE)

# Loads the design matrix
designMatrix <- 
  read_tsv(design_matrix)


# Loads the isoform count matrix, 
#  normalizes col names with the design matrix, 
#  filters columns, and
#  Renames "TXNAME" to "isoform_id"
countsMatrix <- 
  read_tsv(counts) %>%
  filter(GENEID %in% genes_to_keep$gene_id) %>%
  rename_all(function(x) gsub("_filtered_mapq_10$", '', x)) %>%
  select(any_of(c("TXNAME", designMatrix$sampleID))) %>%
  rename(isoform_id = TXNAME)

write_tsv(countsMatrix, "z_data/counts_transcript_FILTERED.txt")

# Loads the isoform_rep_expression matrix,
#  normalizes col names with the design matrix, 
#  filters columns, and
#  Renames "TXNAME" to "isoform_id"
cpmMatrix <- 
  read_tsv(rel_expr) %>%
  filter(GENEID %in% genes_to_keep$gene_id) %>%
  rename_all(function(x) gsub("_filtered_mapq_10$", '', x)) %>%
  select(any_of(c("TXNAME", designMatrix$sampleID))) %>%
  rename(isoform_id = TXNAME)

write_tsv(cpmMatrix, "z_data/CPM_transcript_FILTERED.txt")


