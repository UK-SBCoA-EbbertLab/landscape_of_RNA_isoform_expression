#!/usr/bin/Rscript

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

cpat_file <- args[1]

# Fix CPAT3 format

cpat <- read_tsv(cpat_file) %>%
	select(seq_ID, mRNA, ORF, Fickett, Hexamer, Coding_prob) %>%
	rename(id=seq_ID, mRNA_size = mRNA, ORF_size=ORF, Fickett_score = Fickett, Hexamer_score=Hexamer, coding_prob=Coding_prob)

write_tsv(cpat, "z_output/cpat3_output.formatted.tsv")

