#!/usr/bin/Rscript

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)

deeploc_file <- args[1]

# Fix deeploc2

# Read all CSV files, concatenate them, and select desired columns
deeploc2_combined <- read_csv(deeploc_file) %>%
	select(-c(`Membrane types`, Peripheral, Transmembrane, `Lipid anchor`, Soluble)) %>%
	write_tsv("z_output/deeploc2_output_formatted.tsv")
