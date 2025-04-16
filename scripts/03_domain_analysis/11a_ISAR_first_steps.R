#!/usr/bin/Rscript

library(IsoformSwitchAnalyzeR)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

counts <- args[1]
rel_expr <- args[2]
design_matrix <- args[3]
annotation <- args[4]
transcriptome <- args[5]

# Loads the design matrix
designMatrix <- as.data.frame(read_tsv(design_matrix))

# Loads the isoform count matrix, 
countsMatrix <- as.data.frame(read_tsv(counts))

# Loads the isoform_rep_expression matrix,
cpmMatrix <- as.data.frame(read_tsv(rel_expr))

# import all the data and prep the switch object
mySwitchList <- importRdata(countsMatrix,
			  cpmMatrix, 
			  designMatrix,
                          annotation,
                          transcriptome,
                          ignoreAfterPeriod = TRUE)

# filter out irrelevant isoforms to improve performance
mySwitchList <- preFilter(
				mySwitchList,
				geneExpressionCutoff = NULL,
				#isoformExpressionCutoff = 0.5,
				isoformExpressionCutoff = 1, # original value
				IFcutoff = NULL,
				removeSingleIsoformGenes = TRUE,
				keepIsoformInAllConditions = TRUE)

# Identify isoform switches
mySwitchList <- isoformSwitchTestSatuRn(
				switchAnalyzeRlist = mySwitchList, 
				reduceToSwitchingGenes = FALSE,
				reduceFurtherToGenesWithConsequencePotential = FALSE)


# Analyze open reading frame for new isoforms (ISAR does this even though we get this info from CPAT3 later....)
mySwitchList <- analyzeNovelIsoformORF(
				mySwitchList,
				analysisAllIsoformsWithoutORF = TRUE,
				minORFlength = 50,
				orfMethod = 'longest.AnnotatedWhenPossible')

mySwitchList <- extractSequence(
    mySwitchList,
    onlySwitchingGenes = FALSE,
    pathToOutput = 'z_references/.'
)


extractSwitchSummary(mySwitchList)

save(mySwitchList, file="z_output/ISAR_part1_object.RData")

