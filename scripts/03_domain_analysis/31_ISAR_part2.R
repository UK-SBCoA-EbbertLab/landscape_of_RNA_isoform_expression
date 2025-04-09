#!/usr/bin/Rscript

library(IsoformSwitchAnalyzeR)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

rdata_obj <- args[1]
cpat3 <- args[2]
pfam <- args[3]
iupred2a <- args[4]
signalp6 <- args[5]
deeploc2 <- args[6]
deeptmhmm <- args[7]

load(rdata_obj)

mySwitchList <- analyzeCPAT(mySwitchList, 
			    cpat3,
			    codingCutoff = 0.364, # using suggested cutoff from CPAT, ISAR suggests 0.725 for humans
			    removeNoncodinORFs = TRUE)
			    #removeNoncodinORFs = FALSE)

mySwitchList <- analyzePFAM(mySwitchList,
			    pfam)


mySwitchList <- tryCatch({
	analyzeSignalP(mySwitchList, signalp6)
}, error = function(e) {
	if (grepl("No signial peptides were found", e$message)) {
		message("No signal peptides found, continuing with original data.")
		return(mySwitchList)
	} else {
		stop(e)
	}
})

mySwitchList <- analyzeIUPred2A(mySwitchList,
				iupred2a)

mySwitchList <- analyzeDeepLoc2(mySwitchList,
				deeploc2)

if (deeptmhmm != "None") {
	mySwitchList <- analyzeDeepTMHMM(mySwitchList,
					  deeptmhmm)
}

mySwitchListAnalyzed <- analyzeAlternativeSplicing(mySwitchList)
mySwitchListAnalyzed <- analyzeSwitchConsequences(mySwitchListAnalyzed)

switchingISO <- extractTopSwitches(mySwitchListAnalyzed,
				   filterForConsequences = FALSE,
				   n = NA,
				   extractGenes = FALSE,
				   sortByQvals = TRUE)

switchPlotTopSwitches(mySwitchListAnalyzed,
		      n = NA,
		      filterForConsequences = FALSE,
		      fileType = "pdf",
		      pathToOutput = "z_output/plots/")

save(mySwitchList, mySwitchListAnalyzed, switchingISO, file="z_output/ISAR_part2.RData")

