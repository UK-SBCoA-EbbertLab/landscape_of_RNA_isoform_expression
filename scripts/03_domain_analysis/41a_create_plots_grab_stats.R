#!/usr/bin/Rscript

library(IsoformSwitchAnalyzeR)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])

# Create empty data frames to collect results
all_features <- data.frame()
all_switches <- data.frame()
all_consequences <- data.frame()

for (gene_name in c("ALB", "ARPP21", "CACNA1A", "CRELD1", "CSF3R", "DGUOK", "EEF1AKMT1", "MAOB", "MIR9-1HG", "OAZ1", "PAX6", "SNCA", "TNNT2", "TPM1")) {
	isoforms <- unique(mySwitchListAnalyzed$isoformFeatures[mySwitchListAnalyzed$isoformFeatures$gene_id == gene_name, ]$isoform_id)
	n_isoforms <- length(isoforms)

	gene_plot <- switchPlotTranscript(mySwitchListAnalyzed, 
					  gene = gene_name,
					  IFcutoff = 0,
					  reverseMinus = FALSE)
	ggsave(paste0("z_output/plots/", gene_name, "_domain_transcripts_plot.pdf"), gene_plot, width = 10, height = ifelse(n_isoforms <= 4, n_isoforms, floor((n_isoforms - 2) / 2) + 4))
	
	# Subset and append data
	gene_features <- filter(mySwitchListAnalyzed$isoformFeatures, gene_id == gene_name)
	iso_and_gene <- gene_features %>% select(gene_name, isoform_id) %>% distinct()
	gene_switches <- filter(mySwitchListAnalyzed$isoformSwitchAnalysis, isoform_id %in% isoforms) %>% left_join(iso_and_gene)
	gene_consequences <- filter(mySwitchListAnalyzed$switchConsequence, gene_id == gene_name)

	all_features <- bind_rows(all_features, gene_features)
	all_switches <- bind_rows(all_switches, gene_switches)
	all_consequences <- bind_rows(all_consequences, gene_consequences)
}

# Save the final combined tables
write_tsv(all_features, "z_output/isoformFeatures_FILTERED.tsv")
write_tsv(all_switches, "z_output/isoformSwitchAnalysis_FILTERED.tsv")
write_tsv(all_consequences, "z_output/isoformSwitchConsequences_FILTERED.tsv")
