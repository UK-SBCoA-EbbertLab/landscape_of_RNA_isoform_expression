library(tidyverse)

tissues_to_use = c("Brain - Cerebellar Hemisphere",
		   "Brain - Frontal Cortex (BA9)",
		   "Brain - Putamen (basal ganglia)",
		   "Cells - Cultured fibroblasts",
		   "Heart - Atrial Appendage",
		   "Heart - Left Ventricle",
		   "Liver",
		   "Lung",
		   "Muscle - Skeletal")

# samples we are not using due to replicates, low read depth, and poor PCA clustering\n",
filter_out_samples <- c(
    "GTEX-Q2AG-0011-R11A-SM-2EBL2_rep2.FAK44637.bam",
    "GTEX-Q2AG-0011-R11A-SM-2EBL2_rep.FAK49243.bam",
    "GTEX-T5JC-0011-R10A-SM-2TT23.FAK91589.bam",
    "GTEX-QEG5-0008-SM-3QHW2_exp.FAK30166.bam",
    "GTEX-QV44-0008-SM-3QNG7_ctrl1.FAK55556.bam",
    "GTEX-QV44-0008-SM-3QNG7_exp.FAK52124.bam",
    "GTEX-RWS6-0008-SM-3QHWG_rep.FAK49207.bam",
    "GTEX-S4Z8-0008-SM-2Y983_exp1.FAK55723.bam",
    "GTEX-S4Z8-0008-SM-2Y983_exp2.FAK47416.bam",
    "GTEX-S95S-0008-SM-3RQ8B_exp1.FAK55217.bam",
    "GTEX-S95S-0008-SM-3RQ8B_exp2.FAK47088.bam",
    "GTEX-WY7C-0008-SM-3NZB5_ctrl.FAK55679.bam",
    "GTEX-1GN1W-0226-SM-7AGLJ_rep.FAK91654.bam",
    "GTEX-WY7C-1126-SM-3GS2X_rep2.FAK49168.bam",
    "GTEX-WY7C-1126-SM-3GS2X.FAK39149.bam",
    "GTEX-Y5LM-0426-SM-3YX99.FAK52212.bam",
    "GTEX-Y5LM-0426-SM-3YX99_rep2.FAK41279.bam",
    "GTEX-14BMU-0526-SM-5CA2F.FAK44778.bam",
    "GTEX-14BMU-0526-SM-5CA2F_rep.FAK93376.bam",
    "GTEX-13QJ3-0726-SM-7LDHS.FAK49189.bam",
    "GTEX-ZT9X-1826-SM-4V2KV_rep.FAK39773.bam",
    "GTEX-ZT9X-1826-SM-4V2KV.FAK49260.bam",
    "GTEX-WY7C-0726-SM-3GLGQ.FAK46872.bam")



files_flagstat <- list.files('../data/raw/GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/bam_filtering/', pattern = ".*.flagstat", full.names=TRUE)

mapped_reads = tibble(sample = character(), n_mapped_reads_str = character())

for (f in files_flagstat) {
	mapped_reads <- mapped_reads %>%
		add_row(sample = f, 
			n_mapped_reads_str = str_extract(read_file(f), "[0-9]* . [0-9] primary mapped"))
}

mapped_reads <- mapped_reads %>%
		extract(n_mapped_reads_str, 'n_mapped_reads', "([0-9]*) . [0-9] primary mapped", convert=TRUE) %>%
		mutate(sample = str_remove(tools::file_path_sans_ext(basename(sample)), '_filtered_mapq_10'))

mapped_reads <- read_tsv('../data/GTEx_v9_ONT_metadata.txt') %>%
	select(sample_id, tissue_site_detail, bam_file) %>%
#	filter(!(bam_file %in% filter_out_samples)) %>%
	mutate(sample = str_remove(bam_file, '.bam')) %>%
	right_join(mapped_reads) %>%
	select(sample, n_mapped_reads, tissue_site_detail) %>%
#	filter(tissue_site_detail %in% tissues_to_use) %>%
	arrange(tissue_site_detail, sample)
#	arrange(n_mapped_reads)

write_tsv(mapped_reads, '../tables/GTEx_n_mapped_reads_per_sample_all_samples.tsv')
print(median(mapped_reads$n_mapped_reads))
