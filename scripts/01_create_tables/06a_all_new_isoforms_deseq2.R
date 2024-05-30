library("DESeq2")
library(tidyverse)

########## format data for DESeq2 ##################################
# get sample_id into the correct format
get_samp_id <- function(str) {
        name_array <- unlist(strsplit(str, '.', fixed = TRUE))
        last_element <- unlist(strsplit(name_array[length(name_array)], '_', fixed = TRUE))[1]
        string_1 <- paste0(head(name_array, -1), collapse='-')
        string_2 <- paste(last_element, 'bam', sep = ".")
        return(paste(string_1, string_2, sep = "."))
}

# list of the tissues we will be using
tissues_to_use = c('Muscle - Skeletal', 'Lung', 'Liver',
                   'Heart - Left Ventricle', 'Heart - Atrial Appendage', 'Cells - Cultured fibroblasts',
                   'Brain - Putamen (basal ganglia)', 'Brain - Frontal Cortex (BA9)', 'Brain - Cerebellar Hemisphere')

filter_out_samples = c("GTEX-Q2AG-0011-R11A-SM-2EBL2_rep2.FAK44637.bam",
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

# load the gtex_metadata and filter out any 'direct' samples
gtex_metadata <- read_tsv('../../data/GTEx_v9_ONT_metadata.txt') %>%
        filter(!(grepl('direct', sample_id))) %>%
        mutate(sample_id = as.factor(sample_id))

#load the gtex all_counts matrix, fix the sample names, calculate CPM, join with the metadata, filter out tissues/samples that we are not using, and sort the tibble
gtex_counts <- read_tsv('../../data/raw/GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/bambu_quant/counts_transcript.txt') %>%
#	filter(grepl('Bambu', TXNAME)) %>%
        pivot_longer(cols = !c(TXNAME, GENEID), names_to = "SAMPID", values_to = "counts") %>%
        rowwise() %>%
        mutate(SAMPID = get_samp_id(SAMPID)) %>%
        ungroup() %>%
#        group_by(SAMPID) %>%
#        mutate(total_samp_counts = sum(counts)) %>%
#        ungroup() %>%
#        rowwise() %>%
#        mutate(CPM = round((counts/total_samp_counts)*1000000, 2)) %>%
        inner_join(gtex_metadata %>% select(bam_file, tissue_site_detail) %>% mutate(SAMPID = bam_file)) %>%
        filter(tissue_site_detail %in% tissues_to_use) %>%
        filter(!(bam_file %in% filter_out_samples))  %>%
        arrange(TXNAME, tissue_site_detail) %>%
	select(-c(bam_file, tissue_site_detail)) %>%
        pivot_wider(names_from = "SAMPID", values_from = "counts") %>%
	mutate(across(where(is.numeric), round, 0))

sample_order = colnames(gtex_counts)[-c(1:2)]
gtex_metadata <- gtex_metadata %>%
	filter(tissue_site_detail %in% tissues_to_use) %>%
	filter(!(bam_file %in% filter_out_samples)) %>%
	mutate(sample_id = fct_relevel(bam_file, sample_order)) %>%
	arrange(sample_id)

# check to see if the samples are in the same order
if(!all(gtex_metadata$sample_id == sample_order)) {
	quit(status=1)
}

gtex_metadata = as.data.frame(gtex_metadata)

########### DESEQ2 ##############################################################

	combined <- gtex_counts %>% 
			  rename(transcript_id = TXNAME) %>% 
			  rename(gene_id = GENEID) %>%
		unite('id', c(transcript_id, gene_id))
	print(combined)

	# start deseq2
	#hk_counts <- as.data.frame(housekeeping_deseq2_input)
	comb_counts <- as.data.frame(combined)
	# pass the dataframe to deseq2
	#dds <- DESeqDataSetFromMatrix(countData = hk_counts,
	dds <- DESeqDataSetFromMatrix(countData = comb_counts,
				      colData = gtex_metadata,
				      design=~tissue_site_detail, tidy=TRUE)
	print(dds)
	dds <- DESeq(dds)
	print(dds)
	genes_to_keep <- which(grepl('Bambu', rownames(dds)))
	dds <- dds[genes_to_keep, ]

	#added after original submission #############
	write_tsv(as.data.frame(counts(dds, normalized=TRUE)) %>% rownames_to_column('transcript_id'), 'deseq2_normalized_values_bambu_iso.tsv')

	##############################################

	# to write out pairwise 
	all_the_results <- tibble()
	# run pairwide comparisons
	for (i in 1:8){
		print(i)
		n=i+1
		print(n)
		for (j in n:9) {
			print(j)
			t1 = tissues_to_use[i]
			t2 = tissues_to_use[j]
			res <- results(dds, contrast=c("tissue_site_detail", t1,t2))
			res <- as_tibble(as.data.frame(res), rownames = 'TXNAME')
			all_the_results <- bind_rows(all_the_results,
						     res %>%
							     mutate(tis1=t1) %>%
							     mutate(tis2=t2))
		}
	}
	write_tsv(all_the_results %>%
		  arrange(TXNAME), paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_after_submission.tsv'))
		  #arrange(TXNAME), paste0('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms.tsv'))


