library("DESeq2")
library(tidyverse)

# get the current date and create a directory for it if it doesn't already exist
the_date <- as.character(Sys.Date())
print(the_date)
if (!dir.exists(paste0('deseq_out/', the_date))) dir.create(paste0('deseq_out/', the_date))

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


	
	
# load the gtex_metadata and filter out any 'direct' samples
gtex_metadata <- read_tsv('/pscratch/mteb223_uksr/new_RNA_isoform_expression_across_tissues/data/GTEx_v9_ONT_metadata.txt') %>%
	filter(!(grepl('direct', sample_id))) %>%
	mutate(sample_id = as.factor(sample_id))
head(gtex_metadata)


# load the gtex_counts and filter out the tissues we aren't using
gtex_counts <- read_tsv('../data/raw/GRCh38_quant_mapq10_gtf_109_and_high-confidence_GTEx_DATA/bambu_quant/counts_transcript.txt') %>%
	pivot_longer(cols = !c(TXNAME, GENEID), names_to = "SAMPID", values_to = "counts") %>%
	rowwise() %>%
	mutate(SAMPID = get_samp_id(SAMPID)) %>%
	ungroup() %>%
	inner_join(gtex_metadata %>% select(bam_file, tissue_site_detail) %>% mutate(SAMPID = bam_file)) %>%
	filter(tissue_site_detail %in% tissues_to_use) %>%
	arrange(TXNAME, tissue_site_detail) 

# prep tibble to be filtered later
gtex_counts_filtered <- gtex_counts
gtex_metadata_filtered <- gtex_metadata

# pivot the data wider to be a matrix
gtex_counts <- gtex_counts %>% 
	select(TXNAME, SAMPID, counts) %>%
	pivot_wider(names_from = 'SAMPID', values_from = 'counts') %>%
	mutate(across(where(is.numeric), round, 0))

head(gtex_counts)

# get the order of the samples from the tibble
sample_order = colnames(gtex_counts)[-c(1)]

# check to see if the metadata samples and the samples from the tibble are in the same order
all(gtex_metadata$bam_file == sample_order)

#filter out the samples we won't be using, then arrange the samples we will be using in the same order as those from the tibble
gtex_metadata <- gtex_metadata %>%
	filter(tissue_site_detail %in% tissues_to_use) %>%
	mutate(sample_id = fct_relevel(bam_file, sample_order)) %>%
	arrange(sample_id)

# check again that the samples are in the same order
if(!all(gtex_metadata$sample_id == sample_order)) {
	quit(status=1)
}

############# Start DESeq2 #########################################

# turn the tibbles into data frames
gtex_counts = as.data.frame(gtex_counts)
gtex_metadata = as.data.frame(gtex_metadata)

# pass the dataframes to deseq2 
dds <- DESeqDataSetFromMatrix(countData = gtex_counts,
			      colData = gtex_metadata,
			      design=~tissue_site_detail, tidy = TRUE)

print(dds)

# run deseq2
dds <- DESeq(dds)

# create the stuff for the PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("tissue_site_detail"), returnData = TRUE)
print(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))
print(percentVar)
a_plot <- ggplot(pcaData, aes(PC1, PC2, color = tissue_site_detail)) +
	geom_point() +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	coord_fixed()

ggsave(paste0('deseq_out/', the_date, "/samps_all_isoforms_PCA_gtex.pdf"), plot=a_plot)
ggsave("deseq_out/samps_all_isoforms_PCA_gtex.pdf", plot=a_plot)

# list of samples that we are NOT using, due to low depth, replicates, or bad PCA clustering
filter_out_samples = c(
    "GTEX-Q2AG-0011-R11A-SM-2EBL2_rep2.FAK44637.bam", # brain - cerebellar hemisphere
    "GTEX-Q2AG-0011-R11A-SM-2EBL2_rep.FAK49243.bam", # brain - cerebellar hemisphere
    "GTEX-T5JC-0011-R10A-SM-2TT23.FAK91589.bam", # brain - frontal cortex
    "GTEX-QEG5-0008-SM-3QHW2_exp.FAK30166.bam", # cells - cultured fibroblasts
    "GTEX-QV44-0008-SM-3QNG7_ctrl1.FAK55556.bam", # cells - cultures fibroblasts
    "GTEX-QV44-0008-SM-3QNG7_exp.FAK52124.bam", # cells - cultures fibroblasts
    "GTEX-RWS6-0008-SM-3QHWG_rep.FAK49207.bam", # cells - cultures fibroblasts
    "GTEX-S4Z8-0008-SM-2Y983_exp1.FAK55723.bam", # cells - cultures fibroblasts
    "GTEX-S4Z8-0008-SM-2Y983_exp2.FAK47416.bam", # cells - cultures fibroblasts
    "GTEX-S95S-0008-SM-3RQ8B_exp1.FAK55217.bam", # cells - cultures fibroblasts
    "GTEX-S95S-0008-SM-3RQ8B_exp2.FAK47088.bam", # cells - cultures fibroblasts
    "GTEX-WY7C-0008-SM-3NZB5_ctrl.FAK55679.bam", # cells - cultures fibroblasts
    "GTEX-1GN1W-0226-SM-7AGLJ_rep.FAK91654.bam", # heart - atrial appendage
    "GTEX-WY7C-1126-SM-3GS2X_rep2.FAK49168.bam", # heart - atrial appendage
    "GTEX-WY7C-1126-SM-3GS2X.FAK39149.bam", # heart - atrial appendage
#    "GTEX-WY7C-1126-SM-3GS2X_rep.FAK49218.bam", # heart - atrial appendage 
    "GTEX-Y5LM-0426-SM-3YX99.FAK52212.bam", # liver
    "GTEX-Y5LM-0426-SM-3YX99_rep2.FAK41279.bam", # liver
    "GTEX-14BMU-0526-SM-5CA2F.FAK44778.bam", #lung
    "GTEX-14BMU-0526-SM-5CA2F_rep.FAK93376.bam", #lung
    "GTEX-13QJ3-0726-SM-7LDHS.FAK49189.bam", # muscle - skeletal
    "GTEX-ZT9X-1826-SM-4V2KV_rep.FAK39773.bam", # muscle - skeletal
    "GTEX-ZT9X-1826-SM-4V2KV.FAK49260.bam" # muscle - skeletal

    )

filtered_samples_plus_clustered_bad = append(filter_out_samples, "GTEX-WY7C-0726-SM-3GLGQ.FAK46872.bam") # liver

############# Start DESeq2 #########################################

# format the data from deseq
gtex_counts_filtered_2 <- gtex_counts_filtered %>%
	filter(!(bam_file %in% filter_out_samples))  %>%
	select(TXNAME, SAMPID, counts) %>%
	pivot_wider(names_from = 'SAMPID', values_from = 'counts') %>%
	mutate(across(where(is.numeric), round, 0))


# get the order of the samples from the tibble
sample_order = colnames(gtex_counts_filtered_2)[-c(1)]

# check to see if the metadata samples and the samples from the tibble are in the same order
all(gtex_metadata_filtered$bam_file == sample_order)

#filter out the samples we won't be using, then arrange the samples we will be using in the same order as those from the tibble
gtex_metadata_filtered_2 <- gtex_metadata_filtered %>%
	filter(tissue_site_detail %in% tissues_to_use) %>%
	filter(!(bam_file %in% filter_out_samples)) %>%
	mutate(sample_id = fct_relevel(bam_file, sample_order)) %>%
	arrange(sample_id)

# check again that the samples are in the same order
if(!all(gtex_metadata_filtered_2$sample_id == sample_order)) {
	quit(status=1)
}


# turn the tibbles into data frames
gtex_counts_filtered_2 = as.data.frame(gtex_counts_filtered_2)
gtex_metadata_filtered_2 = as.data.frame(gtex_metadata_filtered_2)



# pass the dataframes to deseq2 
dds <- DESeqDataSetFromMatrix(countData = gtex_counts_filtered_2,
			      colData = gtex_metadata_filtered_2,
			      design=~tissue_site_detail, tidy = TRUE)

print(dds)

# run deseq2
dds <- DESeq(dds)

# create the stuff for the PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("tissue_site_detail"), returnData = TRUE)
print(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))
print(percentVar)
f_plot <- ggplot(pcaData, aes(PC1, PC2, color = tissue_site_detail)) +
	geom_point() +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	coord_fixed()

ggsave(paste0('deseq_out/', the_date, "/filtered_out_exp_and_rep_samps_all_isoforms_PCA_gtex.pdf"), plot=f_plot)
ggsave("deseq_out/filtered_out_exp_and_rep_samps_all_isoforms_PCA_gtex.pdf", plot=f_plot)


############# Start DESeq2 - filter out all 'bad' samples #########################################

gtex_counts_filtered <- gtex_counts_filtered %>%
	filter(!(bam_file %in% filtered_samples_plus_clustered_bad))  %>%
	select(TXNAME, SAMPID, counts) %>%
	pivot_wider(names_from = 'SAMPID', values_from = 'counts') %>%
	mutate(across(where(is.numeric), round, 0))


# get the order of the samples from the tibble
sample_order = colnames(gtex_counts_filtered)[-c(1)]
#sample_order <- sapply(sample_order[-c(1:2)], get_samp_id)

# check to see if the metadata samples and the samples from the tibble are in the same order
all(gtex_metadata$bam_file == sample_order)

#filter out the samples we won't be using, then arrange the samples we will be using in the same order as those from the tibble
gtex_metadata_filtered <- gtex_metadata_filtered %>%
	filter(tissue_site_detail %in% tissues_to_use) %>%
	filter(!(bam_file %in% filtered_samples_plus_clustered_bad)) %>%
	mutate(sample_id = fct_relevel(bam_file, sample_order)) %>%
	arrange(sample_id)

# check again that the samples are in the same order
if(!all(gtex_metadata_filtered$sample_id == sample_order)) {
	quit(status=1)
}


# turn the tibbles into data frames
gtex_counts_filtered = as.data.frame(gtex_counts_filtered)
gtex_metadata_filtered = as.data.frame(gtex_metadata_filtered)



# pass the dataframes to deseq2 
dds <- DESeqDataSetFromMatrix(countData = gtex_counts_filtered,
			      colData = gtex_metadata_filtered,
			      design=~tissue_site_detail, tidy = TRUE)

print(dds)

# run deseq2
dds <- DESeq(dds)

# create the stuff for the PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("tissue_site_detail"), returnData = TRUE)
print(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))
print(percentVar)
f_plot <- ggplot(pcaData, aes(PC1, PC2, color = tissue_site_detail)) +
	geom_point() +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	coord_fixed()

ggsave(paste0('deseq_out/', the_date, "/filtered_samps_all_isoforms_PCA_gtex.pdf"), plot=f_plot)
ggsave("deseq_out/filtered_samps_all_isoforms_PCA_gtex.pdf", plot=f_plot)

