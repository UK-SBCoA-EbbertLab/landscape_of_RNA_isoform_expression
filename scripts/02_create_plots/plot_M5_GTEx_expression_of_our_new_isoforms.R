library(tidyverse)

# number of isoforms in each category
all_nfk = 428
all_nfn = 267
all_nfmito = 5
all_nfmed_rel = 54

# load in data
new_isoforms <- read_tsv("../../tables/GTEx_expression_our_new_isoforms/GTEx_expression_of_new_isoforms.tsv") %>%
	mutate(n_tissues = as.character(n_tissues)) %>%
	mutate(threshold = as.character(threshold))

tissue_levels <- unique(new_isoforms %>% pull(tissue))

# set the colors for the tissues 
colorLevels <- setNames(c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"), levels(as.factor(tissue_levels)))


# new from known plots
pdf("../../figures/M5_our_new_isoforms_in_GTEx/gtex_new_from_known_plots.pdf", 10,10,onefile=TRUE)
# create the number of tissues plot with 4 different CPM thresholds
ggplot(new_isoforms %>% 
       filter(status == 'nfk') %>%
       select(transcript_id, gene_id, n_tissues, status, threshold) %>%
       distinct(), 
	aes(x = n_tissues, fill = threshold)) +
	geom_bar(position='dodge') +
#	theme(legend.position = 'none') + 
	xlab("Number of tissues in which a new RNA isoform from a known gene was seen") +
	ylab("Number of new RNA isoforms")

# create the number of tissues plot with just CPM > 1
ggplot(new_isoforms %>% 
       		filter(threshold == 1.01) %>%
       		filter(status == 'nfk') %>% 
       		select(transcript_id, gene_id, n_tissues, status) %>%
	        distinct(), 
	aes(x = n_tissues)) +
	geom_bar(position='dodge', fill = '#789EBF') +
	theme(legend.position = 'none') + 
	xlab("Number of tissues in which a new RNA isoform from a known gene was seen") +
	ylab("Number of new RNA isoforms")

# when isoform is only expressed in a single tissue, what tissues are the ones expressed?
ggplot(new_isoforms %>% 
       		filter(threshold == 1.01) %>%
       		filter(status == 'nfk') %>% 
		filter(n_tissues == 1) %>%
		drop_na(median_CPM) %>%
       		select(transcript_id, gene_id, tissue) %>%
	        distinct(), 
       	aes(x = tissue, fill = tissue)) +
	geom_bar() +
	theme(
	  axis.text.x = element_blank(),
	  legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx tissue") +
	ylab("Number of new RNA isoforms from a known gene, (median CPM > 1)") +
	labs(fill=NULL)

# display the proportion of new isoforms expressed
proportion_tissue_nfk <- new_isoforms %>%
	filter(threshold == 1.01) %>%
	filter(status == 'nfk') %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfk + all_nfmito))
	
	
ggplot(proportion_tissue_nfk, aes(x=tissue, y=proportion, fill = tissue)) +
	geom_bar(stat='identity') +
	geom_text(aes(label = round(proportion,2)), vjust = -0.5) + 
	theme(
          axis.text.x = element_blank(),
          legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx Tissue") +
	ylab("Proportion of new RNA isoforms from known genes validated (median CPM > 1)") + 
	labs(fill = NULL)
dev.off()

# new from new plots
pdf("../../figures/M5_our_new_isoforms_in_GTEx/gtex_new_from_new_plots.pdf", 10,10,onefile=TRUE)
# create the number of tissues plot with 4 different CPM thresholds
ggplot(new_isoforms %>% 
       filter(status == 'nfn') %>%
       select(transcript_id, gene_id, n_tissues, status, threshold) %>%
       distinct(), 
	aes(x = n_tissues, fill = threshold)) +
	geom_bar(position='dodge') +
#	theme(legend.position = 'none') + 
	xlab("Number of tissues in which a new RNA isoform from a new gene body was seen") +
	ylab("Number of new RNA isoforms")

# create the number of tissues plot with just CPM > 1
ggplot(new_isoforms %>% 
       		filter(threshold == 1.01) %>%
       		filter(status == 'nfn') %>% 
       		select(transcript_id, gene_id, n_tissues, status) %>%
	        distinct(), 
	aes(x = n_tissues)) +
	geom_bar(position='dodge', fill = '#960019') +
	theme(legend.position = 'none') + 
	xlab("Number of tissues in which a new RNA isoform from a new gene body was seen") +
	ylab("Number of new RNA isoforms")

# when isoform is only expressed in a single tissue, what tissues are the ones expressed?
ggplot(new_isoforms %>% 
       		filter(threshold == 1.01) %>%
       		filter(status == 'nfn') %>% 
		filter(n_tissues == 1) %>%
		drop_na(median_CPM) %>%
       		select(transcript_id, gene_id, tissue) %>%
	        distinct(), 
       	aes(x = tissue, fill = tissue)) +
	geom_bar() +
	theme(
	  axis.text.x = element_blank(),
	  legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx tissue") +
	ylab("Number of new RNA isoforms from a new gene body, (median CPM > 1)") +
	labs(fill=NULL)

# display the proportion of new isoforms expressed
proportion_tissue_nfn <- new_isoforms %>%
	filter(threshold == 1.01) %>%
	filter(status == 'nfn') %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfn))
	
	
ggplot(proportion_tissue_nfn, aes(x=tissue, y=proportion, fill = tissue)) +
	geom_bar(stat='identity') +
	geom_text(aes(label = round(proportion,2)), vjust = -0.5) + 
	theme(
          axis.text.x = element_blank(),
          legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx Tissue") +
	ylab("Proportion of new RNA isoforms from new gene bodies validated (median CPM > 1)") + 
	labs(fill = NULL)

dev.off()

# new from medically relevant plots
pdf("../../figures/M5_our_new_isoforms_in_GTEx/gtex_new_from_medically_relevant_plots.pdf", 10,10,onefile=TRUE)
# create the number of tissues plot with 4 different CPM thresholds
ggplot(new_isoforms %>% 
       filter(is_med_rel == TRUE) %>%
       select(transcript_id, gene_id, n_tissues, status, threshold) %>%
       distinct(), 
	aes(x = n_tissues, fill = threshold)) +
	geom_bar(position='dodge') +
#	theme(legend.position = 'none') + 
	xlab("Number of tissues in which a new RNA isoform from a medically relevant gene was seen") +
	ylab("Number of new RNA isoforms")

# create the number of tissues plot with just CPM > 1
ggplot(new_isoforms %>% 
       		filter(threshold == 1.01) %>%
       		filter(is_med_rel == TRUE) %>% 
       		select(transcript_id, gene_id, n_tissues, status) %>%
	        distinct(), 
	aes(x = n_tissues)) +
	geom_bar(position='dodge', fill = '#3D550C') +
	theme(legend.position = 'none') + 
	xlab("Number of tissues in which a new RNA isoform from a medically relevant gene was seen") +
	ylab("Number of new RNA isoforms")

# when isoform is only expressed in a single tissue, what tissues are the ones expressed?
ggplot(new_isoforms %>% 
       		filter(threshold == 1.01) %>%
       		filter(is_med_rel == TRUE) %>% 
		filter(n_tissues == 1) %>%
		drop_na(median_CPM) %>%
       		select(transcript_id, gene_id, tissue) %>%
	        distinct(), 
       	aes(x = tissue, fill = tissue)) +
	geom_bar() +
	theme(
	  axis.text.x = element_blank(),
	  legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx tissue") +
	ylab("Number of new RNA isoforms from a medically relevant gene, (median CPM > 1)") +
	labs(fill=NULL)

# display the proportion of new isoforms expressed
proportion_tissue_nfmr <- new_isoforms %>%
	filter(threshold == 1.01) %>%
	filter(is_med_rel == TRUE) %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfmed_rel))
	
	
ggplot(proportion_tissue_nfmr, aes(x=tissue, y=proportion, fill = tissue)) +
	geom_bar(stat='identity') +
	geom_text(aes(label = round(proportion,2)), vjust = -0.5) + 
	theme(
          axis.text.x = element_blank(),
          legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx Tissue") +
	ylab("Proportion of new RNA isoforms from medically relevant genes validated (median CPM > 1)") + 
	labs(fill = NULL)

dev.off()

q()
#TODO: SKIPPING FOR NOW
pdf('../../figures/M5_our_new_isoforms_in_GTEx/gtex_median_cpm_per_tissue.pdf', 15, 15, onefile=TRUE)
#pdf('../../figures/R/M4_new_tx_from_our_data_in_GTEx.pdf ', onefile=TRUE)

# look at the expression of all the new (bambu) transcripts across the tissues
ggplot(combined %>% mutate(log10_med_cpm = log10(median_cpm)), aes(x = tissue_site_detail, y = log10_med_cpm, fill=tissue_site_detail)) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = colorLevels) +
	geom_jitter(height = 0, width = 0.15)

ggplot(combined %>% filter(median_cpm > 1) %>% mutate(log10_med_cpm = log10(median_cpm)), aes(x = tissue_site_detail, y = log10_med_cpm, fill=tissue_site_detail)) +
	geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = colorLevels) +
	geom_jitter(height = 0, width = 0.15)

dev.off()
