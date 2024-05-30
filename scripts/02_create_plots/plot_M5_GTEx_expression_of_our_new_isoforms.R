library(tidyverse)
set.seed(21)

# number of isoforms in each category
all_nfk = 428
all_nfn = 267
all_nfmito = 5
all_nfmed_rel = 54 

h=10
w=11

# load in data
new_isoforms <- read_tsv("../../tables/GTEx_expression_our_new_isoforms/GTEx_expression_of_new_isoforms.tsv") %>%
	mutate(n_tissues = as.character(n_tissues)) %>%
	mutate(threshold = as.character(threshold))

tissue_levels <- unique(new_isoforms %>% pull(tissue))

# set the colors for the tissues 
colorLevels <- setNames(c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3"), levels(as.factor(tissue_levels)))


# new from known plots
pdf("../../figures/M5_our_new_isoforms_in_GTEx/gtex_new_from_known_plots.pdf", w,h,onefile=TRUE)
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
	geom_text(aes(label = after_stat(count) ),stat='count', vjust = -0.5) + 
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
	geom_text(aes(label = after_stat(count) ),stat='count', vjust = -0.5) + 
	theme(
	  axis.text.x = element_blank(),
	  legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx tissue") +
	ylab("Number of new RNA isoforms from a known gene, (median CPM > 1)") +
	labs(fill=NULL)

# display the proportion of new isoforms expressed
proportion_tissue_nfk_gt_0 <- new_isoforms %>%
	filter(threshold == 0.01) %>%
	filter(status == 'nfk') %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfk + all_nfmito))
	
proportion_tissue_nfk_gt_1 <- new_isoforms %>%
	filter(threshold == 1.01) %>%
	filter(status == 'nfk') %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfk + all_nfmito))
	
	
ggplot() +
	geom_bar(proportion_tissue_nfk_gt_0, mapping=aes(x=tissue, y=proportion, fill=tissue), stat='identity', alpha = 0.25) +
	geom_text(proportion_tissue_nfk_gt_0, mapping=aes(x=tissue, y=proportion, label = round(proportion,2)), alpha = 0.25, vjust = -0.5) + 
	geom_bar(proportion_tissue_nfk_gt_1, mapping=aes(x=tissue, y=proportion, fill=tissue), stat='identity') +
	geom_text(proportion_tissue_nfk_gt_1, mapping=aes(x=tissue, y=proportion, label = round(proportion,2)), vjust = -0.5) + 
	theme(
          axis.text.x = element_blank(),
          legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx Tissue") +
	ylab("Proportion of new RNA isoforms from known genes validated") + 
	labs(fill = NULL)
dev.off()

# new from new plots
pdf("../../figures/M5_our_new_isoforms_in_GTEx/gtex_new_from_new_plots.pdf", w,h,onefile=TRUE)
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
	geom_text(aes(label = after_stat(count) ),stat='count', vjust = -0.5) + 
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
	geom_text(aes(label = after_stat(count) ),stat='count', vjust = -0.5) + 
	theme(
	  axis.text.x = element_blank(),
	  legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx tissue") +
	ylab("Number of new RNA isoforms from a new gene body, (median CPM > 1)") +
	labs(fill=NULL)

# display the proportion of new isoforms expressed
proportion_tissue_nfn_gt_0 <- new_isoforms %>%
	filter(threshold == 0.01) %>%
	filter(status == 'nfn') %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfn))
	
proportion_tissue_nfn_gt_1 <- new_isoforms %>%
	filter(threshold == 1.01) %>%
	filter(status == 'nfn') %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfn))
	
	
ggplot() +
	geom_bar(proportion_tissue_nfn_gt_0, mapping=aes(x=tissue, y=proportion, fill=tissue), stat='identity', alpha = 0.25) +
	geom_text(proportion_tissue_nfn_gt_0, mapping=aes(x=tissue, y=proportion, fill=tissue, label = round(proportion,2)), alpha = 0.25, vjust = -0.5) + 
	geom_bar(proportion_tissue_nfn_gt_1, mapping=aes(x=tissue, y=proportion, fill=tissue), stat='identity') +
	geom_text(proportion_tissue_nfn_gt_1, mapping=aes(x=tissue, y=proportion, fill=tissue, label = round(proportion,2)), vjust = -0.5) + 
	geom_bar(stat='identity') +
	theme(
          axis.text.x = element_blank(),
          legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx Tissue") +
	ylab("Proportion of new RNA isoforms from new gene bodies validated (median CPM > 1)") + 
	labs(fill = NULL)

dev.off()

# new from medically relevant plots
pdf("../../figures/M5_our_new_isoforms_in_GTEx/gtex_new_from_medically_relevant_plots.pdf", w,h,onefile=TRUE)
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
	geom_text(aes(label = after_stat(count) ),stat='count', vjust = -0.5) + 
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
	geom_text(aes(label = after_stat(count) ),stat='count', vjust = -0.5) + 
	theme(
	  axis.text.x = element_blank(),
	  legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx tissue") +
	ylab("Number of new RNA isoforms from a medically relevant gene, (median CPM > 1)") +
	labs(fill=NULL)

# display the proportion of new isoforms expressed
proportion_tissue_nfmr_gt_0 <- new_isoforms %>%
	filter(threshold == 0.01) %>%
	filter(is_med_rel == TRUE) %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfmed_rel))
	
proportion_tissue_nfmr_gt_1 <- new_isoforms %>%
	filter(threshold == 1.01) %>%
	filter(is_med_rel == TRUE) %>%
	drop_na(median_CPM) %>%
	select(transcript_id, gene_id, tissue) %>%
	distinct() %>%
	group_by(tissue) %>%
	summarize(n_isoforms = n()) %>%
	rowwise() %>%
	mutate(proportion = n_isoforms / (all_nfmed_rel))
	
ggplot() +
	geom_bar(proportion_tissue_nfmr_gt_0, mapping=aes(x=tissue, y=proportion, fill=tissue), stat='identity', alpha = 0.25) +
	geom_text(proportion_tissue_nfmr_gt_0, mapping=aes(x=tissue, y=proportion, fill=tissue, label = round(proportion,2)), alpha = 0.25, vjust = -0.5) + 
	geom_bar(proportion_tissue_nfmr_gt_1, mapping=aes(x=tissue, y=proportion, fill=tissue), stat='identity') +
	geom_text(proportion_tissue_nfmr_gt_1, mapping=aes(x=tissue, y=proportion, fill=tissue, label = round(proportion,2)), vjust = -0.5) + 
	geom_bar(stat='identity') +
	theme(
          axis.text.x = element_blank(),
          legend.position = 'none') +
	scale_fill_manual(values = colorLevels) +
	xlab("GTEx Tissue") +
	ylab("Proportion of new RNA isoforms from medically relevant genes validated (median CPM > 1)") + 
	labs(fill = NULL)

dev.off()
# look at expression of isoforms from new genes across tissues
new_genes <- new_isoforms %>%
	filter(threshold == 0.01) %>%
	filter(grepl('Bambu', gene_id)) #%>%
#	select(!c(gene_biotype, gene_name, threshold)) %>%
#	pivot_longer(!c(transcript_id, gene_id), names_to='tissue', values_to='median_CPM') %>%
#	replace(is.na(.), 0)

#ggplot(new_genes, aes(x=tissue, y=median_CPM, fill = tissue)) +
#	geom_boxplot(outlier.shape=NA)+
#	geom_jitter(height=0, width = 0.15, alpha = 0.2) +
#	geom_text(aes(label=ifelse((median_CPM > 10), gene_id, '')), size=3, vjust=-0.2) +
#	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
#	theme(axis.text.x = element_text(angle = 45, hjust=1),
#	      legend.position='none')
#ggsave('../../figures/M5_our_new_isoforms_in_GTEx/median_CPM_of_isoforms_from_new_genes.pdf')

new_gene_mass_spec <- c('BambuTx1361',
			'BambuTx2524',
			'BambuTx1697')
n_iso_above_5_cpm <- length(unique(new_genes %>%
      filter(median_CPM > 5) %>%
      pull(transcript_id)))
n_iso_above_1_cpm <- length(unique(new_genes %>%
				   filter(median_CPM > 1) %>%
				   pull(transcript_id)))

print("Isoforms from new genes above 5 cpm in at least 1 tissue")
print(n_iso_above_5_cpm)
print(n_iso_above_5_cpm / all_nfn)
print("Isoforms from new genes above 1 cpm in at least 1 tissue")
print(n_iso_above_1_cpm)
print(n_iso_above_1_cpm / all_nfn)

new_genes <- new_genes %>% filter(threshold == 0.01)
ggplot(new_genes %>% mutate(med_cpm_log2 = log2(median_CPM +1)), aes(x=tissue, y=med_cpm_log2, fill = tissue)) +
	geom_boxplot(outlier.shape=NA)+
	geom_jitter(height=0, width = 0.15, aes(
                                                color=ifelse(transcript_id == 'BambuTx1361', 'red', 'black'),
                                                alpha=ifelse(transcript_id == 'BambuTx1361', 0.9, 0.2))) +
	scale_color_identity() +
        geom_text(aes(label=ifelse((median_CPM > 10), gene_id, '')), size=3, vjust=-0.2) +
        geom_text(aes(label=ifelse(transcript_id == 'BambuTx1361', paste0(gene_id,": ", transcript_id), '')), size=2, vjust=-0.2) +

	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	geom_hline(yintercept=log2(2), linetype='dotted', col = 'darkred') +
	#theme(axis.text.x = element_text(angle = 45, hjust=1),
	theme(axis.text.x = element_blank(),
	      legend.position='none')
ggsave('../../figures/M5_our_new_isoforms_in_GTEx/median_log2_CPM+1_of_isoforms_from_new_genes.pdf', width=w, height=h)

med_rel_genes <- read_tsv("../../references/medically_relevant_genes.tsv")
med_rel_ids <- med_rel_genes %>% pull(gene_id)

ids_of_tx_for_mass_spec <- c('BambuTx1646', 
			     'BambuTx370', 
			     'BambuTx572', 
			     'BambuTx559', 
			     'BambuTx1879', 
			     'BambuTx1188',
			     'BambuTx1164',
			     'BambuTx703',
			     'BambuTx17',
			     'BambuTx1758',
			     'BambuTx2189')

genes_CPM <- new_isoforms %>%
        filter(grepl('Bambu', transcript_id)) %>%
        filter(threshold == 0.01)


ggplot(genes_CPM %>% mutate(med_cpm_log2 = log2(median_CPM +1)), aes(x=tissue, y=med_cpm_log2, fill = tissue)) +
        geom_boxplot(outlier.shape=NA)+
        geom_jitter(height=0, width = 0.15, aes(
                                                color=ifelse(transcript_id %in% ids_of_tx_for_mass_spec, 'red','black'),
                                                alpha=ifelse(transcript_id %in% ids_of_tx_for_mass_spec, 0.9, 0.2))) +
        scale_color_identity() +
        geom_text(aes(label=ifelse(transcript_id %in% ids_of_tx_for_mass_spec, paste0(gene_name,": ", transcript_id), '')), size=2, vjust=-0.2) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
        geom_hline(yintercept=log2(2),linetype='dotted', col = 'darkred') +
        #theme(axis.text.x = element_text(angle = 45, hjust=1),
        theme(axis.text.x = element_blank(),
              legend.position='none')

ggsave('../../figures/M5_our_new_isoforms_in_GTEx/median_log2_CPM+1_of_isoforms_from_genes.pdf', width=w, height=h)







med_rel_genes_CPM <- new_isoforms %>%
	filter(gene_id %in% med_rel_ids) %>%
	filter(grepl('Bambu', transcript_id)) %>%
#	select(-threshold) %>%
#	pivot_longer(!c(transcript_id, gene_id, gene_name, gene_biotype), names_to='tissue', values_to='median_CPM') %>%
#	replace(is.na(.), 0)
	filter(threshold == 0.01)

#ggplot(med_rel_genes_CPM, aes(x=tissue, y=median_CPM, fill = tissue)) +
#        geom_boxplot(outlier.shape=NA)+
#        geom_jitter(height=0, width = 0.15, alpha = 0.2) +
#        #geom_text(aes(label=ifelse((median_CPM > 20), gene_name, '')), size=3, vjust=-0.2) +
#        geom_text(aes(label=ifelse(transcript_id %in% ids_of_tx_for_mass_spec, gene_name, '')), size=2, vjust=-0.2) +
#        scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
#        theme(axis.text.x = element_text(angle = 45, hjust=1),
#              legend.position='none')
#ggsave('../../figures/M5_our_new_isoforms_in_GTEx/median_CPM_of_isoforms_from_medically_relevant_genes.pdf')

ggplot(med_rel_genes_CPM %>% mutate(med_cpm_log2 = log2(median_CPM +1)), aes(x=tissue, y=med_cpm_log2, fill = tissue)) +
        geom_boxplot(outlier.shape=NA)+
        geom_jitter(height=0, width = 0.15, aes(
						color=ifelse(transcript_id == 'BambuTx1879', 'blue', ifelse(transcript_id == 'BambuTx703', 'yellow' ,'black')),
						alpha=ifelse(transcript_id %in% ids_of_tx_for_mass_spec, 0.9, 0.2))) +
	scale_color_identity() +
        #geom_text(aes(label=ifelse((median_CPM > 20), gene_name, '')), size=3, vjust=-0.2) +
        geom_text(aes(label=ifelse(transcript_id %in% ids_of_tx_for_mass_spec, paste0(gene_name,": ", transcript_id), '')), size=2, vjust=-0.2) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	geom_hline(yintercept=log2(2),linetype='dotted', col = 'darkred') +
        #theme(axis.text.x = element_text(angle = 45, hjust=1),
        theme(axis.text.x = element_blank(),
              legend.position='none')

ggsave('../../figures/M5_our_new_isoforms_in_GTEx/median_log2_CPM+1_of_isoforms_from_medically_relevant_genes.pdf', width=w, height=h)

## Preferential expression

p_expr <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_1.01_preferential_filtered.tsv')

ggplot(p_expr, aes(x=n_tis_pref)) +
	geom_bar()


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
