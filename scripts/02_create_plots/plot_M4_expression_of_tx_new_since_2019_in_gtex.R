library(tidyverse)
library(scales)
library(biomaRt)
library(ggupset)

tissues_to_use = c("Brain - Cerebellar Hemisphere",
                   "Brain - Frontal Cortex (BA9)",
                   "Brain - Putamen (basal ganglia)",
                   "Cells - Cultured fibroblasts",
                   "Heart - Atrial Appendage",
                   "Heart - Left Ventricle",
                   "Liver",
                   "Lung",
                   "Muscle - Skeletal")

##########################################################################################
#### Looking at just isoforms that are new since 2019 ############

# load the numbers
exclusive_values_at_cpm_threshold_gtex <- read_tsv("../../tables/GTEx_expression_since_2019/gtex_values_at_cpm_thresholds_new_since_2019.tsv")

## plot the number of isoforms new since 2019 across different CPM thresholds
ggplot(exclusive_values_at_cpm_threshold_gtex, aes(x = cpm_threshold, y = all_transcript_median, color=tissue)) +
        geom_line(linewidth=0.25) +
        geom_vline(aes(xintercept = 1), color = '#4d4d4d', linetype='dashed', linewidth=0.25) + #, show.legend=FALSE) +
	annotate("text", label='CPM = 1', x =1.75, y = 4000, size=1) +
        xlab("Median CPM Threshold") +
        ylab("Number of new isoforms since 2019") +
        scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10)) +
	theme(legend.position = c(0.75, 0.75),
	      text = element_text(size=5),
	      axis.title = element_text(size=5),
	      legend.key.size = unit(1, 'mm'),
	      legend.title = element_text(size=4),
	      legend.text = element_text(size=3))

ggsave('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_cpm_thresholds.pdf', width=58, height=58, units='mm')


# plot the number of protein_coding isoforms WITH NEW CDS SEQUENCE across 4 different thresholds
ggplot() +
        geom_bar(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 0.01), mapping=aes(x=tissue, y=new_cds_transcript_median), fill = '#00BFC4', stat='identity') +
        geom_text(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 0.01), mapping=aes(x=tissue, y=new_cds_transcript_median, label = new_cds_transcript_median), vjust = -0.2, size=1) +
        geom_bar(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 1.01), mapping=aes(x=tissue, y=new_cds_transcript_median), fill = '#F8766D', stat='identity') +
        geom_text(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 1.01), mapping=aes(x=tissue, y=new_cds_transcript_median, label = new_cds_transcript_median), vjust = -0.2, size=1) +
        geom_bar(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 5.01), mapping=aes(x=tissue, y=new_cds_transcript_median), fill = '#7CAE00', stat='identity') +
        geom_text(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 5.01), mapping=aes(x=tissue, y=new_cds_transcript_median, label = new_cds_transcript_median), vjust = -0.2, size=1) +
        geom_bar(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 10.01), mapping=aes(x=tissue, y=new_cds_transcript_median), fill = '#C77CFF', stat='identity') +
        geom_text(exclusive_values_at_cpm_threshold_gtex %>% filter(cpm_threshold == 10.01), mapping=aes(x=tissue, y=new_cds_transcript_median, label = new_cds_transcript_median), vjust = -0.2, size=1) +
	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text=element_text(size=3),
              text = element_text(size=4),
              legend.title = element_text(size=3),
              legend.key.size=unit(1, 'mm'),
              legend.margin=margin(0,0)) +
        xlab('GTEx tissue') +
        ylab('Number of new protein coding isoforms with new cds sequence since 2019')
ggsave('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_cpm_thresholds_pc_bar.pdf', width=58, height=58, units='mm')


# plot density of CPM for all tissues for reads that are new since 2019
exclusive_tx_cpm <- read_tsv('../../tables/GTEx_expression_since_2019/gtex_values_cpm_tx_new_since_2019.tsv') %>%
	filter(median_CPM > 0) %>%
	mutate(log2_cpm = log2(median_CPM))

ggplot(exclusive_tx_cpm, aes(x=log2_cpm, color = tissue_site_detail)) +
	geom_density(linewidth=0.25) +
	geom_vline(aes(xintercept = 0), color='#4d4d4d', linetype='dashed', linewidth=0.25) +
	xlab("Log2 median CPM (only including transcripts with median CPM > 0)") +
        ylab("Number of new isoforms since 2019") +
	theme(axis.text = element_text(size=3),
              text = element_text(size=4),
              legend.title = element_text(size=3),
              legend.key.size=unit(1, 'mm'),
              legend.margin=margin(1.5,1.5,1.5,1.5),
              legend.position=c(0.75,0.75))

ggsave('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_density_plot.pdf', width=58, height=58, units='mm')


# barplots showing the number of isoforms in each category expressed in each tissue
bar_cat <- c("Med-relevant", "New protein coding sequence", "Med-relevant new protein coding sequence", "Brain disease relevant", "Brain disease relevant new protein coding sequence" )

barplot_df <- exclusive_values_at_cpm_threshold_gtex %>% 
	filter(cpm_threshold == 1.01) %>%
	rename(`Med-relevant` = med_relevant_transcript_median) %>%
	rename(`New protein coding sequence` = new_cds_transcript_median) %>%
	rename(`Med-relevant new protein coding sequence` = med_relevant_new_cds_transcript_median) %>%
	rename(`Brain disease relevant` = brain_relevant_transcript_median) %>%
	rename(`Brain disease relevant new protein coding sequence` = brain_relevant_new_cds_transcript_median) %>%
	pivot_longer(!c(cpm_threshold, tissue), names_to = 'category', values_to = 'n_isoforms') %>%
	filter(category %in% bar_cat) %>%
	mutate(category = factor(category, levels = bar_cat))
	

ggplot(barplot_df, aes(x=category, y=n_isoforms, fill=tissue)) +
	geom_bar(stat='identity', position='dodge') +
	geom_text(aes(label = n_isoforms, group=tissue), position=position_dodge(width=0.9), vjust= -0.2, size=1.5) + 
	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
#	theme(axis.text.x = element_text(angle = 45, hjust=1)) +
	theme(legend.position = 'none',
	      text = element_text(size=5)) +
	ylab("Number of new isoforms since 2019")
ggsave('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_by_category.pdf', width=118, height=58, units='mm')


# plot showing the relative abundance of new isoforms of the genes in each category
perc_expr_df <- read_tsv('../../tables/GTEx_expression_since_2019/gtex_cpm_gt_1_relative_abundance_new_since_2019.tsv') %>%
	mutate(Category = factor(Category, levels = bar_cat)) %>%
	rename(rel_abund = `Relative Abundance of New Transcripts Since 2019 (%)`)
	
ggplot(perc_expr_df, aes(x=Tissue, y=rel_abund, fill=Tissue)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter(height=0, width = 0.15, alpha = 0.2) +
	theme(legend.position = c(0.85, 0.30), axis.text.x = element_blank()) +
	facet_wrap(~Category) 
#	ggtitle(cat)
ggsave('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_relative_abundance.pdf', width=11, height=8.5, units='in')

pdf('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_relative_abundance_by_category.pdf', width=58/25.4, height=58/25.4)
for (cat in bar_cat) {
	print(ggplot(perc_expr_df %>% filter(Category == cat), aes(x=Tissue, y=rel_abund, fill=Tissue)) +
		geom_boxplot(outlier.shape = NA, linewidth=0.25) +
		geom_jitter(height=0, width = 0.15, alpha = 0.2, size=0.25) +
		ggtitle(cat) + 
		scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
		theme(axis.text.x = element_text(angle = 45, hjust=1), 
		      axis.text = element_text(size=3),
 	              text = element_text(size=4),
        	      legend.title = element_text(size=3),
	              legend.key.size=unit(1, 'mm'),
	              legend.margin=margin(1.5,1.5,1.5,1.5),
		      legend.position='none')
 )
}

dev.off()

tiss_color <- setNames(hue_pal()(length(tissues_to_use)), levels(as.factor(tissues_to_use)))

####poster children ###########################################################
rel_abund_poster_children <- perc_expr_df

pdf('../../figures/M4_new_tx_since_2019_in_GTEx/M4_new_since_2019_in_gtex_relative_abundance_poster_children.pdf', width=58/25.4, height=58/25.4)

ggplot(rel_abund_poster_children %>% 
       filter(Category == 'Med-relevant new protein coding sequence') %>%
       filter(gene_name == "HSPA9"), aes(x = rel_abund, y = Tissue, fill = Tissue)) +
	geom_bar(stat = 'identity') +
	scale_fill_manual(values = tiss_color) +
	theme(legend.position = 'none',
 	              text = element_text(size=4),
	      axis.text = element_text(size=3)) +
	ggtitle("HSPA9")
	
ggplot(rel_abund_poster_children %>% 
       filter(Category == 'Med-relevant new protein coding sequence') %>%
       filter(gene_name == "DYM"), aes(x = rel_abund, y = Tissue, fill = Tissue)) +
	geom_bar(stat = 'identity') +
	scale_fill_manual(values = tiss_color) +
	theme(legend.position = 'none',
 	              text = element_text(size=4),
	      axis.text = element_text(size=3)) +
	ggtitle("DYM")
	
ggplot(rel_abund_poster_children %>% 
       filter(Category == 'Brain disease relevant new protein coding sequence') %>%
       filter(gene_name == "PAM"), aes(x = rel_abund, y = Tissue, fill = Tissue)) +
	geom_bar(stat = 'identity') +
	scale_fill_manual(values = tiss_color) +
	theme(legend.position = 'none',
 	              text = element_text(size=4),
	      axis.text = element_text(size=3)) +
	ggtitle("PAM")
	
ggplot(rel_abund_poster_children %>% 
       filter(Category == 'Med-relevant new protein coding sequence') %>%
       filter(gene_name == "KIF5A"), aes(x = rel_abund, y = Tissue, fill = Tissue)) +
	geom_bar(stat = 'identity') +
	scale_fill_manual(values = tiss_color) +
	theme(legend.position = 'none',
 	              text = element_text(size=4),
	      axis.text = element_text(size=3)) +
	ggtitle("KIF5A")

dev.off()


