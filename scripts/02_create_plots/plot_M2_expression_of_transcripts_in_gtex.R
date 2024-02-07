library(tidyverse)
library(scales)
#library(biomaRt)
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
###### Module 2 - TX IN GTEX DATA ########################################
##########################################################################################

# plot the number of isoforms across CPM thresholds for each tissue
gtex_CPM_thresholds <- read_tsv('../../tables/GTEx_expression/gtex_values_at_cpm_thresholds.tsv')

ggplot(gtex_CPM_thresholds, aes(x=cpm_threshold, y=all_transcript_median, color=tissue))+
	geom_line(linewidth=.25) +
	geom_vline(aes(xintercept = 1), color = '#4d4d4d', linetype='dashed', linewidth=.25) +
	annotate("text", label='CPM = 1', x =1.75, y = 30000, size=1) +
	theme(legend.position = c(0.75, 0.75),
	      text = element_text(size=5),
	      axis.title = element_text(size=5),
	      legend.key.size = unit(1, 'mm'),
	      legend.title = element_text(size=4),
	      legend.text = element_text(size=3)) +
	xlab('CPM thresholds') +
	ylab('Number of isoforms expressed at the threshold')
	
ggsave('../../figures/M2_GTEx_expression/M2_gtex_cpm_thresholds.pdf', width = 58, height = 58, units = 'mm')

##########################################################################################
# create a stacked bar plot of the protein coding isoforms in each tissue by CPM threshold. Because we 
# rounded the CPM's to 2 decimal places, taking for example CPM threshold == 1.01 is essentially saying
# median CPM > 1
ggplot() +
        geom_bar(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 0.01),
                 mapping = aes(x = tissue, y = cds_transcript_median, fill = "> 0"), stat = 'identity') +
        geom_text(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 0.01),
                  mapping = aes(x = tissue, y = cds_transcript_median, label = cds_transcript_median), vjust = -0.2, size=1) +

        geom_bar(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 1.01),
                 mapping = aes(x = tissue, y = cds_transcript_median, fill = "> 1"), stat = 'identity') +
        geom_text(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 1.01),
                  mapping = aes(x = tissue, y = cds_transcript_median, label = cds_transcript_median), vjust = -0.2, size=1) +

        geom_bar(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 5.01),
                 mapping = aes(x = tissue, y = cds_transcript_median, fill = "> 5"), stat = 'identity') +
        geom_text(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 5.01),
                  mapping = aes(x = tissue, y = cds_transcript_median, label = cds_transcript_median), vjust = -0.2, size=1) +

        geom_bar(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 10.01),
                 mapping = aes(x = tissue, y = cds_transcript_median, fill = "> 10"), stat = 'identity') +
        geom_text(data = gtex_CPM_thresholds %>% filter(cpm_threshold == 10.01),
                  mapping = aes(x = tissue, y = cds_transcript_median, label = cds_transcript_median), vjust = -0.2, size=1) +
        scale_fill_manual(values = c("> 0" = '#00BFC4', "> 1" = '#F8766D', "> 5" = '#7CAE00', "> 10" = '#C77CFF'),
			  breaks = c("> 0", "> 1", "> 5", "> 10"),
                          name = "CPM Threshold") +

        scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
	      axis.text=element_text(size=3),
	      text = element_text(size=4),
	      legend.title = element_text(size=3),
	      legend.key.size=unit(1, 'mm'),
	      legend.margin=margin(0,0)) +
        xlab('GTEx tissue') +
        ylab('Number of protein coding transcripts expressed at different thresholds')

ggsave('../../figures/M2_GTEx_expression/M2_gtex_cpm_thresholds_pc_bar.pdf', width=58, height=58, units='mm')

##########################################################################################

# Bar plot that showcases the difference between the number of isoforms vs protein-coding isoforms at CPM > 1
ggplot() +
	geom_bar(gtex_CPM_thresholds %>% filter(cpm_threshold == 1.01), mapping=aes(x=tissue, y=all_transcript_median, fill=tissue), stat = 'identity') +
	geom_text(gtex_CPM_thresholds %>% filter(cpm_threshold == 1.01), mapping=aes(x=tissue, y=all_transcript_median, label = all_transcript_median), vjust = -0.2, size=1) +
	geom_bar(gtex_CPM_thresholds %>% filter(cpm_threshold == 1.01), mapping=aes(x=tissue, y=cds_transcript_median), fill = '#4d4d4d', stat='identity', alpha = 0.5) +
	geom_text(gtex_CPM_thresholds %>% filter(cpm_threshold == 1.01), mapping=aes(x=tissue, y=cds_transcript_median, label = cds_transcript_median), vjust = -0.2, size=1) +
	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	theme(axis.text.x = element_text(angle = 45, hjust=1),
	      axis.text = element_text(size=3),
	      text = element_text(size=4)) +
	theme(legend.position = 'none') +
	xlab('GTEx tissue')  +
	ylab('Number of isoforms')

ggsave('../../figures/M2_GTEx_expression/M2_gtex_cpm_thresholds_all_v_pc_values.pdf', width=58, height=58, units='mm')

##########################################################################################

# bar plot that shows the propotion of protein-coding isoforms vs non protein-coding isoforms
all_v_pc_tmp <- gtex_CPM_thresholds %>% 
	filter(cpm_threshold == 1.01) %>%
	dplyr::select(tissue, all_transcript_median, cds_transcript_median) %>%
	rowwise() %>%
	mutate(leftover = all_transcript_median - cds_transcript_median) %>%
	dplyr::select(tissue, cds_transcript_median, leftover) %>%
	pivot_longer(!tissue, names_to = 'type', values_to = 'n_tx') %>%
	mutate(type = factor(type, levels = c("cds_transcript_median", "leftover")))

# Reorder the levels of 'type' factor
all_v_pc_tmp$type <- factor(all_v_pc_tmp$type, levels = c("cds_transcript_median", "leftover"))

ggplot(all_v_pc_tmp, aes(x=tissue, y=n_tx, fill=type)) +
	geom_bar(position = 'fill', stat = 'identity') +
#	theme(legend.position = 'top') +
	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	scale_fill_manual(
            values = c("cds_transcript_median" = '#18bdc2', "leftover" = '#f3766e'),
            #values = c("leftover" = '#f3766e', "cds_transcript_median" = '#18bdc2'),
            labels = c("Protein-coding isoforms", "Other isoforms")
        ) +
	theme(axis.text.x = element_text(angle = 45, hjust=1),
	      axis.text = element_text(size=3),
              text = element_text(size=4), 
	      legend.title = element_text(size=3),
              legend.key.size=unit(1, 'mm'),
              legend.margin=margin(1.5,1.5,1.5,1.5),
	      legend.position=c(0.77,0.15)) +
	xlab('GTEx tissue') +
	ylab('Proportion of isoforms by type expressed')

ggsave('../../figures/M2_GTEx_expression/M2_gtex_cpm_thresholds_all_v_pc_proportion.pdf', width=58, height=58, units='mm')

##########################################################################################
##### upset plot #######################################################################
# like a really really big venn diagram
# looking at what isoforms are unique/shared between tissues, only looking at the top 20 categories
iso_to_tiss <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv') %>% 
	select(!c(threshold)) %>%
	pivot_longer(!c(transcript_id, gene_id, gene_name, gene_biotype), names_to = 'tissue', values_to = 'med_CPM') %>% 
	drop_na(med_CPM) %>% 
	group_by(transcript_id) %>%
        summarize(tissue = list(unique(tissue))) %>%
        ungroup()

upset_plt <- ggplot(iso_to_tiss, aes(x=tissue)) + 
	geom_bar() +
#	geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
	scale_x_upset(n_intersections = 20) +
	theme_combmatrix(
			 combmatrix.panel.line.size = 0.5,
			 combmatrix.panel.point.size = 1) +
	xlab('Tissue combinations') +
	ylab('Number of isoforms') +
	theme(axis.text.y = element_text(size=2),
	      axis.title = element_text(size=5))
ggsave('../../figures/M2_GTEx_expression/M2_GTEx_tissue_upset_plot.pdf', upset_plt, height = 58, width = 58 , units='mm')

upset_plt <- ggplot(iso_to_tiss, aes(x=tissue)) + 
	geom_bar() +
	geom_text(stat='count', aes(label=after_stat(count)), vjust=-1, size=2) +
	scale_x_upset(n_intersections = 20) +
	xlab('Tissue combinations') +
	ylab('Number of isoforms') 
ggsave('../../figures/M2_GTEx_expression/M2_GTEx_tissue_upset_plot_with_numbers.pdf', upset_plt)

# upset plot but for genes
gene_to_tiss <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv') %>%
        select(!c(threshold)) %>%
        pivot_longer(!c(transcript_id, gene_id, gene_name, gene_biotype), names_to = 'tissue', values_to = 'med_CPM') %>%
        drop_na(med_CPM) %>%
	select(gene_id, gene_name, gene_biotype, tissue) %>%
        group_by(gene_id) %>%
        summarize(tissue = list(unique(tissue))) %>%
        ungroup()

upset_plt <- ggplot(gene_to_tiss, aes(x=tissue)) +
        geom_bar() +
        geom_text(stat='count', aes(label=after_stat(count)), vjust=-1, size=2) +
        scale_x_upset(n_intersections = 20) +
        xlab('Tissue combinations') +
        ylab('Number of genes')
ggsave('../../figures/M2_GTEx_expression/M2_GTEx_tissue_gene_upset_plot_with_numbers.pdf', upset_plt)


########################################################################################

# plot the isoforms 'unique' to each tissue, and then the number of those that are protein_coding
unique_to_tiss <- read_tsv('../../tables/GTEx_expression/GTEx_tissue_unique_isoforms_by_type_CPM_gt_1.tsv')
ggplot() +
	geom_bar(unique_to_tiss %>% filter(type=='all'), mapping=aes(x=tissue, y=n_tx, fill=tissue), stat = 'identity') +
	geom_text(unique_to_tiss %>% filter(type=='all'), mapping=aes(x=tissue, y=n_tx, label = n_tx), vjust = -0.2, size=1) +
	geom_bar(unique_to_tiss %>% filter(type=='pc'), mapping=aes(x=tissue, y=n_tx), fill = '#4d4d4d', stat='identity', alpha = 0.5) +
	geom_text(unique_to_tiss %>% filter(type =='pc'), mapping=aes(x=tissue, y=n_tx, label = n_tx), vjust = -0.2, size=1) +
	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	theme(axis.text.x = element_text(angle = 45, hjust=1),
              axis.text = element_text(size=3),
              text = element_text(size=4)) +
	theme(legend.position = 'none') +
	xlab('GTEx tissue') +
	ylab('Number of isoforms')

ggsave('../../figures/M2_GTEx_expression/M2_gtex_tx_unique_to_tissue_all_v_pc_values.pdf', width=58, height=58, units='mm')

########################################################################################

# plot the proportion of the unique isoforms that are protein-coding vs non protein-coding
all_v_pc_unique <- unique_to_tiss %>%
	pivot_wider(names_from = 'type', values_from = 'n_tx') %>%
        rowwise() %>%
        mutate(leftover = all - pc) %>%
        dplyr::select(tissue, pc, leftover) %>%
        pivot_longer(!tissue, names_to = 'type', values_to = 'n_tx') %>%
        mutate(type = factor(type, levels = c("pc", "leftover")))


ggplot(all_v_pc_unique, aes(x=tissue, y=n_tx, fill=type)) +
        geom_bar(position = 'fill', stat = 'identity') +
	scale_fill_manual(
            values = c("pc" = '#18bdc2', "leftover" = '#f3766e'),
            labels = c("Protein-coding isoforms", "Other isoforms")
        ) +
	scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
	theme(axis.text.x = element_text(angle = 45, hjust=1),
              axis.text = element_text(size=3),
              text = element_text(size=4),
              legend.title = element_text(size=3),
              legend.key.size=unit(1, 'mm'),
              legend.margin=margin(1.5,1.5,1.5,1.5),
              legend.position=c(0.77,0.15)) +
	xlab('GTEx tissue') +
	ylab('Proportion of isoforms')

ggsave('../../figures/M2_GTEx_expression/M2_gtex_tx_unique_to_tissue_all_v_pc_proportion.pdf', width=58, height=58, units='mm')

########################################################################################
tiss_color <- setNames(hue_pal()(length(tissues_to_use)), levels(as.factor(tissues_to_use)))

# For each tissue, plot the number of isoforms expressed per gene (or protein-coding isoforms per gene) as a histogram
n_tx_per_gene <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_tx_per_gene_passing_thresholds_2023.tsv') %>%
	pivot_longer(!c('gene_id', 'gene_biotype', 'gene_name'), names_to = 'tissue', values_to = 'n_tx') %>%
	drop_na(n_tx)
n_tx_per_pc_gene <- read_tsv('../../tables/GTEx_expression/GTEx_number_of_protein_coding_tx_per_gene_passing_thresholds_2023.tsv') %>%
	pivot_longer(!c('gene_id', 'gene_biotype', 'gene_name'), names_to = 'tissue', values_to = 'n_tx') %>%
	drop_na(n_tx)

for (tx in tissues_to_use) {
	print(tx)
	print('gene bodies')
	tmp = n_tx_per_gene %>% 
		filter(tissue == tx)
	# print out the genes that have more than 30 isoforms expressed
	print(nrow(tmp %>% filter(n_tx >= 30)))
	t_bins = max(tmp$n_tx)
	ggplot(tmp, aes(x=n_tx, fill=tissue)) +
		geom_histogram(bins=t_bins) +
		theme(legend.position = 'none',
		      text=element_text(size=7),
		      axis.text=element_text(size=5)) + 
		ggtitle(tx) + 
		scale_fill_manual(values = tiss_color) +
		xlab('Number of isoforms') +
		ylab('Number of genes')
	ggsave(paste0("../../figures/M2_GTEx_expression/M2_n_tx_per_gene_", tx, ".pdf"), width = 88*2, height = 58*2, units='mm')

	# looking at the number of protein-coding isoforms per gene
	tmp_pc = n_tx_per_pc_gene %>% filter(tissue == tx)
	print('pc isoforms from gene bodies')
	print(tmp_pc %>% filter(n_tx >= 8))
	# print the number of genes that express more than 10 protein-coding isoforms
	print(nrow(tmp_pc %>% filter(n_tx >= 10)))
	t_bins_pc = max(tmp_pc$n_tx)
	ggplot(tmp_pc, aes(x=n_tx, fill=tissue)) +
		geom_histogram(bins=t_bins_pc) +
		theme(legend.position = 'none', 
		      text=element_text(size=7),
		      axis.text=element_text(size=5)) + 
		scale_fill_manual(values = tiss_color) +
		xlab('Number of protein-coding isoforms') +
		ylab('Number of genes')
	ggsave(paste0("../../figures/M2_GTEx_expression/M2_n_pc_tx_per_gene_", tx, ".pdf"), width = 58*2, height = 26*2, units='mm')
}
######################################################################################


##########TODO: STILL NEED TO SIZE THESE #########################

######################################################################################
# plot the number of isoforms per gene, across all tissue samples, based on unique counts,
# using 4 different thresholds: total unique counts > 1, 5, 10, and 20
unique_counts_n_tx <- read_tsv('../../tables/GTEx_expression/unique_counts_all_samples_across_thresholds_n_tx.tsv')

y_axis <- c(125, 20, 15, 7)
unique_counts_threshold <- c(1,5,10,20)
for (i in c(1:4)){
	plt <- ggplot(unique_counts_n_tx %>% filter(threshold == unique_counts_threshold[i]), aes(x=n_tx)) +
		geom_histogram() +
		xlab('Number of isoforms') +
		ylab('Number of genes')

	ggsave(paste0('../../figures/M2_GTEx_expression/n_tx_unique_counts_gt_', unique_counts_threshold[i], '.pdf'), plt, width = 88*2, height = 58*2, units='mm')

	zplt <- ggplot(unique_counts_n_tx %>% filter(threshold == unique_counts_threshold[i]), aes(x=n_tx)) +
		geom_histogram() +
		xlim(20,120) +
		ylim(0,y_axis[i]) +
		xlab('Number of isoforms') +
		ylab('Number of genes')
	ggsave(paste0('../../figures/M2_GTEx_expression/n_tx_unique_counts_gt_', unique_counts_threshold[i], '_zoomed.pdf'), zplt, width = 58*2, height = 26*2, units='mm')

}

##########################################################################################
############################ Transcript Biotype ##########################################

biotype <- read_tsv('../../tables/GTEx_expression/GTEx_transcript_biotype_expr_by_tissue_2023.tsv') 

ggplot(biotype, aes(x=reorder(transcript_biotype, -count), y=count, fill=transcript_biotype)) +
	geom_bar(stat='identity') +
	geom_text(mapping=aes(x=transcript_biotype, y=count, label = count), vjust = -0.2) +
	facet_wrap(vars(tissue)) +
	xlab('Transcript biotype') +
	ylab('Number of isoforms') +
	theme(legend.position='bottom',
	      axis.text.x=element_blank())
ggsave('../../figures/M2_GTEx_expression/GTEx_tissue_top_5_biotype.pdf', width=10,height=10)
