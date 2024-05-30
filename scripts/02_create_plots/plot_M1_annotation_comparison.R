library(tidyverse)
library(ggpubr)

# load tables created from jupyter notebook
n_genes_per_year <- read_tsv('../../tables/annotation_comparison/number_of_genes_per_year.tsv')
n_pc_genes_per_year <- read_tsv('../../tables/annotation_comparison/number_of_protein_coding_genes_per_year.tsv')
n_transcript_per_year <- read_tsv('../../tables/annotation_comparison/number_of_transcripts_per_year.tsv')
n_pc_transcript_per_year <- read_tsv('../../tables/annotation_comparison/number_of_protein_coding_transcripts_per_year.tsv')
tx_by_biotype <- read_tsv('../../tables/annotation_comparison/tx_and_gene_biotype_for_2023.tsv')

###### MODULE 1 -- ANNOTATION COMPARISONS ################

# plot number of genes per year
ggplot(n_genes_per_year, aes(x=Year,y=n_genes)) +
	theme(text=element_text(size = 5)) +
	geom_bar(stat='identity', fill="#00BFC4") +
	xlab('Year of Ensembl annotation') +
	ylab('Number of annotated genes') +
	coord_cartesian(ylim = c(56000, 64000)) +
	geom_text(aes(label = n_genes), vjust = -0.2, size = 1) +
	theme(legend.position = "none") + 
        scale_x_continuous(breaks=c(2014,2015,2016,2017,2018,2019,2020,2021,2022,2023))
ggsave('../../figures/M1_annotation_comparison/M1_genes_per_year_annotation.pdf', width = 58, height = 58, units='mm')

# plot number of protein-coding genes per year
ggplot(n_pc_genes_per_year, aes(x=Year,y=n_pc_genes)) +
	theme(text=element_text(size = 5)) +
        geom_bar(stat='identity', fill="#00BFC4") +
        xlab('Year of Ensembl annotation') +
        ylab('Number of annotated protein-coding genes') +
        coord_cartesian(ylim = c(19700, 20200)) +
        geom_text(aes(label = n_pc_genes), vjust = -0.2, size = 1) +
	theme(legend.position = "none") + 
        scale_x_continuous(breaks=c(2014,2015,2016,2017,2018,2019,2020,2021,2022,2023))
ggsave('../../figures/M1_annotation_comparison/M1_pc_genes_per_year_annotation.pdf', width = 58, height = 58, units='mm')

# plot number of transcripts per year
ggplot(n_transcript_per_year, aes(x=Year,y=`Number of transcripts`)) +
	geom_bar(stat='identity', fill="#F8766D") +
	theme(text=element_text(size = 5)) +
	xlab('Year of Ensembl annotation') +
	ylab('Number of annotated isoforms') +
	coord_cartesian(ylim = c(190000, 260000)) +
	geom_text(aes(label = `Number of transcripts`), vjust = -0.2, size=1) +
	theme(legend.position = "none") + 
        scale_x_continuous(breaks=c(2014,2015,2016,2017,2018,2019,2020,2021,2022,2023))
ggsave('../../figures/M1_annotation_comparison/M1_isoforms_per_year_annotation.pdf', width = 58, height = 58, units='mm')

# plot number of protein-coding transcripts per year
ggplot(n_pc_transcript_per_year, aes(x=Year,y=`Number of transcripts`, fill="#F8766D")) +
        geom_bar(stat='identity') +
	theme(text=element_text(size = 5)) +
        xlab('Year of Ensembl annotation') +
        ylab('Number of annotated protein-coding isoforms') +
        coord_cartesian(ylim = c(77500, 91000)) +
        geom_text(aes(label = `Number of transcripts`), vjust = -0.2, size=1) +
	theme(legend.position = "none") + 
        scale_x_continuous(breaks=c(2014,2015,2016,2017,2018,2019,2020,2021,2022,2023))
ggsave('../../figures/M1_annotation_comparison/M1_pc_isoforms_per_year_annotation.pdf', width = 58, height = 58, units='mm')


# plot the number of transcripts by their biotype
ggplot(tx_by_biotype, aes(x=forcats::fct_rev(forcats::fct_infreq(transcript_biotype)), fill=forcats::fct_infreq(transcript_biotype))) +
	theme(text=element_text(size = 5)) +
	geom_bar() +
	xlab('Transcript Biotype') +
	ylab('Number of transcripts') +
	#geom_text(aes(label = after_stat(count)), stat = "count", hjust=-.05, angle=90) +
	geom_text(aes(label = after_stat(count)), stat = "count", hjust=-.05, size=1) +
	ylim(0, 98000) +
	theme(legend.position = "none",
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	coord_flip()
ggsave('../../figures/M1_annotation_comparison/M1_tx_biotype_2023.pdf', width = 58, height=58*2+2, units='mm')
 

# plot the number of isoforms per gene
n_tx_per_gene <- read_tsv('../../tables/annotation_comparison/n_tx_per_gene_2023.tsv')
n_tx_per_pc_gene <- read_tsv('../../tables/annotation_comparison/n_tx_per_pc_gene_2023.tsv')
n_tx_per_npc_gene <- read_tsv('../../tables/annotation_comparison/n_tx_per_npc_gene_2023.tsv')


#### all genes (number of isoforms per gene) #########################################################
ggplot(n_tx_per_gene, aes(x=n_tx)) + 
	theme(text=element_text(size = 5)) +
	geom_histogram(bins=150) +
	xlab('Number of isoforms per gene') +
	ylab('Number of genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_gene_annotation.pdf", width = 88, height = 58, units='mm')

ggplot(n_tx_per_gene, aes(x=n_tx)) + 
	theme(text=element_text(size = 4)) +
	geom_histogram(bins=240, fill = '#584B9F') + 
	xlim(60, 300) + 
	ylim(0, 8) +
	xlab('Number of isoforms per gene') +
	ylab('Number of genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_gene_annotation_60-300.pdf", width = 30, height = 20, units='mm')

ggplot(n_tx_per_gene, aes(x=n_tx)) + 
	theme(text=element_text(size = 4)) +
	geom_histogram(bins=30, fill = '#F9C25C') + 
	xlim(30, 60) + 
	ylim(0, 75) +
	xlab('Number of isoforms per gene') +
	ylab('Number of genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_gene_annotation_30-60.pdf", width = 30, height = 20, units='mm')

ggplot(n_tx_per_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=20, fill = '#A71B4B') +
        xlim(10, 30) +
        ylim(0, 1700) +
	xlab('Number of isoforms per gene') +
	ylab('Number of genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_gene_annotation_10-30.pdf", width = 30, height = 20, units='mm')
########################################################################


##### just pc genes (isoforms per gene) ####################################################
ggplot(n_tx_per_pc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 5)) +
        geom_histogram(bins=50) +
	xlab('Number of isoforms per protein-coding gene') +
	ylab('Number of protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_pc_gene_annotation.pdf", width = 88, height = 58, units='mm')

ggplot(n_tx_per_pc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=140, fill = '#584B9F') +
        xlim(60, 200) +
        ylim(0, 8) +
	xlab('Number of isoforms per protein-coding gene') +
	ylab('Number of protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_pc_gene_annotation_60-200.pdf", width = 30, height = 20, units='mm')

ggplot(n_tx_per_pc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=30, fill = '#F9C25C') +
        xlim(30, 60) +
        ylim(0, 50) +
	xlab('Number of isoforms per protein-coding gene') +
	ylab('Number of protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_pc_gene_annotation_30-60.pdf", width = 30, height = 20, units='mm')

ggplot(n_tx_per_pc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=20, fill = '#A71B4B') +
        xlim(10, 30) +
        ylim(0, 1500) +
	xlab('Number of isoforms per protein-coding gene') +
	ylab('Number of protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_pc_gene_annotation_10-30.pdf", width = 30, height = 20, units='mm')
#####################################################################

########## non protein coding genes (isoforms per gene)#################################
ggplot(n_tx_per_npc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 5)) +
        geom_histogram(bins=150) +
	xlab('Number of isoforms per non protein-coding gene') +
	ylab('Number of non protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_npc_gene_annotation.pdf", width = 88, height = 58, units='mm')

ggplot(n_tx_per_npc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=240, fill = '#584B9F') +
        xlim(60, 300) +
        ylim(0, 3) +
	xlab('Number of isoforms per non protein-coding gene') +
	ylab('Number of non protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_npc_gene_annotation_60-300.pdf", width = 30, height = 20, units='mm')

ggplot(n_tx_per_npc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=30, fill = '#F9C25C') +
        xlim(30, 60) +
        ylim(0, 15) +
	xlab('Number of isoforms non protein-coding gene') +
	ylab('Number of non protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_npc_gene_annotation_30-60.pdf", width = 30, height = 20, units='mm')

ggplot(n_tx_per_npc_gene, aes(x=n_tx)) +
	theme(text=element_text(size = 4)) +
        geom_histogram(bins=20, fill = '#A71B4B') +
        xlim(10, 30) +
        ylim(0, 250) +
	xlab('Number of isoforms per non protein-coding gene') +
	ylab('Number of non protein-coding genes')

ggsave(filename="../../figures/M1_annotation_comparison/number_tx_per_npc_gene_annotation_10-30.pdf", width = 30, height = 20, units='mm')
#######################################################################
