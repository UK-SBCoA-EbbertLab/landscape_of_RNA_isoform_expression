library(tidyverse)
set.seed(21)

h=10
w=11

# load in data
preferential_isoforms <- read_tsv("../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_1.01_preferential_filtered.tsv") %>%
	separate_wider_delim(TXNAME, delim = '_', names = c('TXNAME', 'GENEID')) %>%
	mutate(Type = if_else(grepl("Bambu",GENEID), 'From new gene', 'From known gene')) %>%
#	extract(GENEID, 'Type', "([a-zA-Z]+)", remove=FALSE) %>%
	group_by(Type)

ggplot(preferential_isoforms, aes(x=n_tis_pref, fill=Type)) +
	geom_bar() +
	theme(
	      legend.position = 'bottom', 
	      text=element_text(size=5), 
	      legend.key.size = unit(1.5, 'mm'), 
	      plot.margin = unit(c(0.5,0.5,0.5,0.5), 'mm'),
	      legend.margin=margin(t = 0, unit='mm')) +
	xlab('Number of preferentially expressed tissues') +
	ylab('Number of isoforms') 


ggsave('../../figures/M5_our_new_isoforms_in_GTEx/preferentially_expressed_isoforms.pdf', w=58 , h=58 , units='mm')

gini_gtex_genes <- read_tsv('gini_genes_gtex.tsv') %>%
	rename('gene_name'= Genes) %>%
	mutate(is_gini_gene = TRUE)

housekeeping_isoforms <- read_tsv("../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_thresh_1.01_lfc_2_housekeeping_filtered.tsv") %>%
	filter(n_not_sig == 36) %>%
	separate_wider_delim(TXNAME, delim = '_', names = c('TXNAME', 'GENEID')) %>%
        mutate(Type1 = if_else(grepl("Bambu",GENEID), 'From new gene', 'From known gene')) %>%
	left_join(gini_gtex_genes) %>%
	mutate(Type = case_when(
				is_gini_gene == TRUE ~ 'Matches a gini gene', 
				is.na(is_gini_gene) & grepl("ENSG", GENEID) ~ 'From known gene', 
				TRUE ~ 'From new gene')) %>%
	mutate(housekeeping = TRUE)

ggplot(housekeeping_isoforms, aes(x=housekeeping, fill=Type)) +
	geom_bar() +
	theme(
              legend.position = 'bottom',
              text=element_text(size=5),
              legend.key.size = unit(1.5, 'mm'),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), 'mm'),
              legend.margin=margin(t = 0, unit='mm')) +
	xlab('Housekeeping isoforms') +
	ylab('Number of isoforms')

ggsave('../../figures/M5_our_new_isoforms_in_GTEx/housekeeping_isoforms_stacked.pdf', w=58 , h=58 , units='mm')


ggplot(housekeeping_isoforms, aes(x=Type1, fill=Type)) +
	geom_bar() +
	theme(
              legend.position = 'bottom',
              text=element_text(size=5),
              legend.key.size = unit(1.5, 'mm'),
              plot.margin = unit(c(0.5,0.5,0.5,0.5), 'mm'),
              legend.margin=margin(t = 0, unit='mm')) +
	xlab('Housekeeping isoforms') +
	ylab('Number of isoforms')

ggsave('../../figures/M5_our_new_isoforms_in_GTEx/housekeeping_isoforms_spread.pdf', w=58 , h=58 , units='mm')






