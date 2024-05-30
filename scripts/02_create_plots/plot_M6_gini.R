library(tidyverse)

housekeeping_ginic_overlap_lfc <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/housekeeping_ginic_overlap_with_lfc.tsv') %>%
	filter(n_not_sig == 36 | is.na(n_not_sig)) %>%
	mutate(`housekeeping detection method` = case_when(!is.na(giniC) & is.na(n_not_sig) ~ 'Gini',
			       !is.na(giniC) ~ 'Gini and lfc',
			       is.na(giniC) ~ 'lfc'))


ggplot(housekeeping_ginic_overlap_lfc, aes(x=isGiniGene, fill=`housekeeping detection method`)) +
	geom_bar() +
	geom_text(aes(label = after_stat(count)), stat='count', position = position_stack(vjust = .5))

ggsave('../../figures/M5_our_new_isoforms_in_GTEx/housekeeping_detection_method_breakdown.pdf')


