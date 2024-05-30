library(tidyverse)

# load in relative abundance data and rename the relative abundance column
all_rel_abund <- read_tsv('../../tables/GTEx_expression_since_2019/gtex_cpm_gt_1_relative_abundance_new_since_2019.tsv') %>%
	replace(is.na(.), 0)

#grab all the categories
all_categories <- unique(all_rel_abund %>% pull(Category))

output_values=tibble(
		     Category = character(),
		     total_genes = numeric(),
		     gt_25 = numeric(),
		     gt_50 = numeric(),
		     gt_75 = numeric(),
		     lt_or_eq_25 = numeric()
		     )

# for each category, calculate the total number of genes that have isoforms in the category and determine the number of genes for the thresholds
for (cat in all_categories) {
	tmp <- all_rel_abund %>%
		filter(Category == cat)
	total <- nrow(tmp)
	gt_25 <- nrow(tmp %>%
		      pivot_longer(!c(gene_id, gene_name, Category), names_to='tissue', values_to='rel_abund') %>% 
		      filter(rel_abund > 25) %>%
		      select(gene_id, gene_name, Category) %>%
		      distinct())
	gt_50 <- nrow(tmp %>%
		      pivot_longer(!c(gene_id, gene_name, Category), names_to='tissue', values_to='rel_abund') %>% 
		      filter(rel_abund > 50) %>%
		      select(gene_id, gene_name, Category) %>%
		      distinct())
	gt_75 <- nrow(tmp %>%
		      pivot_longer(!c(gene_id, gene_name, Category), names_to='tissue', values_to='rel_abund') %>% 
		      filter(rel_abund > 75) %>%
		      select(gene_id, gene_name, Category) %>%
		      distinct())
	lt_or_eq_25 <- nrow(tmp %>%
			    pivot_longer(!c(gene_id, gene_name, Category), names_to='tissue', values_to='rel_abund') %>%
			    filter(rel_abund <= 25) %>%
			    group_by(gene_id, gene_name, Category) %>%
			    summarise(n_tiss = n()) %>%
			    filter(n_tiss == 9) %>%
			    select(gene_id, gene_name, Category) %>%
			    distinct())

	output_values <- output_values %>% add_row(Category = cat,
						   total_genes = total, 
						   gt_25 = gt_25,
						   gt_50 = gt_50, 
						   gt_75 = gt_75,
						   lt_or_eq_25 = lt_or_eq_25)

}

write_tsv(output_values, '../../tables/GTEx_expression_since_2019/gtex_n_genes_in_relative_abundance_categories.tsv')


