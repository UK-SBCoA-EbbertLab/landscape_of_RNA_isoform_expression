library(tidyverse)

expressed_isoforms <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv')

brain_disease_genes <- read_tsv('../../references/brain_disease_genes_with_disease.tsv')

overlap <- left_join(brain_disease_genes, expressed_isoforms)

write_tsv(overlap, '../../tables/expression_of_brain_disease_genes.tsv')

