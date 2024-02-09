## File and subdirectory description

`./cemetery/` - Directory containing old scripts no longer being used

`./01_make_tables_from_annotations.ipynb` - Ipython notebook used to analyse the difference in annotations across the years

`./02_make_tables_GTEx_all_isoforms.ipynb` - Ipython notebook used to analyse isoform expression across nine GTEx tissues

`./02a_create_unique_counts_tables.R` - R script to create additional table(s) with unique counts

`./02b_calculate_pvalue.R` - R script to calculate the pvalue of a comparison used in the manuscript

`./03_heatmap_n_tx_per_gene.R` - R script used to create the heatmap showcasing the number of isoforms expressed for a subset of the isoforms

`./03x_interrogate_specific_genes.R` - R script used to look at specific genes from the heatmap

`./04_make_tables_GTEx_new_since_2019.ipynb` - Ipython notebook to assess isoform expression across tissues of isoforms that are present in the Ensembl v109 annotation compared to the v96 annotation 

`./04a_create_combined_relative_expr_table.R` - R script used to create a relative expression table for the isoforms that are new since Ensembl v96 annotation

`./05_GTEx_expression_of_new_isoforms.R` - R script to verify and quantify isoforms discovered in Aguzzoli-Heberle et al. expressed across the nine GTEx tissues

`./05a_all_new_isoforms_deseq2.R` - R script to run deseq2 on all isoforms and then perform pairwise comparisons for all new isoforms

`./05b_determine_housekeeping_and_preferential.R` - R script to assess whether a new isoform is preferenentially expressed in a tissue or few, or if it could be considered a housekeeping isoform

`./run_0*.sh` - Bash scripts to run the various R scripts
