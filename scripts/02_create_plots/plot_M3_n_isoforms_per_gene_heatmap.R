library(tidyverse)
library(pheatmap)
library(scales)

## Functions Mark Ebbert wrote for the heatmap
adj_contrast<-function(x, contrast){
  
  if(is.na(contrast)){
    return(x)
  }
  
  x[x > contrast] <- contrast
  x[x < -1 * contrast] <- -1 * contrast
  return(x)
}


standardize<-function(x){
	annAll<-dimnames(x)
	x<-scale(x)
	dimnames(x)<-annAll
	return(x)
}

medianCtr<-function(x){
	annAll<-dimnames(x)
	medians <- apply(x,1,median,na.rm=T)
	x <- t(scale(t(x),center=medians,scale=F))
	dimnames(x) <- annAll
	return(x)
}

#####LOAD DATA############################################

#Load n isoforms based on CPM threshold
n_tx_per_gene <- read_tsv('../../tables/GTEx_expression/GTEx_any_n_tx_by_tissue_gt_5.tsv') %>%
	arrange(gene_id)

############################################
######### HEATMAPS OF RAW VALUES ###########
############################################
# look at n isoforms based on CPM threshold
# set the expression matrix
expression.matrix <- n_tx_per_gene
x <- as.data.frame(expression.matrix[,-c(1,2,3,13,14)]) # remove gene_id, gene_biotype, min iso, max iso, and the hgnc columns

# Sort columns for x alphabetically
x <- as_tibble(x) %>%
  dplyr::select(order(colnames(x)))
x<-as.data.frame(x)

# set the rownames
rownames(x) <- paste(expression.matrix$gene_id, expression.matrix$gene_name, expression.matrix$gene_biotype, sep = ",")

# Identify rows with non-zero standard deviation in 'x'
rows_to_keep <- which(apply(x, 1, sd) != 0)
x <- x[rows_to_keep, ]

# CODE FOR THE HEATMAP -- WE WILL CREATE 2, ONE WITH GENE NAMES AND ONE WITHOUT GENE NAMES
print(x)
# set the annotations for the heatmap (gene biotype and tissues)
annos = data.frame(gene_biotype = sapply(str_split(rownames(x),','), function(y) y[3]))
rownames(annos) = rownames(x)

annos2 = data.frame(Tissue = colnames(x))
rownames(annos2) = colnames(x)
tissue_colors <- setNames(hue_pal()(9), annos2$Tissue)

w=14
h=56
pdf('../../figures/M3_GTEx_n_isoforms_heatmap/M3_n_tx_by_gene_CPM_heatmap_raw_values_7_contrast.pdf', width = w, height = h, onefile=TRUE)

# create the heatmap
tmp.heat <- pheatmap(adj_contrast(x, 7),
         clustering_distance_rows = as.dist(1-cor(t(x), method = "pearson")),
         clustering_distance_cols = dist(1-cor(x, method = "spearman")),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
	 color = colorRampPalette(c('#66FFFF', '#FFFFFF', '#47005c'))(30),
	 #color = colorRampPalette(c('#00FFF6', '#FFFFFF', '#4B0082'))(30),
         clustering_method = "complete",
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_row = annos,
	 annotation_col = annos2,
	 annotation_colors = list(Tissue = tissue_colors)
         )
# heatmap with names
tmp2.heat <- pheatmap(adj_contrast(x, 7),
         clustering_distance_rows = as.dist(1-cor(t(x), method = "pearson")),
         clustering_distance_cols = dist(1-cor(x, method = "spearman")),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
	 color = colorRampPalette(c('#66FFFF', '#FFFFFF', '#47005c'))(30),
	 #color = colorRampPalette(c('#00FFF6', '#FFFFFF', '#4B0082'))(30),
         clustering_method = "complete",
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_row = annos,
	 annotation_col = annos2,
	 annotation_colors = list(Tissue = tissue_colors)
         )
dev.off()

############## HEATMAPS USING STANDARDIZED AND MEDIAN CENTERED VALUES ###################
# standardized and median centered heatmaps
x <- standardize(medianCtr(x))

pdf('../../figures/M3_GTEx_n_isoforms_heatmap/M3_n_tx_by_gene_CPM_heatmap_standardized_and_centered_100_contrast.pdf', width = w, height = h, onefile=TRUE)

# create standardized heatmap
tmp.heat <- pheatmap(adj_contrast(x,1),
         clustering_distance_rows = as.dist(1-cor(t(x), method = "pearson")),
         clustering_distance_cols = dist(1-cor(x, method = "spearman")),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_method = "complete",
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_row = annos,
         annotation_col = annos2,
	 annotation_colors = list(Tissue = tissue_colors)
         )

# heatmap without gene names
tmp2.heat <- pheatmap(adj_contrast(x, 1),
         clustering_distance_rows = as.dist(1-cor(t(x), method = "pearson")),
         clustering_distance_cols = dist(1-cor(x, method = "spearman")),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_method = "complete",
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_row = annos,
         annotation_col = annos2,
	 annotation_colors = list(Tissue = tissue_colors)
         )

dev.off()

