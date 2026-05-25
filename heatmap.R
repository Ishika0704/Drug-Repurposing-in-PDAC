install.packages(c("pheatmap","dplyr"))
library(pheatmap)
library(dplyr)
final_deg <- read.csv("C:/Users/Ishika/Downloads/FINAL DATASET - Final curated dataset (1).csv", stringsAsFactors = FALSE)
colnames(final_deg) <- make.names(colnames(final_deg))

# Remove duplicated header row
final_deg <- final_deg[-1, ]

gene_list <- final_deg$GENE.ID
tcga <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/TCGA_primary_tumor_only.csv", stringsAsFactors = FALSE)
geo  <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/GSE205154_pancreas_primary_expression.csv", stringsAsFactors = FALSE)
gtex <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/GTEx_pancreas_final_clean_noZero.csv", stringsAsFactors = FALSE)

#diagnostics
str(tcga[1:5, 1:5])
head(tcga[, 1:5])
colnames(tcga)[1:5]

str(geo[1:5, 1:5])
head(geo[, 1:5])
colnames(geo)[1:5]

str(gtex[1:5, 1:5])
head(gtex[, 1:5])
colnames(gtex)[1:5]
# Set gene names as rownames
rownames(tcga) <- tcga$sample
tcga <- tcga[, -1]

rownames(geo) <- geo$gene
geo <- geo[, -1]

rownames(gtex) <- gtex$GeneSymbol
gtex <- gtex[, -1]
common_genes <- Reduce(intersect, list(
  gene_list,
  rownames(tcga),
  rownames(geo),
  rownames(gtex)
))
length(common_genes)
tcga_sub <- tcga[common_genes, ]
geo_sub  <- geo[common_genes, ]
gtex_sub <- gtex[common_genes, ]
combined_matrix <- cbind(tcga_sub, geo_sub, gtex_sub)

combined_matrix <- apply(combined_matrix, 2, as.numeric)
rownames(combined_matrix) <- common_genes
annotation <- data.frame(
  Tissue = c(
    rep("TCGA_Tumour", ncol(tcga_sub)),
    rep("GEO_Tumour", ncol(geo_sub)),
    rep("GTEx_Normal", ncol(gtex_sub))
  )
)

rownames(annotation) <- colnames(combined_matrix)
heat_matrix_scaled <- t(scale(t(combined_matrix)))
pheatmap(
  heat_matrix_scaled,
  annotation_col = annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  fontsize = 8,
  color = colorRampPalette(c("blue","white","red"))(100),
  main = "Final Gene Signature Across TCGA, GEO and GTEx"
)
pca <- prcomp(t(heat_matrix_scaled), scale = FALSE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Tissue = annotation$Tissue
)

library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Final Gene Signature")
