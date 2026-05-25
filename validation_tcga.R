install.packages(c("limma","pheatmap","dplyr","ggplot2"))
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)
geo  <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/GSE205154_pancreas_primary_expression.csv", stringsAsFactors = FALSE)
gtex <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/GTEx_pancreas_final_clean_noZero.csv", stringsAsFactors = FALSE)

rownames(geo)  <- geo$gene
rownames(gtex) <- gtex$GeneSymbol

geo  <- geo[, -1]
gtex <- gtex[, -1]
common_genes <- intersect(rownames(geo), rownames(gtex))

geo  <- geo[common_genes, ]
gtex <- gtex[common_genes, ]
combined <- cbind(geo, gtex)

group <- factor(c(
  rep("Tumour", ncol(geo)),
  rep("Normal", ncol(gtex))
))

design <- model.matrix(~ group)
colnames(design) <- c("Intercept","Tumour_vs_Normal")
fit <- lmFit(combined, design)
fit <- eBayes(fit)

res_geo <- topTable(fit,
                    coef = "Tumour_vs_Normal",
                    adjust.method = "BH",
                    number = Inf)
sig_genes_strict <- subset(res_geo,
                           adj.P.Val < 0.01 &
                             abs(logFC) > 2)

nrow(sig_genes_strict)

nrow(sig_genes_strict)
tcga <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/TCGA_primary_tumor_only.csv", stringsAsFactors = FALSE)
rownames(tcga) <- tcga$sample
tcga <- tcga[, -1]
validation_genes <- intersect(rownames(tcga), rownames(sig_genes_strict))

tcga_sub <- tcga[validation_genes, ]
sig_genes_strict$Direction <- sign(sig_genes_strict$logFC)

tcga_mean <- rowMeans(tcga_sub)

validation_df <- data.frame(
  Gene = validation_genes,
  GEO_logFC = sig_genes[validation_genes,"logFC"],
  TCGA_mean = tcga_mean
)

head(validation_df)
cor(validation_df$GEO_logFC,
    validation_df$TCGA_mean,
    use="complete.obs")
