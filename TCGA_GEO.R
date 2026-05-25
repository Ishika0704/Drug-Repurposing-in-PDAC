library(limma)
library(dplyr)
library(ggplot2)
tcga <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/TCGA_noZero_genes.csv", stringsAsFactors = FALSE)

rownames(tcga) <- tcga[,1]
tcga <- tcga[,-1]
colnames(tcga)[1:20]
group <- ifelse(substr(colnames(tcga), 14, 15) == "11",
                "Normal",
                "Tumour")

table(group)
group <- factor(group, levels = c("Normal","Tumour"))

design <- model.matrix(~ group)
colnames(design) <- c("Intercept","Tumour_vs_Normal")
fit <- lmFit(tcga, design)
fit <- eBayes(fit)

res_tcga <- topTable(fit,
                     coef = "Tumour_vs_Normal",
                     adjust.method = "BH",
                     number = Inf)
sig_tcga <- subset(res_tcga,
                   adj.P.Val < 0.05 &
                     abs(logFC) > 1)

nrow(sig_tcga)
geo  <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/GSE205154_log2TPM_noZero_TCGA_like.csv", stringsAsFactors = FALSE)
rownames(geo) <- geo[,1]
geo  <- geo[,-1]
validation_genes <- intersect(rownames(sig_tcga),
                              rownames(geo))

length(validation_genes)
geo_sub <- geo[validation_genes, ]

geo_mean <- rowMeans(geo_sub)

validation_df <- data.frame(
  Gene = validation_genes,
  TCGA_logFC = sig_tcga[validation_genes,"logFC"],
  GEO_mean = geo_mean
)

cor(validation_df$TCGA_logFC,
    validation_df$GEO_mean,
    use="complete.obs")
validation_df$Direction_TCGA <- sign(validation_df$TCGA_logFC)

# Standardize GEO means to center them
geo_centered <- geo_mean - mean(geo_mean)
validation_df$Direction_GEO <- sign(geo_centered)

robust_genes <- subset(validation_df,
                       Direction_TCGA == Direction_GEO)

nrow(robust_genes)
sig_tcga_strict <- subset(res_tcga,
                          adj.P.Val < 0.01 &
                            abs(logFC) > 1.5)

length(intersect(rownames(sig_tcga_strict),
                 rownames(robust_genes)))
final_genes <- intersect(rownames(sig_tcga_strict),
                         rownames(robust_genes))

geo_sub_final <- geo[final_genes, ]
geo_mean_final <- rowMeans(geo_sub_final)

validation_final <- data.frame(
  Gene = final_genes,
  TCGA_logFC = sig_tcga_strict[final_genes,"logFC"],
  GEO_mean = geo_mean_final
)

cor(validation_final$TCGA_logFC,
    validation_final$GEO_mean,
    use="complete.obs")
genes_124 <- robust_genes$Gene

write.csv(
  data.frame(Gene = genes_124),
  "TCGA_GEO_robust_124_genes.csv",
  row.names = FALSE
)
write.csv(
  robust_genes,
  "TCGA_GEO_robust_124_genes_full_info.csv",
  row.names = FALSE
)
final_signature_df <- data.frame(
  Gene = final_genes,
  TCGA_logFC = sig_tcga_strict[final_genes,"logFC"]
)

write.csv(
  final_signature_df,
  "TCGA_GEO_strict_45_gene_signature_full_info.csv",
  row.names = FALSE
)
