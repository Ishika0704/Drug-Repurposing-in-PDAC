# =====================================================
# Z-SCORE CALCULATION USING GTEx AS NORMAL REFERENCE
# TCGA + GSE = Tumour dataset
# =====================================================

# -----------------------------
# 1. Load libraries
# -----------------------------
library(dplyr)

cat("Step 1: Libraries loaded\n\n")


# -----------------------------
# 2. Load datasets
# -----------------------------
tcga <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/TCGA_primary_tumor_only.csv",
                 row.names = 1, check.names = FALSE)

gse <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/GSE205154_pancreas_primary_expression.csv",
                row.names = 1, check.names = FALSE)

gtex <- read.csv("C:/Users/Ishika/OneDrive/Desktop/Capstone Project/GTEx_pancreas_final_clean_noZero.csv",
                 row.names = 1, check.names = FALSE)

cat("Step 2: Data loaded\n")
cat("TCGA dimensions:", dim(tcga), "\n")
cat("GSE dimensions:", dim(gse), "\n")
cat("GTEx dimensions:", dim(gtex), "\n\n")


# -----------------------------
# 3. Convert to matrices
# -----------------------------
tcga <- as.matrix(tcga)
gse  <- as.matrix(gse)
gtex <- as.matrix(gtex)

cat("Step 3: Converted to matrices\n\n")


# -----------------------------
# 4. Check gene overlap
# -----------------------------
tcga_genes <- rownames(tcga)
gse_genes  <- rownames(gse)
gtex_genes <- rownames(gtex)

common_genes <- Reduce(intersect, list(tcga_genes, gse_genes, gtex_genes))

cat("Step 4: Gene overlap diagnostics\n")
cat("Genes in TCGA:", length(tcga_genes), "\n")
cat("Genes in GSE:", length(gse_genes), "\n")
cat("Genes in GTEx:", length(gtex_genes), "\n")
cat("Common genes:", length(common_genes), "\n\n")


# -----------------------------
# 5. Subset datasets
# -----------------------------
tcga <- tcga[common_genes, ]
gse  <- gse[common_genes, ]
gtex <- gtex[common_genes, ]

cat("Step 5: Subsetting complete\n")
cat("TCGA:", dim(tcga), "\n")
cat("GSE:", dim(gse), "\n")
cat("GTEx:", dim(gtex), "\n\n")


# -----------------------------
# 6. Combine tumour datasets
# -----------------------------
tumor <- cbind(tcga, gse)

cat("Step 6: Tumour datasets merged\n")
cat("Tumour matrix dimensions:", dim(tumor), "\n")
cat("Total tumour samples:", ncol(tumor), "\n\n")


# -----------------------------
# 7. Calculate GTEx statistics
# -----------------------------
normal_mean <- rowMeans(gtex)
normal_sd   <- apply(gtex, 1, sd)

cat("Step 7: Normal statistics calculated\n")
cat("SD summary:\n")
print(summary(normal_sd))
cat("\n")


# -----------------------------
# 8. Filter low variance genes
# -----------------------------
low_var_genes <- normal_sd < 0.1

cat("Step 8: Variance filtering\n")
cat("Genes with SD < 0.1:", sum(low_var_genes), "\n")

tcga  <- tcga[!low_var_genes, ]
gse   <- gse[!low_var_genes, ]
gtex  <- gtex[!low_var_genes, ]
tumor <- tumor[!low_var_genes, ]

normal_mean <- normal_mean[!low_var_genes]
normal_sd   <- normal_sd[!low_var_genes]

cat("Genes remaining after filtering:", length(normal_mean), "\n\n")


# -----------------------------
# 9. Z-score calculation
# -----------------------------
z_scores <- sweep(tumor, 1, normal_mean, "-")
z_scores <- sweep(z_scores, 1, normal_sd, "/")

cat("Step 9: Z-score calculation complete\n")
cat("Z-score matrix dimensions:", dim(z_scores), "\n\n")


# -----------------------------
# 10. Z-score diagnostics
# -----------------------------
cat("Step 10: Z-score summary statistics\n")
print(summary(as.vector(z_scores)))

cat("\nExample Z-score values:\n")
print(z_scores[1:5,1:5])
cat("\n")


# -----------------------------
# 11. Save output file
# -----------------------------
z_scores_df <- as.data.frame(z_scores)

write.csv(z_scores_df,
          "PDAC_Tumor_Zscores_GTEx_reference_filtered.csv")

cat("Step 11: File saved successfully\n")
cat("Output file: PDAC_Tumor_Zscores_GTEx_reference_filtered.csv\n")
quantile(z_scores, probs = c(0.9,0.95,0.99))
getwd()
# -----------------------------
# 12. Calculate combined gene Z-score
# -----------------------------

gene_mean_z <- rowMeans(z_scores)

combined_z <- data.frame(
  Gene = rownames(z_scores),
  Mean_Zscore = gene_mean_z
)

cat("Combined gene Z-scores created\n")
print(head(combined_z))


# -----------------------------
# 13. Save combined gene scores
# -----------------------------

write.csv(combined_z,
          "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/PDAC_Gene_Mean_Zscores.csv",
          row.names = FALSE)

cat("File saved: PDAC_Gene_Mean_Zscores.csv\n")
head(sort(gene_mean_z, decreasing=TRUE),20)
head(sort(gene_mean_z),20)
