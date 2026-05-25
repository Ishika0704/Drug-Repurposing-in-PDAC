# =========================================
# GSE Validation + Intersection with TCGA
# =========================================

suppressPackageStartupMessages({
  library(limma)
})

# ---------- FILE PATHS (EDIT THESE)
gse_path      <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/geo_clean_expression_matrix.csv"
tcga_de_path  <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/TCGA_GTEx_DE_results.csv"
output_gse_de <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/GSE_DE_results.csv"
output_intersection <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/HighConfidence_Signature.csv"

# ---------- LOAD GSE DATA
gse <- read.csv(gse_path, row.names = 1, check.names = FALSE)

gse[] <- lapply(gse, as.numeric)

cat("GSE dimensions:", dim(gse), "\n")

# -------------------------------------------------------
# ⚠ EDIT THIS LINE BASED ON YOUR SAMPLE NAMING
# Example assumes column names contain "tumor" or "normal"
# -------------------------------------------------------

group <- ifelse(grepl("tumor", colnames(gse), ignore.case = TRUE),
                "Tumor", "Normal")

group <- factor(group, levels = c("Normal", "Tumor"))

table(group)

# ---------- Log transform if needed
max_val <- max(as.matrix(gse), na.rm = TRUE)

if (max_val > 50) {
  gse <- log2(gse + 1)
  cat("Log2 transform applied to GSE.\n")
}

# ---------- Normalize
gse <- normalizeBetweenArrays(gse)

# ---------- Filter low-expression genes
keep <- rowMeans(gse) > 1
gse <- gse[keep, ]

cat("Genes after filtering:", nrow(gse), "\n")

# ---------- Differential Expression
design <- model.matrix(~ group)

fit <- lmFit(gse, design)
fit <- eBayes(fit)

gse_results <- topTable(
  fit,
  coef = "groupTumor",
  number = Inf,
  adjust.method = "BH"
)

write.csv(gse_results, output_gse_de)

cat("GSE DE complete.\n")
cat("Significant genes (FDR < 0.05):",
    sum(gse_results$adj.P.Val < 0.05), "\n")

# =========================================
# INTERSECTION WITH TCGA RESULTS
# =========================================

tcga_results <- read.csv(tcga_de_path, row.names = 1)

common_genes <- intersect(rownames(tcga_results),
                          rownames(gse_results))

cat("Common genes:", length(common_genes), "\n")

# Subset
tcga_sub <- tcga_results[common_genes, ]
gse_sub  <- gse_results[common_genes, ]

# Direction agreement
direction_match <- sign(tcga_sub$logFC) ==
  sign(gse_sub$logFC)

cat("Direction agreement proportion:",
    mean(direction_match), "\n")

# High-confidence signature
high_conf <- common_genes[
  tcga_sub$adj.P.Val < 0.05 &
    gse_sub$adj.P.Val < 0.05 &
    direction_match
]

cat("High-confidence replicated genes:",
    length(high_conf), "\n")

final_signature <- data.frame(
  Gene = high_conf,
  TCGA_logFC = tcga_sub[high_conf, "logFC"],
  GSE_logFC  = gse_sub[high_conf, "logFC"]
)

write.csv(final_signature, output_intersection, row.names = FALSE)

cat("Validation complete.\n")