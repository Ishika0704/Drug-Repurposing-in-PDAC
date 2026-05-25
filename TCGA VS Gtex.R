# =========================================
# TCGA vs GTEx Differential Expression
# Fully optimized & statistically correct
# =========================================

suppressPackageStartupMessages({
  library(limma)
})

# ---------- File paths
tcga_path  <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/TCGA_primarytumor_only.csv"
gtex_path  <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/GTEx_noZero_genes.csv"
output_path <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/TCGA_GTEx_DE_results.csv"

# ---------- Load
tcga <- read.csv(tcga_path, row.names = 1, check.names = FALSE)
gtex <- read.csv(gtex_path, row.names = 1, check.names = FALSE)

tcga[] <- lapply(tcga, as.numeric)
gtex[] <- lapply(gtex, as.numeric)

# ---------- Gene intersection
genes <- intersect(rownames(tcga), rownames(gtex))
tcga <- tcga[genes, ]
gtex <- gtex[genes, ]

# ---------- Log transform if needed
max_val <- max(
  max(as.matrix(tcga), na.rm = TRUE),
  max(as.matrix(gtex), na.rm = TRUE)
)

if (max_val > 50) {
  tcga <- log2(tcga + 1)
  gtex <- log2(gtex + 1)
}

# ---------- Merge
expr <- cbind(tcga, gtex)

group <- factor(
  c(rep("Tumor", ncol(tcga)),
    rep("Normal", ncol(gtex))),
  levels = c("Normal", "Tumor")
)

# ---------- Normalize between samples
expr <- normalizeBetweenArrays(expr)

# ---------- Filter low-expression genes (CRITICAL)
keep <- rowMeans(expr) > 1
expr <- expr[keep, ]

# ---------- limma DE
design <- model.matrix(~ group)

fit <- lmFit(expr, design)
fit <- eBayes(fit)

results <- topTable(
  fit,
  coef = "groupTumor",
  number = Inf,
  adjust.method = "BH"
)

# ---------- Means
results$TumorMean  <- rowMeans(expr[, group == "Tumor"])
results$NormalMean <- rowMeans(expr[, group == "Normal"])

# ---------- Save
write.csv(results, output_path)

cat("Genes tested:", nrow(results), "\n")
cat("Significant (FDR < 0.05):", sum(results$adj.P.Val < 0.05), "\n")
results$adj.P.Val
min(results$P.Value[results$P.Value > 0])
results$log10P <- -log10(results$P.Value + 1e-300)
summary(results$adj.P.Val)
