# ===============================
# LIBRARIES
# ===============================
library(ggplot2)
library(ggrepel)
library(dplyr)

# ===============================
# DATA LOADING & CLEANING
# ===============================
deg <- read.csv("C:/Users/Ishika/Downloads/FINAL DATASET - Final curated dataset (1).csv",
                stringsAsFactors = FALSE)
colnames(deg) <- make.names(colnames(deg))
deg <- deg[-1, ]

gene_col  <- "GENE.ID"
logfc_col <- "logFC..Tumour.Normal."
pval_col  <- "pVAlue"

deg$logFC  <- as.numeric(deg[[logfc_col]])
deg$plot_p <- as.numeric(deg[[pval_col]])
deg$Gene   <- deg[[gene_col]]
deg <- deg[!is.na(deg$logFC) & !is.na(deg$plot_p), ]
deg$plot_p[deg$plot_p == 0] <- 1e-300
deg$negLogP <- -log10(deg$plot_p)

deg$Significance <- "Not Significant"
deg$Significance[deg$plot_p < 0.05 & deg$logFC >  1] <- "Upregulated"
deg$Significance[deg$plot_p < 0.05 & deg$logFC < -1] <- "Downregulated"

# ===============================
# TOP 5 UP AND DOWN
# ===============================
top_up <- deg %>%
  filter(Significance == "Upregulated") %>%
  arrange(plot_p) %>%
  slice(1:5)

top_down <- deg %>%
  filter(Significance == "Downregulated") %>%
  arrange(plot_p) %>%
  slice(1:5)

top_genes <- bind_rows(top_up, top_down)

# ===============================
# 20 NETWORK GENES
# ===============================
network_genes <- c(
  "KRAS","TP53","EGFR","MET","PIK3CA","CDKN2A","MYC",
  "NOTCH1","CTNNB1","SMAD4","MDM2","BRCA1","ATM","NF1",
  "MMP11","FAP","HK2","JAK1","PTP4A1","NQO1"
)

network_in_deg <- deg %>% filter(Gene %in% network_genes)

missing <- setdiff(network_genes, network_in_deg$Gene)
if (length(missing) > 0) {
  cat("\nNetwork genes NOT found in DEG dataset:\n"); print(missing)
}

genes_to_label <- bind_rows(top_genes, network_in_deg) %>%
  distinct(Gene, .keep_all = TRUE)

# ===============================
# PLOT
# ===============================
p <- ggplot(deg, aes(x = logFC, y = negLogP, color = Significance)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  scale_color_manual(
    values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not Significant" = "grey70"
    ),
    breaks = c("Upregulated", "Downregulated", "Not Significant")
  ) +
  
  # linewidth >= 2pt → in ggplot2 'linewidth' is in mm; 2pt ≈ 0.71 mm
  geom_vline(xintercept = c(-1, 1),
             linetype  = "dashed",
             color     = "grey40",
             linewidth = 0.75) +
  geom_hline(yintercept = -log10(0.05),
             linetype  = "dashed",
             color     = "grey40",
             linewidth = 0.75) +
  
  geom_text_repel(
    data               = genes_to_label,
    aes(label          = Gene,
        color          = Significance),
    size               = 2.8,
    fontface           = "bold",
    max.overlaps       = 40,
    box.padding        = 0.4,
    point.padding      = 0.3,
    segment.alpha      = 0.7,
    segment.size       = 0.75,
    min.segment.length = 0.2,
    show.legend        = FALSE
  ) +
  
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.minor  = element_blank(),
    legend.title      = element_blank(),
    legend.position   = "right",
    plot.title        = element_blank(),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 8),
    legend.text       = element_text(size = 8)
  ) +
  
  labs(
    title = NULL,
    x     = "log2 Fold Change (Tumour vs Normal)",
    y     = expression(-log[10](p-value))
  )

# ===============================
# EXPORT — ALL THREE FORMATS
# ===============================

# TIFF (primary submission format — lossless)
ggsave(
  filename    = "Figure1.tiff",
  plot        = p,
  width       = 180,
  height      = 140,
  units       = "mm",
  dpi         = 300,
  device      = "tiff",
  compression = "lzw"
)

# JPEG (accepted alternate format)
ggsave(
  filename = "Figure1.jpg",
  plot     = p,
  width    = 180,
  height   = 140,
  units    = "mm",
  dpi      = 300,
  device   = "jpeg",
  quality  = 95
)

# PNG (for personal reference only — not a Frontiers-accepted format)
ggsave(
  filename = "Figure1.png",
  plot     = p,
  width    = 180,
  height   = 140,
  units    = "mm",
  dpi      = 300,
  device   = "png"
)

cat("Saved: Figure1.tiff, Figure1.jpg, Figure1.png\n")
cat("\nManuscript caption:\n")
cat("Figure 1. Volcano plot of differentially expressed genes in PDAC.\n")
cat("Genes with log2FC > 1 and p < 0.05 are shown in red (upregulated);\n")
cat("genes with log2FC < -1 and p < 0.05 are shown in blue (downregulated);\n")
cat("grey points are not significant. Vertical dashed lines indicate |log2FC| = 1;\n")
cat("horizontal dashed line indicates p = 0.05.\n")
cat("Selected genes of network pharmacological relevance are labelled.\n")