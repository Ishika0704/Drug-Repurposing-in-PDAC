# =========================================================
# FRONTIERS PUBLICATION-GRADE FIGURE 3
# Docking Poses + Binding Energy Box Plots
# Panel (A): Docking poses image
# Panels (B)-(F): Binding energy box plots per protein
# =========================================================

# install.packages(c("tidyverse","patchwork","png","grid"))

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(tidyverse)
library(patchwork)
library(png)
library(grid)

# =========================================================
# SHARED SETTINGS
# =========================================================

protein_order <- c("MMP11", "PTP4A1", "FAP", "NQO1", "HK2")

# Shared drug color palette — same across all panels
drug_colors <- c(
  "Oxandrolone"  = "#E69F00",
  "Netupitant"   = "#56B4E9",
  "BGT-226"      = "#009E73",
  "Tepotinib"    = "#CC79A7",
  "Carfilzomib"  = "#D55E00"
)

# =========================================================
# SHARED PUBLICATION THEME
# All text >= 8pt per Frontiers requirement
# =========================================================

publication_theme <- theme_classic(base_size = 9) +
  theme(
    plot.title      = element_blank(),       # removed per Frontiers rules
    axis.title      = element_text(face = "bold", size = 9),
    axis.text.x     = element_text(
      angle = 45,
      hjust = 1,
      size  = 8
    ),
    axis.text.y     = element_text(size = 8),
    legend.position = "none",                # no legend; colors consistent
    plot.margin     = margin(4, 4, 4, 4)
  )

# =========================================================
# FUNCTION: Load and clean one CSV
# =========================================================

load_docking_data <- function(filepath, protein_name) {
  
  df <- read.csv(filepath, header = TRUE)
  
  # Rename first column to Drug
  colnames(df)[1] <- "Drug"
  
  # Remove junk first row if it exists
  if (is.na(suppressWarnings(as.numeric(df$Drug[1]))) &&
      df$Drug[1] %in% c("Drug", "drug", "")) {
    df <- df[-1, ]
  }
  
  rownames(df) <- NULL
  
  # Convert all columns except Drug to numeric
  df[-1] <- lapply(df[-1], function(x) as.numeric(as.character(x)))
  
  # Wide to long
  df_long <- df %>%
    pivot_longer(
      cols      = -Drug,
      names_to  = "Pose",
      values_to = "Binding_Energy"
    ) %>%
    na.omit() %>%
    mutate(Protein = protein_name)
  
  return(df_long)
}

# =========================================================
# FUNCTION: Build one box plot panel
# =========================================================

make_boxplot <- function(df_long, protein_name) {
  
  ggplot(
    df_long,
    aes(
      x    = reorder(Drug, Binding_Energy, median),
      y    = Binding_Energy,
      fill = Drug
    )
  ) +
    
    geom_boxplot(
      outlier.color = "red",
      outlier.size  = 1.2,
      outlier.shape = 16,
      linewidth     = 0.75,    # ~2.1pt — satisfies Frontiers >= 2pt rule
      alpha         = 0.80
    ) +
    
    geom_jitter(
      width = 0.18,
      size  = 0.8,
      alpha = 0.65,
      color = "grey25"
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.75,       # ~2.1pt — satisfies Frontiers >= 2pt rule
      color      = "grey50"
    ) +
    
    scale_fill_manual(values = drug_colors) +
    
    # Protein name as subtitle inside plot area
    annotate(
      "text",
      x        = Inf,
      y        = Inf,
      label    = protein_name,
      hjust    = 1.1,
      vjust    = 1.5,
      fontface = "bold",
      size     = 3.2,          # ~9pt
      color    = "grey20"
    ) +
    
    publication_theme +
    
    labs(
      title = NULL,
      x     = "Drug",
      y     = "Binding Energy\n(kcal/mol)"
    )
}

# =========================================================
# LOAD ALL FIVE DATASETS
# =========================================================

data_MMP11  <- load_docking_data(
  "C:/Users/Ishika/Downloads/MMP11-iLINCS.csv",  "MMP11")
data_PTP4A1 <- load_docking_data(
  "C:/Users/Ishika/Downloads/PTP4A1-iLINCS.csv", "PTP4A1")
data_FAP    <- load_docking_data(
  "C:/Users/Ishika/Downloads/FAP-iLINCS.csv",    "FAP")
data_NQO1   <- load_docking_data(
  "C:/Users/Ishika/Downloads/NQO1-iLINCS.csv",   "NQO1")
data_HK2    <- load_docking_data(
  "C:/Users/Ishika/Downloads/HK2-iLINCS.csv",    "HK2")

# =========================================================
# BUILD FIVE BOX PLOT PANELS (B)-(F)
# =========================================================

panel_B <- make_boxplot(data_MMP11,  "MMP11")
panel_C <- make_boxplot(data_PTP4A1, "PTP4A1")
panel_D <- make_boxplot(data_FAP,    "FAP")
panel_E <- make_boxplot(data_NQO1,   "NQO1")
panel_F <- make_boxplot(data_HK2,    "HK2")

# =========================================================
# PANEL (A) — Load docking poses image
# =========================================================

docking_img <- readPNG(
  "C:/Users/Ishika/Downloads/Figure 4 (a) (1).png"
)

panel_A <- ggplot() +
  annotation_custom(
    rasterGrob(
      docking_img,
      width  = unit(1, "npc"),
      height = unit(1, "npc")
    )
  ) +
  theme_void() +          # completely blank background
  theme(
    plot.margin = margin(2, 2, 2, 2)
  )

# =========================================================
# COMBINE ALL PANELS WITH PATCHWORK
#
# Layout:
#   (A) docking poses — full width top
#   (B)(C)(D) — three box plots middle row
#   (E)(F)    — two box plots bottom row (centred)
# =========================================================

boxplot_row1 <- panel_B | panel_C | panel_D
boxplot_row2 <- panel_E | panel_F

final_figure <- panel_A /
  boxplot_row1 /
  boxplot_row2 +
  
  plot_layout(
    heights = c(2.2, 1, 1)   # panel A taller; box plots equal height
  ) +
  
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 9)
  )

# =========================================================
# DISPLAY
# =========================================================

print(final_figure)

# =========================================================
# EXPORT — ALL THREE FORMATS
# 180mm wide, 300dpi per Frontiers spec
# Height set to 280mm to fit all panels on one page
# =========================================================

# TIFF — primary Frontiers submission format
ggsave(
  filename    = "C:/Users/Ishika/Downloads/Figure4.tiff",
  plot        = final_figure,
  width       = 180,
  height      = 280,
  units       = "mm",
  dpi         = 300,
  device      = "tiff",
  compression = "lzw",
  bg          = "white"
)

# JPEG — accepted alternate format
ggsave(
  filename = "C:/Users/Ishika/Downloads/Figure4.jpg",
  plot     = final_figure,
  width    = 180,
  height   = 280,
  units    = "mm",
  dpi      = 300,
  device   = "jpeg",
  quality  = 95,
  bg       = "white"
)

# PNG — personal reference only
ggsave(
  filename = "C:/Users/Ishika/Downloads/Figure4.png",
  plot     = final_figure,
  width    = 180,
  height   = 280,
  units    = "mm",
  dpi      = 300,
  device   = "png",
  bg       = "white"
)

cat("\n✅ Saved: Figure3.tiff, Figure3.jpg, Figure3.png\n")
cat("\n📝 Manuscript caption:\n")
cat("Figure 3. Molecular docking analysis of shortlisted drug candidates against PDAC-relevant protein targets.\n")
cat("(A) Docking poses of five candidate drugs (rows) against five protein targets (columns).\n")
cat("    Protein structures shown in green ribbon representation; ligands shown in cyan stick representation.\n")
cat("(B)-(F) Binding energy distributions (kcal/mol) for each drug across all docking poses\n")
cat("    against MMP11, PTP4A1, FAP, NQO1, and HK2 respectively.\n")
cat("    Box plots show median, interquartile range, and individual pose energies (grey points).\n")
cat("    Red points indicate outlier poses. Dashed line marks 0 kcal/mol reference.\n")