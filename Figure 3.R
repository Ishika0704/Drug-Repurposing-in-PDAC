# =========================================================
# FRONTIERS PUBLICATION-GRADE FIGURE 2
# PDAC Drug Repurposing Landscape — SigCom LINCS
# =========================================================

# install.packages(c("readr","dplyr","ggplot2","ggrepel","scales","forcats","patchwork"))

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(forcats)
library(patchwork)

# =========================================================
# LOAD DATA
# =========================================================

data <- read_csv(
  "C:/Users/Ishika/Downloads/SigCom - final reversers.csv",
  show_col_types = FALSE
)
colnames(data) <- make.names(colnames(data))

mimic_data <- read_csv(
  "C:/Users/Ishika/Downloads/SigCom - final mimickers.csv",
  show_col_types = FALSE
)
colnames(mimic_data) <- make.names(colnames(mimic_data))

# =========================================================
# ASSIGN COLUMNS
# z.score..sum. used directly — no bias correction
# preserves original ranking exactly
# =========================================================

data$zScore    <- as.numeric(data$z.score..sum.)
data$pValue    <- as.numeric(data$p.value..up.)

mimic_data$zScore  <- as.numeric(mimic_data$z.score..sum.)
mimic_data$pValue  <- as.numeric(mimic_data$p.value..up.)

# =========================================================
# REMOVE MISSING VALUES
# =========================================================

data <- data %>%
  filter(!is.na(zScore), !is.na(pValue), !is.na(Perturbagen))

mimic_data <- mimic_data %>%
  filter(!is.na(zScore), !is.na(pValue), !is.na(Perturbagen))

# =========================================================
# SAFE P-VALUE FLOOR
# =========================================================

data$pValue       <- pmax(data$pValue,       1e-50)
mimic_data$pValue <- pmax(mimic_data$pValue, 1e-50)

# =========================================================
# DEFINE CONTROL DRUGS
# =========================================================

control_names <- c("Curcumin", "Bosutinib", "Azacitidine")

# =========================================================
# SELECT TOP REVERSAL DRUGS (30)
# ranked by most negative z.score..sum.
# =========================================================

top_reversers <- data %>%
  arrange(zScore) %>%
  slice(1:30) %>%
  mutate(Effect = "Suppresses PDAC Signature")

# =========================================================
# SELECT TOP MIMICKER DRUGS (20)
# ranked by most positive z.score..sum.
# =========================================================

top_mimickers <- mimic_data %>%
  arrange(desc(zScore)) %>%
  slice(1:20) %>%
  mutate(Effect = "Promotes PDAC Signature")

# =========================================================
# EXTRACT CONTROLS FROM FULL UNSLICED DATA
# =========================================================

controls <- data %>%
  filter(tolower(Perturbagen) %in% tolower(control_names)) %>%
  mutate(Effect = "Suppresses PDAC Signature")

# =========================================================
# ASSIGN DIRECTIONAL SIGN FOR VISUALIZATION
# Reversers/controls → negative x
# Mimickers          → positive x
# =========================================================

top_reversers$zScore <- -abs(top_reversers$zScore)
controls$zScore      <- -abs(controls$zScore)
top_mimickers$zScore <-  abs(top_mimickers$zScore)

# =========================================================
# COMBINE & DEDUPLICATE
# =========================================================

plot_drugs <- bind_rows(
  top_mimickers,
  top_reversers,
  controls
) %>%
  distinct(Perturbagen, .keep_all = TRUE)

# =========================================================
# PDAC CLINICAL EVIDENCE LABEL
# =========================================================

plot_drugs <- plot_drugs %>%
  mutate(
    PDAC_Evidence = ifelse(
      tolower(Perturbagen) %in% tolower(control_names),
      "Tested",
      "Not Tested"
    )
  )

# =========================================================
# COMPUTE -log10(p), capped at 50
# =========================================================

plot_drugs$negLog10P <- pmin(-log10(plot_drugs$pValue), 50)

# =========================================================
# BUBBLE SIZE — based on absolute z-score magnitude
# larger score = larger bubble; clean and interpretable
# =========================================================

plot_drugs$BubbleSize <- abs(plot_drugs$zScore)

# =========================================================
# COLOR PALETTE
# =========================================================

effect_colors <- c(
  "Suppresses PDAC Signature" = "red",
  "Promotes PDAC Signature"   = "blue"
)

# =========================================================
# SHARED PUBLICATION THEME
# All text >= 8pt per Frontiers requirement
# =========================================================

publication_theme <- theme_classic(base_size = 9) +
  theme(
    axis.title      = element_text(face = "bold", size = 9),
    axis.text       = element_text(color = "black", size = 8),
    legend.title    = element_text(face = "bold", size = 8),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    plot.margin     = margin(6, 6, 6, 6),
    plot.title      = element_blank()    # removed per Frontiers rules
  )

# =========================================================
# PANEL (A) — HORIZONTAL BAR PLOT
# =========================================================

plot_drugs_bar <- plot_drugs %>%
  arrange(zScore) %>%
  mutate(Perturbagen = factor(Perturbagen, levels = Perturbagen))

panel_A <- ggplot(
  plot_drugs_bar,
  aes(x = zScore, y = Perturbagen)
) +
  
  geom_col(
    aes(fill = Effect),
    width     = 0.82,
    color     = "grey30",
    linewidth = 0.2
  ) +
  
  # Thick border overlay for PDAC-tested control drugs
  geom_col(
    data      = filter(plot_drugs_bar, PDAC_Evidence == "Tested"),
    aes(fill  = Effect),
    width     = 0.82,
    color     = "black",
    linewidth = 1.0             # ~2.8pt satisfies Frontiers >= 2pt rule
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.75,          # ~2.1pt satisfies Frontiers >= 2pt rule
    color      = "grey40"
  ) +
  
  scale_fill_manual(values = effect_colors) +
  
  scale_x_continuous(
    breaks = seq(-12, 12, 2),
    expand = expansion(mult = 0.05)
  ) +
  
  publication_theme +
  theme(
    axis.text.y     = element_text(size = 7),
    legend.position = "right"
  ) +
  
  labs(
    title = NULL,
    x     = "Connectivity Score",
    y     = "Drug Compounds",
    fill  = "Drug Effect"
  )

# =========================================================
# PANEL (B) — BUBBLE PLOT (SINGLE PANEL, REVISED)
# =========================================================

set.seed(42)

panel_B <- ggplot(
  plot_drugs,
  aes(x = zScore, y = negLog10P)
) +
  
  geom_point(
    aes(
      size  = BubbleSize,
      fill  = Effect,
      color = PDAC_Evidence
    ),
    shape  = 21,
    stroke = 0.8,
    alpha  = 0.85
  ) +
  
  # Label only the most significant drugs to reduce clutter
  # Top 8 by negLog10P from each group
  geom_text_repel(
    data = bind_rows(
      plot_drugs %>%
        filter(Effect == "Suppresses PDAC Signature") %>%
        arrange(desc(negLog10P)) %>%
        slice(1:8),
      plot_drugs %>%
        filter(Effect == "Promotes PDAC Signature") %>%
        arrange(desc(negLog10P)) %>%
        slice(1:8)
    ),
    aes(label = Perturbagen),
    size               = 2.5,
    fontface           = "bold",
    box.padding        = 0.6,
    point.padding      = 0.4,
    segment.color      = "grey50",
    segment.size       = 0.4,
    segment.curvature  = 0.1,
    force              = 10,
    force_pull         = 0.3,
    max.overlaps       = 100,
    min.segment.length = 0.1,
    max.iter           = 20000,
    # Push labels toward edges away from center
    xlim = c(-18, 18),
    ylim = c(0, 58)
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.75,
    color      = "grey30"
  ) +
  
  # Shaded background regions to distinguish suppressors vs promoters
  annotate(
    "rect",
    xmin = -Inf, xmax = 0,
    ymin = -Inf, ymax = Inf,
    fill = "red", alpha = 0.04
  ) +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "blue", alpha = 0.04
  ) +
  scale_fill_manual(values = effect_colors) +
  
  scale_color_manual(
    values = c(
      "Tested"     = "red3",
      "Not Tested" = "black"
    ),
    drop = FALSE
  ) +
  
  # Much smaller bubble size range
  scale_size_continuous(
    range  = c(1, 5),
    breaks = c(3, 6, 9, 12),
    name   = "Connectivity\nStrength"
  ) +
  
  scale_x_continuous(
    breaks = seq(-16, 16, 4),
    expand = expansion(mult = 0.08),
    limits = c(-16, 16)
  ) +
  
  scale_y_continuous(
    breaks = c(5, 10, 20, 30, 40, 50),
    expand = expansion(mult = 0.08),
    limits = c(0, 58)
  ) +
  
  publication_theme +
  theme(
    legend.position = "right"
  ) +
  
  labs(
    title = NULL,
    x     = "Connectivity Score (Reversal ← → Mimicry)",
    y     = expression(-log[10](italic(P))),
    fill  = "Drug Effect",
    color = "PDAC Evidence",
    size  = "Connectivity\nStrength"
  )
# =========================================================
# COMBINE PANELS WITH PATCHWORK
# Uppercase bold (A)(B) tags per Frontiers requirement
# =========================================================

final_figure <- (panel_A / panel_B) +
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
# 180mm wide, 260mm tall, 300dpi
# =========================================================

# TIFF — primary Frontiers submission format
ggsave(
  filename    = "C:/Users/Ishika/Downloads/Figure3.tiff",
  plot        = final_figure,
  width       = 180,
  height      = 260,
  units       = "mm",
  dpi         = 300,
  device      = "tiff",
  compression = "lzw",
  bg          = "white"
)

# JPEG — accepted alternate format
ggsave(
  filename = "C:/Users/Ishika/Downloads/Figure3.jpg",
  plot     = final_figure,
  width    = 180,
  height   = 260,
  units    = "mm",
  dpi      = 300,
  device   = "jpeg",
  quality  = 95,
  bg       = "white"
)

# PNG — personal reference only
ggsave(
  filename = "C:/Users/Ishika/Downloads/Figure3.png",
  plot     = final_figure,
  width    = 180,
  height   = 260,
  units    = "mm",
  dpi      = 300,
  device   = "png",
  bg       = "white"
)

cat("\n✅ Saved: Figure2.tiff, Figure2.jpg, Figure2.png\n")
cat("\n📝 Manuscript caption:\n")
cat("Figure 2. Drug perturbation landscape for PDAC using SigCom LINCS connectivity analysis.\n")
cat("(A) Horizontal bar chart of connectivity scores for top candidate drugs.\n")
cat("    Orange bars: PDAC signature suppressors; blue bars: signature promoters.\n")
cat("    Thick-bordered bars indicate drugs with prior PDAC clinical evidence.\n")
cat("(B) Connectivity-significance bubble plot stratified by drug effect.\n")
cat("    Bubble size reflects absolute connectivity strength; border color indicates\n")
cat("    PDAC clinical evidence (red = tested, black = not tested).\n")
cat("    Dashed vertical line marks connectivity score = 0.\n")