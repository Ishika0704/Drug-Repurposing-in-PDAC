##############################
# 1️⃣ Load Libraries
##############################

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(forcats)
library(patchwork)

##############################
# 2️⃣ Load Reversal Data
##############################

data <- read_csv(
  "C:/Users/Ishika/Downloads/Reversers-iLINCS.csv",
  show_col_types = FALSE
)

colnames(data) <- make.names(colnames(data))

##############################
# 3️⃣ Load Mimicker Data
##############################

mimic_data <- read_csv(
  "C:/Users/Ishika/Downloads/Mimickers-iLINCS.csv",
  show_col_types = FALSE
)

colnames(mimic_data) <- make.names(colnames(mimic_data))

##############################
# 4️⃣ Fix Zero p-values
##############################

data$pValue[data$pValue == 0]             <- 1e-300
mimic_data$pValue[mimic_data$pValue == 0] <- 1e-300

##############################
# 5️⃣ Bias Correction
##############################

data$zScore_adjusted       <- data$zScore       / sqrt(data$NoOfSignatures)
mimic_data$zScore_adjusted <- mimic_data$zScore / sqrt(mimic_data$NoOfSignatures)

##############################
# 6️⃣ Define Control Drugs
##############################

control_names <- c("Curcumin", "Bosutinib", "Azacitidine")

##############################
# 7️⃣ Filter Robust Candidates
##############################

robust_hits  <- data       %>% filter(NoOfSignatures >= 5)
robust_mimic <- mimic_data %>% filter(NoOfSignatures >= 5)

##############################
# 8️⃣ Select Top Reversal Drugs
##############################

top_reversers <- robust_hits %>%
  arrange(desc(zScore_adjusted)) %>%
  slice(1:30) %>%
  mutate(Effect = "Suppresses PDAC Signature")

##############################
# 9️⃣ Extract Controls
##############################

controls <- data %>%
  filter(tolower(Perturbagen) %in% tolower(control_names)) %>%
  mutate(Effect = "Suppresses PDAC Signature")

##############################
# 🔟 Select Top Mimickers
##############################

top_mimickers <- robust_mimic %>%
  arrange(desc(zScore_adjusted)) %>%
  slice(1:20) %>%
  mutate(Effect = "Promotes PDAC Signature")

##############################
# 1️⃣1️⃣ Assign Visualization Sign
##############################

top_reversers$zScore_adjusted <- -abs(top_reversers$zScore_adjusted)
controls$zScore_adjusted      <- -abs(controls$zScore_adjusted)
top_mimickers$zScore_adjusted <-  abs(top_mimickers$zScore_adjusted)

##############################
# 1️⃣2️⃣ Combine & Deduplicate
##############################

plot_drugs <- bind_rows(
  top_mimickers,
  top_reversers,
  controls
) %>%
  arrange(desc(Effect)) %>%
  distinct(Perturbagen, .keep_all = TRUE)

##############################
# 1️⃣3️⃣ Mark PDAC-Tested Controls
##############################

plot_drugs <- plot_drugs %>%
  mutate(
    PDAC_Evidence = ifelse(
      tolower(Perturbagen) %in% tolower(control_names),
      "Tested",
      "Not Tested"
    )
  )

##############################
# 1️⃣4️⃣ Compute -log10(p)
##############################

plot_drugs$negLog10P <- pmin(-log10(plot_drugs$pValue), 50)

##############################
# 1️⃣5️⃣ Verify scales version
##############################

if (packageVersion("scales") < "1.2.0") {
  stop(
    "scales >= 1.2.0 required for pseudo_log_trans().\n",
    "Run: install.packages('scales')"
  )
}

################################################
# 1️⃣6️⃣ Bubble Plot
# GUIDELINE FIXES:
# - smaller readable fonts
# - thinner lines
# - panel-ready formatting
# - no oversized text
################################################

bubble_plot <- ggplot(
  plot_drugs,
  aes(
    x = zScore_adjusted,
    y = negLog10P
  )
) +
  
  geom_point(
    aes(
      size  = NoOfSignatures,
      fill  = Effect,
      color = PDAC_Evidence
    ),
    shape  = 21,
    stroke = 0.8,
    alpha  = 0.9
  ) +
  
  geom_text_repel(
    aes(label = Perturbagen),
    size          = 2.8,
    max.overlaps  = 40,
    segment.color = "grey50",
    segment.size  = 0.3,
    box.padding   = 0.3,
    point.padding = 0.2
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.5
  ) +
  
  scale_x_continuous(
    breaks = seq(-4, 4, 1),
    expand = expansion(mult = 0.05)
  ) +
  
  scale_y_continuous(
    trans  = pseudo_log_trans(base = 10),
    breaks = c(3, 5, 10, 20, 50)
  ) +
  
  scale_fill_manual(values = c(
    "Suppresses PDAC Signature" = "firebrick",
    "Promotes PDAC Signature"   = "steelblue"
  )) +
  
  scale_color_manual(values = c(
    "Tested"     = "red",
    "Not Tested" = "black"
  )) +
  
  scale_size(
    range = c(2, 8),
    trans = "sqrt"
  ) +
  
  theme_classic(base_size = 9) +
  
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 8),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    
    axis.title       = element_text(size = 9, face = "bold"),
    axis.text        = element_text(size = 8),
    
    plot.title       = element_text(
      face = "bold",
      hjust = 0.5,
      size = 10
    )
  ) +
  
  labs(
    title = NULL,
    x     = "Connectivity Score",
    y     = expression(-log[10](p~value)),
    fill  = "Drug Effect",
    color = "PDAC Clinical Evidence",
    size  = "LINCS Signatures"
  )

################################################
# 1️⃣7️⃣ Horizontal Connectivity Bar Plot
# GUIDELINE FIXES:
# - reduced font size
# - thinner lines
# - panel consistency
################################################

plot_drugs_bar <- plot_drugs %>%
  arrange(zScore_adjusted) %>%
  mutate(
    Perturbagen = factor(Perturbagen, levels = Perturbagen)
  )

horizontal_plot <- ggplot(
  plot_drugs_bar,
  aes(
    x = zScore_adjusted,
    y = Perturbagen
  )
) +
  
  geom_col(
    aes(fill = Effect),
    width = 0.9,
    color = "grey30",
    linewidth = 0.2
  ) +
  
  geom_col(
    data = filter(plot_drugs_bar, PDAC_Evidence == "Tested"),
    aes(fill = Effect),
    width     = 0.9,
    color     = "black",
    linewidth = 0.8
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.5
  ) +
  
  scale_fill_manual(values = c(
    "Suppresses PDAC Signature" = "#c41e1e",
    "Promotes PDAC Signature"   = "steelblue"
  )) +
  
  theme_classic(base_size = 9) +
  
  theme(
    plot.title   = element_text(
      hjust = 0.5,
      face = "bold",
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 9,
      face = "bold"
    ),
    
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 6.5),
    
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    
    legend.text  = element_text(size = 7),
    
    legend.key.size = unit(0.45, "cm")
  ) +
  
  labs(
    title = NULL,
    x     = "Connectivity Score",
    y     = "Drug Compounds",
    fill  = "Drug Effect"
  )

################################################
# 1️⃣8️⃣ Combine into Panel Figure
# Frontiers-compatible panel layout
################################################

final_panel <- (
  horizontal_plot /
    bubble_plot
) +
  
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.tag.position = c(0, 1)
  )

################################################
# 1️⃣9️⃣ Display Panel
################################################

print(final_panel)

################################################
# 2️⃣0️⃣ Save Panel Figure
# Frontiers-ready
################################################

ggsave(
  filename = "C:/Users/Ishika/Downloads/Figure2_iLINCS_panel.tiff",
  plot     = final_panel,
  
  width    = 180,
  height   = 230,
  
  units    = "mm",
  
  dpi      = 300,
  
  compression = "lzw",
  
  bg       = "white"
)

ggsave(
  filename = "C:/Users/Ishika/Downloads/Figure2_iLINCS_panel.png",
  plot     = final_panel,
  
  width    = 180,
  height   = 230,
  
  units    = "mm",
  
  dpi      = 300,
  
  bg       = "white"
)

message("✅ Saved: Figure2_iLINCS_panel")
