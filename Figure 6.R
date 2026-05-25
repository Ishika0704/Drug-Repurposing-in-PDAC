# =========================================================
# FRONTIERS PUBLICATION-GRADE FIGURE 6
# UNIFIED PDAC NETWORK FIGURE
#
# PANEL (A):
# Gene–Gene Interaction Network
#
# PANEL (B):
# Drug–Gene Interaction Network
#
# IMPORTANT:
# - Preserves ALL biological logic
# - Preserves ALL layouts
# - Preserves ALL colors
# - Preserves ALL edge styles
# - Preserves ALL calculations
#
# ONLY improves:
# • panel organization
# • typography
# • spacing
# • legend readability
# • Frontiers compliance
# • export quality
# =========================================================

rm(list = ls())

.libPaths("C:/Users/Ishika/Documents/R/win-library/4.5")

# =========================================================
# LIBRARIES
# =========================================================

library(STRINGdb)
library(igraph)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(httr)

set.seed(42)

# =========================================================
# FIX STRING DOWNLOAD TIMEOUT
# =========================================================

options(timeout = 600)

download.file.method <- "libcurl"

options(
  download.file.method = download.file.method
)

# =========================================================
# STEP 1: INPUT GENES
# =========================================================

genes <- c(
  "KRAS","TP53","EGFR","MET","PIK3CA","CDKN2A","MYC",
  "NOTCH1","CTNNB1","SMAD4","MDM2","BRCA1","ATM","NF1",
  "MMP11","FAP","HK2","JAK1","PTP4A1","NQO1"
)

gene_df <- data.frame(
  gene = genes
)

# =========================================================
# STEP 2: STRING MAPPING
# =========================================================

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 400
)

mapped <- string_db$map(
  gene_df,
  "gene",
  removeUnmappedRows = TRUE
)

cat("\nMapped genes:\n")
print(mapped$gene)

# =========================================================
# STEP 3: BUILD NETWORK VIA STRING API
# =========================================================

identifiers <- paste(
  mapped$gene,
  collapse = "%0d"
)

url <- paste0(
  "https://string-db.org/api/tsv/network?",
  "identifiers=", identifiers,
  "&species=9606",
  "&required_score=400",
  "&caller_identity=PDAC_capstone"
)

response <- GET(url)

ppi_raw <- read.delim(
  text = content(response, "text"),
  sep = "\t"
)

cat("\nColumns returned by API:\n")
print(colnames(ppi_raw))

cat("\nEdges found:", nrow(ppi_raw), "\n")

ppi_filtered <- ppi_raw[
  ppi_raw$preferredName_A %in% mapped$gene &
    ppi_raw$preferredName_B %in% mapped$gene,
]

ppi_filtered <- data.frame(
  from = ppi_filtered$preferredName_A,
  to   = ppi_filtered$preferredName_B
)

ppi_filtered <- ppi_filtered[
  !is.na(ppi_filtered$from) &
    !is.na(ppi_filtered$to),
]

cat("Edges after filtering:", nrow(ppi_filtered), "\n")

if (nrow(ppi_filtered) == 0) {
  stop("No edges remain.")
}

# =========================================================
# STEP 4: CREATE GENE NETWORK
# =========================================================

network_gene <- graph_from_data_frame(
  ppi_filtered,
  directed = FALSE
)

isolated_idx <- which(
  degree(network_gene) == 0
)

if (length(isolated_idx) > 0) {
  
  cat("\nRemoved isolated genes:\n")
  
  print(
    V(network_gene)$name[isolated_idx]
  )
  
  network_gene <- delete_vertices(
    network_gene,
    isolated_idx
  )
}

# =========================================================
# STEP 5: ENRICHMENT
# =========================================================

gene_entrez <- bitr(
  V(network_gene)$name,
  
  fromType = "SYMBOL",
  
  toType = "ENTREZID",
  
  OrgDb = org.Hs.eg.db
)

kegg <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  
  organism = "hsa",
  
  pvalueCutoff = 0.05
)

kegg_df <- as.data.frame(kegg)

cat("\nTop KEGG pathways:\n")

print(
  head(kegg_df$Description, 10)
)

# =========================================================
# STEP 6: MAP PATHWAYS
# =========================================================

gene_pathway_map <- data.frame()

for (i in seq_len(nrow(kegg_df))) {
  
  entrez_ids <- unlist(
    strsplit(
      kegg_df$geneID[i],
      "/"
    )
  )
  
  gene_symbols <- gene_entrez$SYMBOL[
    match(
      entrez_ids,
      gene_entrez$ENTREZID
    )
  ]
  
  gene_symbols <- gene_symbols[
    !is.na(gene_symbols)
  ]
  
  gene_symbols <- gene_symbols[
    gene_symbols %in%
      V(network_gene)$name
  ]
  
  if (length(gene_symbols) == 0)
    next
  
  gene_pathway_map <- rbind(
    gene_pathway_map,
    
    data.frame(
      gene = gene_symbols,
      
      pathway = kegg_df$Description[i],
      
      stringsAsFactors = FALSE
    )
  )
}

# =========================================================
# STEP 7: ASSIGN PATHWAYS
# =========================================================

assign_pathway <- function(gene) {
  
  pw <- gene_pathway_map$pathway[
    gene_pathway_map$gene == gene
  ]
  
  if (length(pw) == 0)
    return("Unassigned")
  
  pw_lower <- tolower(pw)
  
  if (any(grepl(
    "glycolysis|gluconeogenesis|carbon metabolism|hif|hypoxia|pentose",
    pw_lower
  ))) {
    
    return(
      "Tumor metabolic reprogramming"
    )
  }
  
  if (any(grepl(
    "mapk|pi3k|erbb|ras signaling|vegf|jak-stat|mtor",
    pw_lower
  ))) {
    
    return(
      "Growth signaling (MAPK/PI3K)"
    )
  }
  
  if (any(grepl(
    "cell cycle|p53|apoptosis|senescence|dna repair|dna replication",
    pw_lower
  ))) {
    
    return(
      "Cell cycle / apoptosis"
    )
  }
  
  return(
    "Cell survival / stress signaling"
  )
}

V(network_gene)$pathway <- sapply(
  V(network_gene)$name,
  assign_pathway
)

# =========================================================
# STEP 8: COLORS
# =========================================================

path_colors <- c(
  "Tumor metabolic reprogramming"    = "#FF9AA2",
  "Growth signaling (MAPK/PI3K)"     = "#85E0C9",
  "Cell cycle / apoptosis"           = "#90C8F0",
  "Cell survival / stress signaling" = "#CBA8E8",
  "Unassigned"                       = "#CCCCCC"
)

V(network_gene)$color <- path_colors[
  V(network_gene)$pathway
]

# =========================================================
# STEP 9: GENE NETWORK LAYOUT
# =========================================================

layout_gene <- layout_with_fr(
  network_gene,
  
  niter = 3000,
  
  grid = "nogrid"
)

layout_gene <- layout_gene * 2.6

# =========================================================
# STEP 10: CREATE DRUG NETWORK
# =========================================================

network_drug <- network_gene

drug_data <- read.csv(
  "C:/Users/Ishika/Downloads/drug_gene_data.csv",
  stringsAsFactors = FALSE
)

colnames(drug_data) <- c(
  "drug",
  "gene"
)

drug_data <- na.omit(drug_data)

drug_data <- drug_data %>%
  filter(
    gene %in%
      V(network_gene)$name
  )

cat("\nDrugs retained:\n")

print(unique(drug_data$drug))

# =========================================================
# STEP 11: ADD DRUG NODES
# =========================================================

n_genes <- vcount(network_drug)

drug_nodes <- unique(
  drug_data$drug
)

network_drug <- add_vertices(
  network_drug,
  
  nv = length(drug_nodes),
  
  name = drug_nodes,
  
  color = "#FFD97D",
  
  type = "drug",
  
  pathway = "Drug",
  
  shape = "square",
  
  size = 9
)

V(network_drug)$type[
  seq_len(n_genes)
] <- "gene"

V(network_drug)$shape[
  seq_len(n_genes)
] <- "circle"

V(network_drug)$size[
  seq_len(n_genes)
] <- 11

# =========================================================
# STEP 12: ADD DRUG EDGES
# =========================================================

for (i in seq_len(nrow(drug_data))) {
  
  d_idx <- which(
    V(network_drug)$name ==
      drug_data$drug[i]
  )
  
  g_idx <- which(
    V(network_drug)$name ==
      drug_data$gene[i]
  )
  
  if (
    length(d_idx) == 1 &&
    length(g_idx) == 1
  ) {
    
    network_drug <- add_edges(
      network_drug,
      c(d_idx, g_idx)
    )
  }
}

# =========================================================
# STEP 13: EDGE STYLING
# =========================================================

n_ppi <- nrow(ppi_filtered)

m <- length(E(network_drug))

E(network_drug)$edge_type <- "ppi"

if (m > n_ppi) {
  
  E(network_drug)$edge_type[
    (n_ppi + 1):m
  ] <- "drug"
}

E(network_drug)$color <- ifelse(
  E(network_drug)$edge_type == "drug",
  "#E0A020AA",
  "#AAAAAA66"
)

E(network_drug)$width <- ifelse(
  E(network_drug)$edge_type == "drug",
  1.6,
  0.7
)

E(network_drug)$lty <- ifelse(
  E(network_drug)$edge_type == "drug",
  2,
  1
)

# =========================================================
# STEP 14: DRUG NETWORK LAYOUT
# =========================================================

total_nodes <- vcount(network_drug)

layout_fixed <- matrix(
  0,
  nrow = total_nodes,
  ncol = 2
)

rownames(layout_fixed) <- V(network_drug)$name

cluster_centers <- list(
  "Tumor metabolic reprogramming"    = c( 0,    2.8),
  "Growth signaling (MAPK/PI3K)"     = c(-2.8,  0),
  "Cell cycle / apoptosis"           = c( 2.8,  0),
  "Cell survival / stress signaling" = c( 0,   -2.8),
  "Unassigned"                       = c( 0,    0)
)

for (pw in names(cluster_centers)) {
  
  nodes <- V(network_drug)$name[
    V(network_drug)$pathway == pw &
      V(network_drug)$type == "gene"
  ]
  
  n <- length(nodes)
  
  if (n == 0)
    next
  
  center <- cluster_centers[[pw]]
  
  angles <- seq(
    0,
    2 * pi,
    length.out = n + 1
  )[-(n + 1)]
  
  radius <- ifelse(
    n == 1,
    0,
    0.9
  )
  
  layout_fixed[nodes, 1] <-
    center[1] +
    radius * cos(angles)
  
  layout_fixed[nodes, 2] <-
    center[2] +
    radius * sin(angles)
}

nd <- length(drug_nodes)

dangles <- seq(
  0,
  2 * pi,
  length.out = nd + 1
)[-1]

layout_fixed[drug_nodes, 1] <-
  5 * cos(dangles)

layout_fixed[drug_nodes, 2] <-
  5 * sin(dangles)

# =========================================================
# STEP 15: EXPORT FINAL FIGURE
# =========================================================

tiff(
  "C:/Users/Ishika/Downloads/Figure6_FINAL_Frontiers.tiff",
  
  width  = 180,
  height = 240,
  
  units  = "mm",
  
  res    = 300,
  
  compression = "lzw",
  
  bg = "white"
)

# =========================================================
# PANEL LAYOUT
# =========================================================

layout(
  matrix(
    c(1,2,
      3,4),
    nrow = 2,
    byrow = TRUE
  ),
  
  widths  = c(4.8, 1.8),
  heights = c(1,1)
)

# =========================================================
# PANEL (A) — GENE NETWORK
# =========================================================

par(
  mar = c(2.5,2.5,3,0.5),
  xpd = FALSE,
  bg  = "white"
)

plot(
  network_gene,
  
  layout = layout_gene,
  
  vertex.label = V(network_gene)$name,
  
  vertex.label.cex = 0.9,
  
  vertex.label.color = "#111111",
  
  vertex.label.font = 2,
  
  vertex.size = 18,
  
  vertex.shape = "circle",
  
  vertex.color = V(network_gene)$color,
  
  vertex.frame.color = "#BBBBBB",
  
  edge.color = "#66666688",
  
  edge.width = 0.9,
  
  rescale = TRUE
)

title(
  main = NULL,
  
  cex.main = 1.15,
  
  font.main = 2
)

mtext(
  "(A)",
  
  side = 3,
  
  line = 1,
  
  adj = 0,
  
  cex = 1.1,
  
  font = 2
)

mtext(
  "Nodes coloured by pathway  |  Edges = STRING PPI interactions (score ≥ 400)",
  
  side = 1,
  
  line = 0.6,
  
  cex = 0.75,
  
  col = "#555555"
)

# =========================================================
# PANEL (A) LEGEND
# =========================================================

par(
  mar = c(2,0,3,2),
  bg  = "white"
)

plot.new()

legend(
  "center",
  
  legend = names(path_colors),
  
  pch = 21,
  
  pt.bg = path_colors,
  
  pt.cex = 1.9,
  
  cex = 0.9,
  
  bty = "n",
  
  title = "Pathway",
  
  title.font = 2,
  
  text.col = "#222222",
  
  y.intersp = 1.2
)

# =========================================================
# PANEL (B) — DRUG NETWORK
# =========================================================

par(
  mar = c(2.5,2.5,3,0.5),
  xpd = FALSE,
  bg  = "white"
)

plot(
  network_drug,
  
  layout = layout_fixed,
  
  vertex.label = V(network_drug)$name,
  
  vertex.label.cex = ifelse(
    V(network_drug)$type == "drug",
    0.72,
    0.88
  ),
  
  vertex.label.color = "#111111",
  
  vertex.label.font = 2,
  
  vertex.size = V(network_drug)$size * 1.5,
  
  vertex.shape = V(network_drug)$shape,
  
  vertex.color = V(network_drug)$color,
  
  vertex.frame.color = "#BBBBBB",
  
  edge.color = E(network_drug)$color,
  
  edge.width = E(network_drug)$width * 1.2,
  
  edge.lty = E(network_drug)$lty,
  
  rescale = TRUE
)

title(
  main = NULL,
  
  cex.main = 1.15,
  
  font.main = 2
)

mtext(
  "(B)",
  
  side = 3,
  
  line = 1,
  
  adj = 0,
  
  cex = 1.1,
  
  font = 2
)

mtext(
  "Circles = genes  |  Squares = drugs  |  Solid edges = PPI  |  Dashed edges = drug–gene",
  
  side = 1,
  
  line = 0.6,
  
  cex = 0.75,
  
  col = "#555555"
)

# =========================================================
# PANEL (B) LEGEND
# =========================================================

par(
  mar = c(2,0,3,2),
  bg  = "white"
)

plot.new()

legend(
  "center",
  
  legend = c(
    names(path_colors),
    "Drug node"
  ),
  
  pch = c(
    rep(21, length(path_colors)),
    22
  ),
  
  pt.bg = c(
    path_colors,
    "#FFD97D"
  ),
  
  pt.cex = 1.9,
  
  cex = 0.9,
  
  bty = "n",
  
  title = "Pathway / Node Type",
  
  title.font = 2,
  
  text.col = "#222222",
  
  y.intersp = 1.2
)

dev.off()

cat("\n====================================\n")
cat("FINAL FIGURE 6 SAVED SUCCESSFULLY\n")
cat("====================================\n")

cat("\nSaved as:\n")
cat("Figure6_FINAL_Frontiers.tiff\n")
