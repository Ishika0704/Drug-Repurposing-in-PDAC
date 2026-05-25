# =========================================
# GEO EXPRESSION PREPROCESSING PIPELINE
# + Ensembl → Gene Symbol (stable)
# =========================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
})

input_path  <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/GSE205154.csv"
output_path <- "C:/Users/Ishika/OneDrive/Desktop/Capstone Project/Data/geo_clean_symbol_matrix.csv"

# ---- Load data
geo <- read.csv(
  input_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Raw matrix:", dim(geo), "\n")

# ---- Remove junk columns
geo <- geo[, sapply(geo, function(x)
  !all(is.na(suppressWarnings(as.numeric(x))))
)]

# ---- Convert numeric
geo[] <- lapply(geo, function(x) as.numeric(as.character(x)))

# ---- Strip Ensembl versions
rownames(geo) <- sub("\\..*$", "", rownames(geo))

# ---- Map Ensembl → Gene Symbol
cat("Mapping Ensembl IDs...\n")

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(geo),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# ---- Drop unmapped
valid <- !is.na(symbols)
geo <- geo[valid, ]
symbols <- symbols[valid]

rownames(geo) <- symbols

# ---- Collapse duplicates
geo_clean <- geo %>%
  rownames_to_column("gene") %>%
  group_by(gene) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  column_to_rownames("gene")

# ---- Validate
cat("Final matrix:", dim(geo_clean), "\n")
stopifnot(
  anyDuplicated(rownames(geo_clean)) == 0,
  !any(is.na(geo_clean))
)

cat("Expression range:", range(geo_clean), "\n")

# ---- Save
write.csv(geo_clean, output_path)

cat("Done.\n")
