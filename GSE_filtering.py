import pandas as pd

# ============================================================
# USER INPUT FILE PATHS (EDIT ONLY IF NEEDED)
# ============================================================

SERIES_MATRIX_PATH = r"C:\Users\Ishika\Downloads\GSE205154_series_matrix.txt"
EXPRESSION_MATRIX_PATH = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_log2TPM_noZero_TCGA_like.csv"

OUTPUT_EXPR_PATH = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GSE205154_pancreas_primary_expression.csv"
OUTPUT_META_PATH = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GSE205154_pancreas_primary_metadata.csv"
OUTPUT_REMOVED_PATH = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GSE205154_removed_samples.csv"

# ============================================================
# STEP 1: READ GEO SERIES MATRIX AS TEXT
# ============================================================

with open(SERIES_MATRIX_PATH, "r", encoding="utf-8") as f:
    lines = f.readlines()

def extract_sample_line(prefix):
    """
    Extract values from a GEO !Sample_* line.
    Returns a list of values in correct sample order.
    """
    matches = [l for l in lines if l.startswith(prefix)]
    if len(matches) == 0:
        raise ValueError(f"{prefix} not found in GEO file")

    values = matches[0].strip().split("\t")[1:]
    return [v.strip().strip('"') for v in values]

# ============================================================
# STEP 2: EXTRACT REQUIRED METADATA
# ============================================================

sample_titles = extract_sample_line("!Sample_title")
gsm_ids = extract_sample_line("!Sample_geo_accession")
tissue_source = extract_sample_line("!Sample_source_name_ch1")
tumor_type_raw = extract_sample_line("!Sample_characteristics_ch1")

# Build metadata table
meta_df = pd.DataFrame({
    "sample_title": sample_titles,
    "GSM_ID": gsm_ids,
    "tissue_source": tissue_source,
    "tumor_type_raw": tumor_type_raw
})

print("Metadata table shape:", meta_df.shape)

# ============================================================
# STEP 3: FILTER METADATA (Pancreas + Primary tumors)
# ============================================================

filtered_meta = meta_df[
    (meta_df["tissue_source"].str.lower() == "pancreas") &
    (meta_df["tumor_type_raw"].str.contains("tumor type: Primary", case=False))
].copy()

print("Samples after Pancreas + Primary filtering:", filtered_meta.shape[0])

# ============================================================
# STEP 4: LOAD EXPRESSION MATRIX
# ============================================================

# Rows = genes
# Columns = Sample Titles (ST-xxxx)
expr_df = pd.read_csv(EXPRESSION_MATRIX_PATH, index_col=0)

print("Expression matrix shape BEFORE filtering:", expr_df.shape)

# ============================================================
# STEP 5: FILTER EXPRESSION MATRIX USING SAMPLE TITLES
# ============================================================

valid_titles = filtered_meta["sample_title"].tolist()

expr_filtered = expr_df.loc[:, expr_df.columns.isin(valid_titles)]

print("Expression matrix shape AFTER filtering:", expr_filtered.shape)

# ============================================================
# STEP 6: SAVE OUTPUT FILES
# ============================================================

# Save filtered expression matrix
expr_filtered.to_csv(OUTPUT_EXPR_PATH)

# Save filtered metadata
filtered_meta.to_csv(OUTPUT_META_PATH, index=False)

# Save removed samples list
removed_samples = expr_df.columns.difference(expr_filtered.columns)
pd.DataFrame({"removed_sample_title": removed_samples}).to_csv(
    OUTPUT_REMOVED_PATH, index=False
)

print("\nFILES SAVED SUCCESSFULLY:")
print("✔ Filtered expression matrix")
print("✔ Filtered metadata")
print("✔ Removed samples list")

# ============================================================
# STEP 7: FINAL SANITY CHECKS
# ============================================================

assert expr_filtered.shape[1] == filtered_meta.shape[0], \
    "Mismatch between metadata and expression sample counts!"

print("\nSanity check passed: sample counts match.")
