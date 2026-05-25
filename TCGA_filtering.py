import pandas as pd

# -------------------------
# 1. Load TCGA expression data
# -------------------------
# Rows = genes
# Columns = TCGA sample barcodes

expr_df = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA_noZero_genes.csv",
    index_col=0
)

print("Original shape (genes × samples):", expr_df.shape)

# -------------------------
# 2. Identify tumor and normal samples
# -------------------------
tumor_mask = expr_df.columns.str.endswith("01")
normal_mask = expr_df.columns.str.endswith("11")

tumor_samples = expr_df.columns[tumor_mask]
normal_samples = expr_df.columns[normal_mask]

# -------------------------
# 3. Keep only tumor samples
# -------------------------
expr_tumor = expr_df[tumor_samples]

print("After keeping only '-01' samples:", expr_tumor.shape)

# -------------------------
# 4. Save filtered data
# -------------------------
expr_tumor.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA_primarytumor_only.csv"
)

print("Saved tumor-only TCGA data successfully.")

# -------------------------
# 5. Save summary table
# -------------------------
summary = pd.DataFrame({
    "Category": ["Original samples", "Tumor kept", "Normals removed"],
    "Count": [
        expr_df.shape[1],
        len(tumor_samples),
        len(normal_samples)
    ]
})

summary.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\TCGA_sample_filtering_summary.csv",
    index=False
)

print("Saved sample filtering summary.")
