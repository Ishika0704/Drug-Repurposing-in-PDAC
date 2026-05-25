import pandas as pd

# ---------------------------------------------------------
# 1. Load TCGA expression matrix
# ---------------------------------------------------------

tcga = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA Dataset.txt",
    sep="\t",
    index_col=0
)

print("TCGA data loaded")
print("Shape (genes x samples):", tcga.shape)

# ---------------------------------------------------------
# 2. Select first 15 samples (columns)
# ---------------------------------------------------------

tcga_15 = tcga.iloc[:, :15]

print("Subset with first 15 samples:")
print(tcga_15.shape)

# ---------------------------------------------------------
# 3. Calculate mean expression per gene
# ---------------------------------------------------------

tcga_mean_15 = tcga_15.mean(axis=1)

# ---------------------------------------------------------
# 4. Create summary DataFrame
# ---------------------------------------------------------

tcga_mean_df = pd.DataFrame({
    "GENE_ID": tcga_mean_15.index,
    "Avg_TCGA_expression_15_samples": tcga_mean_15.values
})

# ---------------------------------------------------------
# 5. Save small summary CSV
# ---------------------------------------------------------

output_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA_mean_expression_15_samples.csv"

tcga_mean_df.to_csv(
    output_path,
    index=False
)

print("Summary file saved at:")
print(output_path)
print("Rows:", tcga_mean_df.shape[0])
