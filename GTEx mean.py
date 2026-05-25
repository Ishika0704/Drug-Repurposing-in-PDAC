import pandas as pd

# ---------------------------------------------------------
# 1. Load GTEx expression matrix
# ---------------------------------------------------------

gtex = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GTEx_pancreas_final_clean_noZero.csv",
    index_col=0
)

print("GTEx data loaded")
print("Shape (genes x samples):", gtex.shape)

# ---------------------------------------------------------
# 2. Select first 15 samples (columns)
# ---------------------------------------------------------

gtex_15 = gtex.iloc[:, :15]

print("Subset with first 15 samples:")
print(gtex_15.shape)

# ---------------------------------------------------------
# 3. Calculate mean expression per gene
# ---------------------------------------------------------

gtex_mean_15 = gtex_15.mean(axis=1)

# ---------------------------------------------------------
# 4. Create summary DataFrame
# ---------------------------------------------------------

gtex_mean_df = pd.DataFrame({
    "GENE_ID": gtex_mean_15.index,
    "Avg_GTEx_expression_15_samples": gtex_mean_15.values
})

# ---------------------------------------------------------
# 5. Save small summary CSV
# ---------------------------------------------------------

output_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GTEx_mean_expression_15_samples.csv"

gtex_mean_df.to_csv(
    output_path,
    index=False
)

print("GTEx summary file saved at:")
print(output_path)
print("Rows:", gtex_mean_df.shape[0])
