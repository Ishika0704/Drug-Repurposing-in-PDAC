import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# -----------------------------
# Load TCGA tumor expression data
# -----------------------------
tcga = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA_noZero_genes.csv",
    sep=",",
    index_col=0
)

# -----------------------------
# Load GTEx normal pancreas data
# -----------------------------
gtex = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GTEx_noZero_genes.csv",
    sep=",",
    index_col=0
)

# -----------------------------
# Find common genes
# -----------------------------
common_genes = tcga.index.intersection(gtex.index)
tcga_c = tcga.loc[common_genes]
gtex_c = gtex.loc[common_genes]

print("Common genes:", len(common_genes))

# -----------------------------
# Mean expression per gene
# -----------------------------
tcga_signature = tcga_c.mean(axis=1)
gtex_signature = gtex_c.mean(axis=1)

# -----------------------------
# LogFC (TCGA - GTEx)
# -----------------------------
logFC = tcga_signature - gtex_signature
logFC = logFC.sort_values(ascending=False)

# -----------------------------
# P-value calculation
# Unpaired Student’s t-test
# -----------------------------
p_values = []

for gene in logFC.index:
    tumor_vals = tcga_c.loc[gene].values
    normal_vals = gtex_c.loc[gene].values

    _, p = ttest_ind(
        tumor_vals,
        normal_vals,
        equal_var=True,      # Student's t-test
        nan_policy="omit"
    )
    p_values.append(p)

p_values = pd.Series(p_values, index=logFC.index, name="p-value")

# -----------------------------
# Adjust p-values (Benjamini–Hochberg FDR)
# -----------------------------
_, adj_pvals, _, _ = multipletests(
    p_values,
    method="fdr_bh"
)

adj_pvals = pd.Series(adj_pvals, index=logFC.index, name="Adjusted p-value")

# -----------------------------
# Final disease signature table
# -----------------------------
contrast_df = pd.DataFrame({
    "TCGA-GTEx": logFC,
    "LogFC": logFC,
    "p-value": p_values,
    "Adjusted p-value": adj_pvals
})

# -----------------------------
# Save output
# -----------------------------
output_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\TCGA_GTEx_signature_p_and_adjP.csv"
contrast_df.to_csv(output_path)

print("File saved:", output_path)
print("Columns:", contrast_df.columns.tolist())
print("Total genes:", contrast_df.shape[0])
print("\nTop 5 genes:")
print(contrast_df.head())
