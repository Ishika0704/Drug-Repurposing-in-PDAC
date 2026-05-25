import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

# -----------------------------
# Load TCGA tumor expression data
# -----------------------------
tcga = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA Dataset.txt",
    sep="\t",
    index_col=0
)

# -----------------------------
# Load GTEx normal pancreas data
# -----------------------------
gtex = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GTEx_noZero_genes.csv",
    index_col=0
)

# -----------------------------
# Find common genes
# -----------------------------
common_genes = tcga.index.intersection(gtex.index)
tcga_c = tcga.loc[common_genes]
gtex_c = gtex.loc[common_genes]

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
# Linear fold change
# -----------------------------
epsilon = 1e-6
fold_change = (tcga_signature + epsilon) / (gtex_signature + epsilon)

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
        equal_var=True,      # Unpaired Student's t-test
        nan_policy="omit"
    )
    p_values.append(p)

p_values = pd.Series(p_values, index=logFC.index)

# -----------------------------
# Final disease signature table
# -----------------------------
contrast_df = pd.DataFrame({
    "TCGA-GTEx": logFC,
    "Fold change": fold_change.loc[logFC.index],
    "LogFC": logFC,
    "p value": p_values
})

# -----------------------------
# Save output
# -----------------------------
contrast_df.to_csv("TCGA GTEx signature.csv")

print("File saved: TCGA_GTEx_disease_signature.csv")
print("Columns:", contrast_df.columns.tolist())
print("Total genes:", contrast_df.shape[0])
print("\nTop 5 genes:")
print(contrast_df.head())
