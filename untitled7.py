import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# -----------------------------
# STEP 1: Load datasets
# -----------------------------
tcga = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA_primary_tumor_only.csv",
    index_col=0
)

gse = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GSE205154_pancreas_primary_expression.csv",
    index_col=0
)

gtex = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GTEx_pancreas_final_clean_noZero.csv",
    index_col=0
)

print("TCGA shape:", tcga.shape)
print("GSE shape :", gse.shape)
print("GTEx shape:", gtex.shape)

print("\nFirst 5 genes TCGA:", tcga.index[:5])
print("First 5 genes GSE :", gse.index[:5])
print("First 5 genes GTEx:", gtex.index[:5])

# -----------------------------
# STEP 2: Ensure numeric data
# -----------------------------
tcga = tcga.apply(pd.to_numeric, errors="coerce")
gse  = gse.apply(pd.to_numeric, errors="coerce")
gtex = gtex.apply(pd.to_numeric, errors="coerce")

# -----------------------------
# STEP 3: Identify common genes
# -----------------------------
genes_tcga = set(tcga.index)
genes_gse  = set(gse.index)
genes_gtex = set(gtex.index)

common_genes = genes_tcga & genes_gse & genes_gtex
common_genes = sorted(list(common_genes))  # ✅ CRITICAL FIX

print("Genes in TCGA:", len(genes_tcga))
print("Genes in GSE :", len(genes_gse))
print("Genes in GTEx:", len(genes_gtex))
print("Common genes across all 3:", len(common_genes))

# -----------------------------
# STEP 4: Track non-common genes
# -----------------------------
tcga_only = genes_tcga - set(common_genes)
gse_only  = genes_gse  - set(common_genes)
gtex_only = genes_gtex - set(common_genes)

non_common_df = pd.concat([
    pd.DataFrame({"Gene": list(tcga_only), "Source": "TCGA"}),
    pd.DataFrame({"Gene": list(gse_only),  "Source": "GSE"}),
    pd.DataFrame({"Gene": list(gtex_only), "Source": "GTEx"})
])

print("TCGA-only genes:", len(tcga_only))
print("GSE-only genes :", len(gse_only))
print("GTEx-only genes:", len(gtex_only))
print("Total non-common genes:", non_common_df.shape[0])

# -----------------------------
# STEP 5: Subset common genes
# -----------------------------
tcga_c = tcga.loc[common_genes]
gse_c  = gse.loc[common_genes]
gtex_c = gtex.loc[common_genes]

# Alignment checks
assert tcga_c.shape[0] == gse_c.shape[0] == gtex_c.shape[0]
assert all(tcga_c.index == gse_c.index)
assert all(tcga_c.index == gtex_c.index)

print("Gene alignment across datasets: PASSED")

# -----------------------------
# STEP 6: Combine tumour datasets
# -----------------------------
tumour = pd.concat([tcga_c, gse_c], axis=1)
normal = gtex_c

print("Tumour matrix shape:", tumour.shape)
print("Normal matrix shape:", normal.shape)
print("Tumour samples:", tumour.shape[1])
print("Normal samples:", normal.shape[1])

# -----------------------------
# STEP 7: Compute logFC
# -----------------------------
tumour_mean = tumour.mean(axis=1)
normal_mean = normal.mean(axis=1)

eps = 1e-6
logFC = np.log2((tumour_mean + eps) / (normal_mean + eps))

print("logFC summary:")
print(logFC.describe())
print("Genes with infinite logFC:", np.isinf(logFC).sum())

# -----------------------------
# STEP 8: Statistical testing
# -----------------------------
pvals = []

for gene in tumour.index:
    t_vals = tumour.loc[gene].values
    n_vals = normal.loc[gene].values
    _, p = ttest_ind(t_vals, n_vals, equal_var=False, nan_policy="omit")
    pvals.append(p)

pvals = np.array(pvals)

print("NaN p-values:", np.isnan(pvals).sum())
print("Min p-value:", np.nanmin(pvals))
print("Max p-value:", np.nanmax(pvals))

# -----------------------------
# STEP 9: FDR correction
# -----------------------------
adj_pvals = multipletests(pvals, method="fdr_bh")[1]

print("Adjusted p-value range:",
      adj_pvals.min(), "to", adj_pvals.max())

# -----------------------------
# STEP 10: Final results table
# -----------------------------
results = pd.DataFrame({
    "Gene": tumour.index,
    "Tumour_mean": tumour_mean.values,
    "Normal_mean": normal_mean.values,
    "logFC": logFC.values,
    "p_value": pvals,
    "adj_p_value": adj_pvals
}).sort_values("adj_p_value")

print(results.head())
print("Total genes tested:", results.shape[0])

# -----------------------------
# STEP 11: Save to Excel
# -----------------------------
with pd.ExcelWriter(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA_GSE_GTEx_DE_analysis.xlsx",
) as writer:
    
    results.to_excel(writer, sheet_name="Common_Genes_DEGs", index=False)
    non_common_df.to_excel(writer, sheet_name="Non_Common_Genes", index=False)

# -----------------------------
# STEP 12: DEG filtering
# -----------------------------
degs = results[
    (results["adj_p_value"] < 0.05) &
    (abs(results["logFC"]) > 1)
]

print("Significant DEGs:", degs.shape[0])
print("Upregulated:", (degs["logFC"] > 1).sum())
print("Downregulated:", (degs["logFC"] < -1).sum())
