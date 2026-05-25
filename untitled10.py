# ======================================================
# STEP 0: Imports
# ======================================================
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# ======================================================
# STEP 1: Load expression matrices
# ======================================================
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

print("Shapes:")
print("TCGA:", tcga.shape)
print("GSE :", gse.shape)
print("GTEx:", gtex.shape)

# ======================================================
# STEP 2: Ensure numeric values
# ======================================================
tcga = tcga.apply(pd.to_numeric, errors="coerce")
gse  = gse.apply(pd.to_numeric, errors="coerce")
gtex = gtex.apply(pd.to_numeric, errors="coerce")

# ======================================================
# STEP 3: Find common genes
# ======================================================
genes_tcga = set(tcga.index)
genes_gse  = set(gse.index)
genes_gtex = set(gtex.index)

common_genes = sorted(list(genes_tcga & genes_gse & genes_gtex))

print("Common genes:", len(common_genes))

tcga = tcga.loc[common_genes]
gse  = gse.loc[common_genes]
gtex = gtex.loc[common_genes]

# ======================================================
# STEP 4: Combine tumour samples
# ======================================================
tumour = pd.concat([tcga, gse], axis=1)
normal = gtex

print("Tumour samples:", tumour.shape[1])
print("Normal samples:", normal.shape[1])

# ======================================================
# STEP 5: Compute gene-wise statistics
# ======================================================
results = []

for gene in common_genes:
    tumour_vals = tumour.loc[gene].dropna().values
    normal_vals = normal.loc[gene].dropna().values

    # Welch's t-test (UNPAIRED)
    t_stat, p_val = ttest_ind(
        tumour_vals,
        normal_vals,
        equal_var=False
    )

    tumour_mean = tumour_vals.mean()
    normal_mean = normal_vals.mean()

    logFC = np.log2((tumour_mean + 1e-6) / (normal_mean + 1e-6))

    results.append([
        gene,
        tumour_mean,
        normal_mean,
        logFC,
        t_stat,
        p_val
    ])

# ======================================================
# STEP 6: Create results DataFrame
# ======================================================
results_df = pd.DataFrame(
    results,
    columns=[
        "Gene",
        "Tumour_mean",
        "Normal_mean",
        "logFC",
        "t_statistic",
        "p_value"
    ]
)

# ======================================================
# STEP 7: Multiple testing correction
# ======================================================
results_df["adj_p_value"] = multipletests(
    results_df["p_value"],
    method="fdr_bh"
)[1]

# ======================================================
# STEP 8: Sort and save
# ======================================================
results_df = results_df.sort_values("adj_p_value")

output_file = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA_GSE_vs_GTEx_geneWise_ttest.csv"
results_df.to_csv(output_file, index=False)

print("Analysis complete.")
print("Results saved to:", output_file)
