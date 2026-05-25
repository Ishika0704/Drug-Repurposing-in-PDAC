# ======================================
# 1. IMPORT LIBRARIES
# ======================================
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ======================================
# 2. LOAD DATA
# ======================================
tumor1 = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA_primary_tumor_only.csv", index_col=0)
tumor2 = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GSE205154_pancreas_primary_expression.csv", index_col=0)
normal = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GTEx_pancreas_final_clean_noZero.csv", index_col=0)

# ======================================
# 3. KEEP COMMON GENES
# ======================================
common_genes = tumor1.index.intersection(tumor2.index)
common_genes = common_genes.intersection(normal.index)

tumor1 = tumor1.loc[common_genes]
tumor2 = tumor2.loc[common_genes]
normal = normal.loc[common_genes]

# ======================================
# 4. MERGE TUMOR DATASETS
# ======================================
tumor = pd.concat([tumor1, tumor2], axis=1)

print("Tumor samples:", tumor.shape[1])
print("Normal samples:", normal.shape[1])

# ======================================
# 5. OPTIONAL FILTER (zero in >50%)
# ======================================
combined = pd.concat([tumor, normal], axis=1)
zero_counts = (combined == 0).sum(axis=1)
threshold = 0.5 * combined.shape[1]
keep = zero_counts <= threshold

tumor = tumor.loc[keep]
normal = normal.loc[keep]

print("Genes after filtering:", tumor.shape[0])

# ======================================
# 6. RUN TESTS
# ======================================
results = []

for gene in tumor.index:

    tumor_vals = tumor.loc[gene].values
    normal_vals = normal.loc[gene].values

    if np.var(tumor_vals) == 0 and np.var(normal_vals) == 0:
        continue

    # Welch (primary)
    t_welch, p_welch = ttest_ind(
        tumor_vals, normal_vals, equal_var=False
    )

    # Student (secondary)
    t_student, p_student = ttest_ind(
        tumor_vals, normal_vals, equal_var=True
    )

    # Mann-Whitney
    try:
        u_stat, p_mw = mannwhitneyu(
            tumor_vals, normal_vals, alternative='two-sided'
        )
    except:
        p_mw = np.nan

    logFC = np.mean(tumor_vals) - np.mean(normal_vals)

    results.append([
        gene,
        logFC,
        p_welch,
        p_student,
        p_mw
    ])

results_df = pd.DataFrame(
    results,
    columns=[
        "Gene",
        "logFC",
        "p_welch",
        "p_student",
        "p_mannwhitney"
    ]
)

# ======================================
# 7. FDR CORRECTION
# ======================================
results_df["adj_p_welch"] = multipletests(
    results_df["p_welch"], method="fdr_bh"
)[1]

results_df["adj_p_student"] = multipletests(
    results_df["p_student"], method="fdr_bh"
)[1]

results_df["adj_p_mannwhitney"] = multipletests(
    results_df["p_mannwhitney"], method="fdr_bh"
)[1]

# ======================================
# 8. SORT & SAVE
# ======================================
results_df = results_df.sort_values("adj_p_welch")

results_df.to_csv("Combined_Tumor_vs_Normal_All_Tests.csv", index=False)

print("Done. Results saved.")