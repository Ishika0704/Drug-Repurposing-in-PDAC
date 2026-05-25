import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
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
tcga = tcga.apply(pd.to_numeric, errors="coerce")
gse  = gse.apply(pd.to_numeric, errors="coerce")
gtex = gtex.apply(pd.to_numeric, errors="coerce")
common_genes = sorted(
    set(tcga.index) & set(gse.index) & set(gtex.index)
)

tcga = tcga.loc[common_genes]
gse  = gse.loc[common_genes]
gtex = gtex.loc[common_genes]
tumour = pd.concat([tcga, gse], axis=1)
normal = gtex

expr = pd.concat([tumour, normal], axis=1)

group = np.array(
    [1] * tumour.shape[1] + [0] * normal.shape[1]
)
results = []

for gene in common_genes:
    y = expr.loc[gene].values

    valid = ~np.isnan(y)
    y = y[valid]
    g = group[valid]

    # Design matrix: intercept + group
    X = sm.add_constant(g)

    model = sm.OLS(y, X)
    fit = model.fit()

    logFC = fit.params[1]          # group coefficient
    t_stat = fit.tvalues[1]
    p_val = fit.pvalues[1]

    tumour_mean = y[g == 1].mean()
    normal_mean = y[g == 0].mean()

    results.append([
        gene,
        tumour_mean,
        normal_mean,
        logFC,
        t_stat,
        p_val
    ])
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
results_df["adj_p_value"] = multipletests(
    results_df["p_value"],
    method="fdr_bh"
)[1]
results_df = results_df.sort_values("adj_p_value")

output_file = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA_GSE_vs_GTEx_linearModel_DE.csv"
results_df.to_csv(output_file, index=False)

print("Linear-model DE analysis complete.")
print("Results saved to:", output_file)
