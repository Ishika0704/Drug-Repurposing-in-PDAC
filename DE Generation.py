# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 22:57:28 2025

@author: Ishika
"""
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# ----------------------------
# 1) Load Expression Matrix
# ----------------------------

expr = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA Dataset.txt", sep="\t", index_col=0)

# Ensure first row is not sample names accidentally
expr.columns = expr.columns.astype(str)

print("Expression matrix loaded with shape:", expr.shape)

# ----------------------------
# 2) Load Metadata
# ----------------------------

meta = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\metadata TCGA.csv")

# Make sure metadata samples match expression columns
meta = meta[meta["sample"].isin(expr.columns)]

# Reorder expression matrix according to metadata order
expr = expr[meta["sample"]]

print("Metadata loaded. Tumor/normal counts:")
print(meta["condition"].value_counts())

# ----------------------------
# 3) Build Design Matrix
# ----------------------------

meta["condition_binary"] = meta["condition"].map({"normal": 0, "tumor": 1})

design = sm.add_constant(meta["condition_binary"])

# ----------------------------
# 4) Run Differential Expression (OLS regression per gene)
# ----------------------------

genes = expr.index
logFC = []
pvalues = []

for gene in genes:
    y = expr.loc[gene].values
    model = sm.OLS(y, design).fit()
    coef = model.params["condition_binary"]   # tumor vs normal
    p = model.pvalues["condition_binary"]
    logFC.append(coef)
    pvalues.append(p)

# ----------------------------
# 5) Multiple Testing Correction (FDR)
# ----------------------------

fdr = multipletests(pvalues, method='fdr_bh')[1]

# ----------------------------
# 6) Create DE results table
# ----------------------------

results = pd.DataFrame({
    "gene": genes,
    "logFC": logFC,
    "pvalue": pvalues,
    "FDR": fdr
})

results = results.sort_values("FDR")

# Save results
results.to_csv("DE_results.csv", index=False)
print("DE_results.csv generated.")

# ----------------------------
# 7) Extract UP and DOWN genes (CMap-ready)
# ----------------------------

up = results[(results["logFC"] > 1) & (results["FDR"] < 0.05)]["gene"]
down = results[(results["logFC"] < -1) & (results["FDR"] < 0.05)]["gene"]

up.to_csv("up_genes.txt", index=False, header=False)
down.to_csv("down_genes.txt", index=False, header=False)

print("up_genes.txt and down_genes.txt created successfully!")
print(f"UP genes: {len(up)}")
print(f"DOWN genes: {len(down)}")

