import pandas as pd
df = pd.read_csv(
    "GTEx_pancreas_geneSymbols_log2TPM.csv",
    index_col=0
)
print("Loaded shape:", df.shape)
print("Total genes:", df.shape[0])

dup_mask = df.index.duplicated(keep=False)
n_dup_rows = dup_mask.sum()
n_dup_genes = df.index[dup_mask].nunique()
print("\nDUPLICATE CHECK")
print("Rows involved in duplicates:", n_dup_rows)
print("Unique duplicated gene symbols:", n_dup_genes)

df.loc[dup_mask].head(10)
df_nodup = (
    df
    .groupby(df.index)
    .mean()
)
print("\nAFTER COLLAPSING DUPLICATES")
print("Shape:", df_nodup.shape)
print("Remaining duplicated genes:",
      df_nodup.index.duplicated().sum())

print("\nSANITY CHECK")
print("Min expression:", df_nodup.to_numpy().min())
print("Max expression:", df_nodup.to_numpy().max())

df_nodup.to_csv("GTEx_pancreas_final_clean.csv")

# Identify genes with zero expression across all samples
zero_mask = (df_nodup == 0).all(axis=1)
n_zero_genes = zero_mask.sum()
total_genes = df_nodup.shape[0]
print("\nSTEP 8: ZERO-EXPRESSION GENE CHECK")
print("Total genes before removal:", total_genes)
print("Genes with zero expression in all samples:", n_zero_genes)
print("Percentage:",
      round(100 * n_zero_genes / total_genes, 2), "%")

df_final = df_nodup.loc[~zero_mask]
print("\nSTEP 9: AFTER ZERO-EXPRESSION REMOVAL")
print("Genes remaining:", df_final.shape[0])
print("Any all-zero genes left?",
      ((df_final == 0).all(axis=1)).any())

print("\nSTEP 10: SANITY CHECK")
print("Min expression:", df_final.to_numpy().min())
print("Max expression:", df_final.to_numpy().max())

df_final.to_csv("GTEx_pancreas_final_clean_noZero.csv")
print("\nSTEP 11: FILE SAVED SUCCESSFULLY")
print("Output file: GTEx_pancreas_final_clean_noZero.csv")
print("Final shape (genes x samples):", df_final.shape)







