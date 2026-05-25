import pandas as pd
import numpy as np
print("Libraries loaded successfully")

df = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\GTEx_Normal_Pancreas_TPM.csv", index_col=0)
print("\nSTEP 1: FILE LOAD CHECK")
print("Shape (genes x samples):", df.shape)
print("First 5 row names:", df.index[:5].tolist())
print("First 5 column names:", df.columns[:5].tolist())

df = df.drop(columns=["Description"], errors="ignore")
print("\nSTEP 2: DESCRIPTION COLUMN REMOVAL")
print("Columns containing 'Description':",
      [c for c in df.columns if "Description" in c])
print("Updated shape:", df.shape)

df = df.apply(pd.to_numeric, errors="coerce")
print("\nSTEP 3: NUMERIC CONVERSION CHECK")
print("Total NaN values:", df.isna().sum().sum())
print("Data types present:\n", df.dtypes.value_counts())

df.index = df.index.str.split(".").str[0]
print("\nSTEP 4: ENSEMBL VERSION REMOVAL")
print("Example gene IDs:", df.index[:5].tolist())
print("Any IDs still containing '.' ?",
      df.index.str.contains("\.").any())

dup_count = df.index.duplicated().sum()
print("\nSTEP 5: DUPLICATE GENE CHECK")
print("Duplicated gene IDs:", dup_count)

mapping = pd.read_csv(r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\ensembl_to_symbol.csv")
print("\nSTEP 6: MAPPING FILE CHECK")
print(mapping.head())
print("Mapping columns:", mapping.columns.tolist())
print("Unique Ensembl IDs:", mapping["gene_id"].nunique())

mapping_dict = dict(
    zip(mapping["gene_id"], mapping["gene_name"])
)
df["GeneSymbol"] = df.index.map(mapping_dict)
print("\nSTEP 7: MAPPING DIAGNOSTICS")
print("Mapped genes:", df["GeneSymbol"].notna().sum())
print("Unmapped genes:", df["GeneSymbol"].isna().sum())

df = df.dropna(subset=["GeneSymbol"])
print("\nSTEP 8: POST-MAPPING SHAPE")
print("Final shape:", df.shape)

df = df.set_index("GeneSymbol")
print("\nSTEP 9: INDEX CHECK")
print("Index example:", df.index[:5].tolist())
print("Any duplicated gene symbols?",
      df.index.duplicated().sum())

max_val = df.to_numpy().max()
min_val = df.to_numpy().min()
print("\nSTEP 10: EXPRESSION SCALE CHECK")
print("Min expression:", min_val)
print("Max expression:", max_val)
if max_val > 1e5:
    expr_type = "Raw counts"
elif max_val > 100:
    expr_type = "TPM (linear scale)"
else:
    expr_type = "Log-transformed"
print("Likely expression type:", expr_type)

df_log = np.log2(df + 1)
print("\nSTEP 11: LOG TRANSFORMATION CHECK")
print("New max:", df_log.to_numpy().max())
print("New min:", df_log.to_numpy().min())

df_log.to_csv("GTEx_pancreas_geneSymbols_log2TPM.csv")
print("\nSTEP 12: FILE SAVED")
print("Output file: GTEx_pancreas_geneSymbols_log2TPM.csv")




