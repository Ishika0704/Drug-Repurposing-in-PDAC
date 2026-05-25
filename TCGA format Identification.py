import pandas as pd
file_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA Dataset.txt"
df = pd.read_csv(
    file_path,
    sep="\t",
    index_col=0
)
print("File loaded successfully")

print("\n DATASET OVERVIEW")
print("-" * 40)
n_genes, n_samples = df.shape
print(f"Number of genes   : {n_genes}")
print(f"Number of samples : {n_samples}")
print("\n First 5 genes:")
print(df.index[:5].tolist())
print("\n First 5 samples:")
print(df.columns[:5].tolist())

# Extract sample type code (01 = Tumor, 11 = Normal)
sample_type_codes = df.columns.str.split("-").str[3].str[:2]
sample_type_map = {
    "01": "Primary Tumor",
    "11": "Solid Tissue Normal"
}
sample_types = sample_type_codes.map(sample_type_map).fillna("Other")
print("\n SAMPLE TYPE DISTRIBUTION")
print(sample_types.value_counts())


print("\n EXPRESSION VALUE CHECK")
value_summary = df.iloc[:1000, :5].describe().loc[['min', 'max']]
print(value_summary)
max_val = df.values.max()
if max_val > 1e5:
    expr_type = "Raw read counts (HTSeq)"
elif max_val > 100:
    expr_type = "FPKM / FPKM-UQ"
else:
    expr_type = "TPM or log-normalized"
print(f"\n➡ Likely expression type: {expr_type}")

# Remove HTSeq technical rows
df_clean = df[~df.index.str.startswith("__")]
print("\n🧹 CLEANED DATASET")
print(f"Genes after cleanup: {df_clean.shape[0]}")

def tcga_summary(txt_file):
    df = pd.read_csv(txt_file, sep="\t", index_col=0)
    df = df[~df.index.str.startswith("__")]
    sample_types = df.columns.str.split("-").str[3].str[:2]
    sample_types = sample_types.map({"01": "Tumor", "11": "Normal"}).fillna("Other")
    summary = {
        "Genes": df.shape[0],
        "Samples": df.shape[1],
        "Tumor samples": (sample_types == "Tumor").sum(),
        "Normal samples": (sample_types == "Normal").sum(),
        "Max expression value": df.values.max()
    }
    return summary
print(tcga_summary(file_path))

