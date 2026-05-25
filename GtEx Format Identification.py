import pandas as pd
file_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154.csv"
df_gtex = pd.read_csv(file_path, index_col=0)
print("GTEx pancreas-only file loaded successfully")

#Step 2
print("\n GSe205154 PANCREAS DATASET OVERVIEW")
print("-" * 40)
print(f"Number of genes   : {df_gtex.shape[0]}")
print(f"Number of samples : {df_gtex.shape[1]}")
print("\n First 5 genes:")
print(df_gtex.index[:5].tolist())
print("\n First 5 samples:")
print(df_gtex.columns[:5].tolist())

#Step 3
df_gtex_numeric = df_gtex.apply(pd.to_numeric, errors="coerce")
print("\n EXPRESSION VALUE CHECK (NUMERIC ONLY)")
value_summary = df_gtex_numeric.iloc[:1000, :5].describe().loc[['min', 'max']]
print(value_summary)
max_val = df_gtex_numeric.to_numpy().max()
if max_val > 1e5:
    expr_type = "Raw counts"
elif max_val > 100:
    expr_type = "TPM / FPKM"
else:
    expr_type = "TPM or log-normalized"
print(f"\n Likely expression type: {expr_type}")

#Step 4
print("\n Missing values:", df_gtex.isna().sum().sum())

#Step 5
print("\n🔁 Duplicated genes:", df_gtex.index.duplicated().sum())

#Step 6
def gtex_pancreas_summary(csv_file):
    df = pd.read_csv(csv_file, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")

    summary = {
        "Genes": df.shape[0],
        "Normal pancreas samples": df.shape[1],
        "Max expression value": df.to_numpy().max(),
        "Missing values": df.isna().sum().sum()
    }
    return summary


