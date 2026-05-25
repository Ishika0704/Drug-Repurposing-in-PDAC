import pandas as pd

# ---------- File paths ----------
input_txt = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_Gene_Level_TPM_Estimates.txt"
output_csv = routput_csv = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154.csv"

# ---------- Read TCGA .txt file ----------
# sep='\t' is critical for TCGA files
df = pd.read_csv(
    input_txt,
    sep="\t",
    header=0,
    index_col=0
)

# ---------- Sanity checks ----------
print("Shape of data (genes x samples):", df.shape)
print("First 5 genes:")
print(df.head())

# ---------- Save as CSV ----------
df.to_csv(output_csv)

print("Conversion successful! CSV saved at:")
print(output_csv)
