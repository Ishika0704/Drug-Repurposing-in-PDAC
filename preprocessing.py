import pandas as pd
import os
# ---------------- PATHS ----------------
EXPR_FILE = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct\GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct"
META_FILE = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
OUT_FILE  = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\GTEx_Normal_Pancreas_TPM.csv"

# ---------------- CHECK FILES ----------------
print("Checking files...")
print("Expression exists:", os.path.exists(EXPR_FILE))
print("Metadata exists:", os.path.exists(META_FILE))

if not os.path.exists(EXPR_FILE):
    raise FileNotFoundError("Expression file not found")
if not os.path.exists(META_FILE):
    raise FileNotFoundError("Metadata file not found")

# ---------------- LOAD METADATA ----------------
print("\nLoading metadata...")
meta = pd.read_csv(META_FILE, sep="\t")

pancreas_meta = meta[meta["SMTSD"].str.contains("Pancreas", na=False)]
pancreas_samples = set(pancreas_meta["SAMPID"])

print("Pancreas samples found:", len(pancreas_samples))

# ---------------- CHUNKED EXPRESSION READ ----------------
print("\nStarting chunked expression read...")

chunk_size = 300
written = False

expr_iter = pd.read_csv(
    EXPR_FILE,
    sep="\t",
    skiprows=2,        # GCT header lines
    chunksize=chunk_size
)

for chunk in expr_iter:
    gene_col = chunk.columns[0]
    keep_cols = [gene_col] + [c for c in chunk.columns if c in pancreas_samples]

    chunk = chunk[keep_cols]

    chunk.to_csv(
        OUT_FILE,
        mode="w" if not written else "a",
        index=False,
        header=not written
    )

    written = True
    print("Processed chunk")

print("\n✅ DONE: GTEx pancreas TPM file saved")
