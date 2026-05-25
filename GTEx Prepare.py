import pandas as pd
import os

# ---------------- PATHS ----------------
EXPR_FILE = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct\GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct"
META_FILE = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
OUT_FILE  = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\GTEx_Normal_Pancreas_TPM.csv"

# ---------------- CHECK FILES ----------------
print("--- Step 1: Checking Files ---")
if not os.path.exists(EXPR_FILE):
    raise FileNotFoundError(f"Expression file not found: {EXPR_FILE}")
if not os.path.exists(META_FILE):
    raise FileNotFoundError(f"Metadata file not found: {META_FILE}")
print("Files found. Proceeding...")

# ---------------- LOAD METADATA ----------------
print("\n--- Step 2: Filtering Metadata ---")
# Reading metadata to find which Sample IDs (SAMPID) belong to Pancreas
meta = pd.read_csv(META_FILE, sep="\t")
pancreas_meta = meta[meta["SMTSD"].str.contains("Pancreas", na=False)]
pancreas_samples = set(pancreas_meta["SAMPID"])

print(f"Total Pancreas samples identified: {len(pancreas_samples)}")

# ---------------- IDENTIFY COLUMNS ----------------
print("\n--- Step 3: Mapping Columns ---")
# We read only the header (skiprows=2 for GCT format) to get column names
header_df = pd.read_csv(EXPR_FILE, sep="\t", skiprows=2, nrows=0)
all_cols = header_df.columns.tolist()

# We want the ID and Name columns (usually index 0 and 1) plus any Pancreas IDs
# Note: GTEx IDs in the GCT file sometimes use '.' instead of '-' 
# so we check for both or partial matches if necessary.
cols_to_keep = [all_cols[0], all_cols[1]] + [c for c in all_cols if c in pancreas_samples]

print(f"Extracting {len(cols_to_keep)} columns (ID/Description + Pancreas Samples)")

# ---------------- CHUNKED EXTRACTION ----------------
print("\n--- Step 4: Extracting Data (Chunked) ---")

chunk_size = 5000  # Larger chunk size for better performance
written = False

# usecols is the secret to speed: it ignores non-pancreas data during the read
expr_iter = pd.read_csv(
    EXPR_FILE,
    sep="\t",
    skiprows=2,
    chunksize=chunk_size,
    usecols=cols_to_keep
)

try:
    for i, chunk in enumerate(expr_iter):
        # Append to CSV: 'w' for the first chunk (write), 'a' for others (append)
        chunk.to_csv(
            OUT_FILE,
            mode="w" if not written else "a",
            index=False,
            header=not written
        )
        written = True
        processed_rows = (i + 1) * chunk_size
        print(f"   > Processed ~{processed_rows} genes...")

    print(f"\n SUCCESS: File saved to {OUT_FILE}")

except Exception as e:
    print(f"\n ERROR during processing: {e}")