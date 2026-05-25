# =========================================================
# Prepare iLINCS input files from TCGA–GTEx contrast
# =========================================================

import pandas as pd

# ---------------------------------------------------------
# 1. Load TCGA–GTEx contrast file
# ---------------------------------------------------------

contrast_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\TCGA_GTEx_signature_p_and_adjP.csv"

contrast_df = pd.read_csv(
    contrast_path,
    index_col=0   # gene names as index
)

print("Total genes in contrast file:", contrast_df.shape[0])
print("Columns:", contrast_df.columns.tolist())

# ---------------------------------------------------------
# 2. Filter significant genes
#    (You can adjust thresholds if sir asks)
# ---------------------------------------------------------

LOGFC_CUTOFF = 1
ADJ_P_CUTOFF = 0.05

sig_genes = contrast_df[
    (contrast_df["Adjusted p-value"] < ADJ_P_CUTOFF) &
    (contrast_df["LogFC"].abs() > LOGFC_CUTOFF)
]

print("Significant genes after filtering:", sig_genes.shape[0])

# ---------------------------------------------------------
# 3. Create UP-regulated gene list (TCGA > GTEx)
# ---------------------------------------------------------

UP_N = 150   # iLINCS recommended: 50–300

up_genes = (
    sig_genes[sig_genes["LogFC"] > LOGFC_CUTOFF]
    .sort_values("LogFC", ascending=False)
    .head(UP_N)
)

print("UP genes selected:", up_genes.shape[0])

# ---------------------------------------------------------
# 4. Create DOWN-regulated gene list (TCGA < GTEx)
# ---------------------------------------------------------

DOWN_N = 150

down_genes = (
    sig_genes[sig_genes["LogFC"] < -LOGFC_CUTOFF]
    .sort_values("LogFC", ascending=True)
    .head(DOWN_N)
)

print("DOWN genes selected:", down_genes.shape[0])

# ---------------------------------------------------------
# 5. Save iLINCS-compatible files
#    (One gene symbol per line, NO header)
# ---------------------------------------------------------

up_output = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\iLINCS\PDAC_UP_genes.txt"
down_output = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\iLINCS\PDAC_DOWN_genes.txt"

up_genes.index.to_series().to_csv(
    up_output,
    index=False,
    header=False
)

down_genes.index.to_series().to_csv(
    down_output,
    index=False,
    header=False
)

print("\nFiles successfully generated:")
print("UP genes file  :", up_output)
print("DOWN genes file:", down_output)

# ---------------------------------------------------------
# 6. Optional: Save ranked tables for record-keeping
# ---------------------------------------------------------

up_genes.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\iLINCS\PDAC_UP_genes_full.csv"
)

down_genes.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\iLINCS\PDAC_DOWN_genes_full.csv"
)

print("\nRanked gene tables also saved for documentation.")
