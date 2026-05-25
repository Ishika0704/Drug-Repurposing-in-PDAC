import pandas as pd
import numpy as np

# --------------------------------
# STEP 1: Load the already formatted TPM matrix
# --------------------------------
expr_numeric = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_usable.csv",
    index_col=0
)

print("Loaded expression matrix shape:", expr_numeric.shape)

# --------------------------------
# STEP 2: Identify zero-expression genes
# --------------------------------
zero_genes = expr_numeric[(expr_numeric == 0).all(axis=1)]
nonzero_genes = expr_numeric[(expr_numeric > 0).any(axis=1)]

print("Total genes:", expr_numeric.shape[0])
print("Zero-expression genes removed:", zero_genes.shape[0])
print("Genes retained:", nonzero_genes.shape[0])

# --------------------------------
# STEP 3: Save removed genes separately
# --------------------------------
zero_genes.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_zero_expression_genes.csv"
)

# --------------------------------
# STEP 4: Log2(TPM + 1) transform remaining genes
# --------------------------------
expr_log2 = np.log2(nonzero_genes + 1)

# Sanity checks
print("Min after log transform:", expr_log2.min().min())
print("Max after log transform:", expr_log2.max().max())
print("Any NaNs?", expr_log2.isna().any().any())
print("Any infinite values?", np.isinf(expr_log2.values).any())

# --------------------------------
# STEP 5: Save final processed matrix
# --------------------------------
expr_log2.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_log2TPM_noZero_TCGA_like.csv"
)

print("Processing complete. Files saved successfully.")
