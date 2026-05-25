import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# STEP 0: Load data (no assumptions)
# -----------------------------
df = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_usable.csv",
    sep=","
)

print("Full dataframe shape:", df.shape)

# -----------------------------
# STEP 1: Identify numeric vs non-numeric columns
# -----------------------------
numeric_df = df.select_dtypes(include=[np.number])
non_numeric_cols = df.columns.difference(numeric_df.columns)

print("\nNumeric (expression) matrix shape:", numeric_df.shape)
print("\nNon-numeric (annotation) columns:")
print(list(non_numeric_cols))

# -----------------------------
# STEP 2: Raw count test (must be integers)
# -----------------------------
is_integer = np.all(
    np.equal(
        numeric_df.values,
        numeric_df.values.astype(int)
    )
)

print("\nAre all expression values integers (raw counts)?", is_integer)

# -----------------------------
# STEP 3: Log-transformation test
# -----------------------------
min_value = numeric_df.min().min()
max_value = numeric_df.max().max()

print("\nMinimum expression value:", min_value)
print("Maximum expression value:", max_value)

if max_value < 20:
    print("Data MAY be log-transformed")
else:
    print("Data is NOT log-transformed")

# -----------------------------
# STEP 4: TPM constant-sum test
# -----------------------------
col_sums = numeric_df.sum(axis=0)

print("\nColumn-wise sum statistics:")
print(col_sums.describe())

# Check if sums are approximately constant (TPM ~ 1e6)
is_tpm = np.allclose(col_sums.mean(), 1e6, rtol=0.1)

print("\nDo column sums approximate 1e6 (TPM-like)?", is_tpm)

# -----------------------------
# STEP 5: Distribution check (visual confirmation)
# -----------------------------
plt.figure()
numeric_df.iloc[:, 0].hist(bins=100)
plt.xlabel("Expression value")
plt.ylabel("Number of genes")
plt.title("Expression distribution (Sample 1)")
plt.show()

# -----------------------------
# STEP 6: Final computational verdict
# -----------------------------
print("\nFINAL COMPUTATIONAL VERDICT:")

if is_integer:
    print("- Data represents raw RNA-seq counts")
elif is_tpm:
    print("- Data represents TPM-normalized expression values")
else:
    print("- Data represents continuous, length-normalized expression values")
    print("  (FPKM-like or similar; not raw counts, not log-transformed)")
