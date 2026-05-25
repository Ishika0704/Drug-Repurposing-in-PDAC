import pandas as pd

# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------

gct_path = r"C:\Users\Ishika\Downloads\Connectivity Drugs.gct"
csv_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\iLINCS\connectivity_results.csv"

# ---------------------------------------------------------
# Read GCT v1.3 file
# ---------------------------------------------------------

# Step 1: Skip first 2 metadata lines
df = pd.read_csv(
    gct_path,
    sep="\t",
    skiprows=2
)

# Step 2: Remove the 'desc' row (first data row)
df = df[df.iloc[:, 0] != "desc"]

# Optional: reset index
df.reset_index(drop=True, inplace=True)

print("GCT file parsed successfully")
print("Shape:", df.shape)
print("First 5 columns:", df.columns[:5].tolist())

# ---------------------------------------------------------
# Save as CSV
# ---------------------------------------------------------

df.to_csv(
    csv_path,
    index=False
)

print("CSV file saved at:")
print(csv_path)
