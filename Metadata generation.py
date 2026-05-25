import pandas as pd
import re

# 🔹 Replace this path with your file path
file_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data 2021.txt"

# Load your expression matrix (TAB-separated)
expr = pd.read_csv(file_path, sep="\t", index_col=0)

# Get the sample names (columns)
samples = expr.columns

# Function to classify tumor vs normal using TCGA barcode suffix
def classify(sample_id):
    if re.search(r'-11$', sample_id):   # -11 = solid tissue normal
        return 'normal'
    else:
        return 'tumor'                  # everything else = tumor

# Build metadata table
metadata = pd.DataFrame({
    "sample": samples,
    "condition": [classify(s) for s in samples]
})

# Save metadata file in the same folder as your script
metadata.to_csv("metadata.csv", index=False)

print("Metadata file created successfully!")