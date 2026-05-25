import pandas as pd
import numpy as np

df = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154.csv",
    sep=","
)
# Separate numeric TPM matrix
expr_numeric = df.select_dtypes(include=[np.number])

# Add gene symbol column
expr_numeric["gene"] = df["hgnc_symbol"]
expr_numeric = expr_numeric[
    expr_numeric["gene"].notna() & (expr_numeric["gene"] != "")
]
expr_numeric = (
    expr_numeric
    .groupby("gene", as_index=False)
    .mean()
)
expr_numeric = expr_numeric.set_index("gene")
expr_numeric = expr_numeric.sort_index()
expr_numeric.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\GSE205154_usable.csv"
)
print(expr_numeric.shape)
expr_numeric.head()

