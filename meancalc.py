import pandas as pd

tcga = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA_primary_tumor_only.csv",
    index_col=0
)

gtex = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GTEx_pancreas_final_clean_noZero.csv",
    index_col=0
)

gse = pd.read_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Excel Files after codes\GSE205154_pancreas_primary_expression.csv",
    index_col=0
)

print("TCGA shape:", tcga.shape)
print("GTEx shape:", gtex.shape)
print("GSE shape:", gse.shape)
tcga_mean = tcga.mean(axis=1)
gtex_mean = gtex.mean(axis=1)
gse_mean = gse.mean(axis=1)
tcga_mean_df = tcga_mean.to_frame(name="TCGA_mean")
gtex_mean_df = gtex_mean.to_frame(name="GTEx_mean")
gse_mean_df  = gse_mean.to_frame(name="GSE_mean")
#tcga_mean_df.to_csv(
   # r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\TCGA_gene_mean_expression.csv"
#)

gtex_mean_df.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GTEx_gene_mean_expression.csv"
)

#gse_mean_df.to_csv(
 #   r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\\GSE205154_gene_mean_expression.csv"
#)
