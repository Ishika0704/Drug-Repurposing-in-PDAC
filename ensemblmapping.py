import pandas as pd

GTF_FILE = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\gencode.v26.chr_patch_hapl_scaff.annotation.gtf\gencode.v26.chr_patch_hapl_scaff.annotation.gtf"

gtf = pd.read_csv(
    GTF_FILE,
    sep="\t",
    comment="#",
    header=None,
    low_memory=False
)

gtf.columns = [
    "chr", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes"
]

genes = gtf[gtf["feature"] == "gene"].copy()

genes["gene_id"] = genes["attributes"].str.extract(r'gene_id "([^"]+)"')
genes["gene_name"] = genes["attributes"].str.extract(r'gene_name "([^"]+)"')

genes["gene_id"] = genes["gene_id"].str.replace(r"\.\d+", "", regex=True)

mapping = genes[["gene_id", "gene_name"]].drop_duplicates()

mapping.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\ensembl_to_symbol.csv",
    index=False
)

print("Saved ensembl_to_symbol.csv")
print(mapping.head())
