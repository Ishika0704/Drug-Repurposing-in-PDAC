# =========================================================
# Unsupervised clustering of TCGA + GTEx (log-normalized)
# =========================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore

# ---------------------------------------------------------
# 1. Load log-normalized expression matrices
#    Rows = genes, Columns = samples
# ---------------------------------------------------------

tcga_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Data\TCGA_Dataset.csv"
gtex_path = r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\GTEx_pancreas_final_clean_noZero.csv"

tcga = pd.read_csv(tcga_path, index_col=0)
gtex = pd.read_csv(gtex_path, index_col=0)

print("TCGA shape:", tcga.shape)
print("GTEx shape:", gtex.shape)

# ---------------------------------------------------------
# 2. Keep only common genes
# ---------------------------------------------------------

common_genes = tcga.index.intersection(gtex.index)

tcga = tcga.loc[common_genes]
gtex = gtex.loc[common_genes]

print("Common genes:", len(common_genes))

# ---------------------------------------------------------
# 3. Merge TCGA + GTEx
# ---------------------------------------------------------

combined = pd.concat([tcga, gtex], axis=1)

print("Combined matrix shape:", combined.shape)

# ---------------------------------------------------------
# 4. Z-score normalization per gene (SAFE VERSION)
# ---------------------------------------------------------

combined_z = pd.DataFrame(
    zscore(combined, axis=1, nan_policy="omit"),
    index=combined.index,
    columns=combined.columns
)

# Drop genes with zero variance (NaNs after z-score)
combined_z = combined_z.dropna()

print("After Z-score (genes retained):", combined_z.shape[0])


# ---------------------------------------------------------
# 5. Select top variable genes
# ---------------------------------------------------------

gene_variance = combined_z.var(axis=1)

TOP_N_GENES = 1000  # try 500 / 1000 / 2000
top_genes = gene_variance.sort_values(ascending=False).head(TOP_N_GENES).index

combined_var = combined_z.loc[top_genes]

print(f"Using top {TOP_N_GENES} variable genes")

# ---------------------------------------------------------
# 6. Create sample annotations (TCGA vs GTEx)
# ---------------------------------------------------------

sample_type = [
    "TCGA" if col.startswith("TCGA") else "GTEx"
    for col in combined_var.columns
]

col_colors = pd.Series(sample_type, index=combined_var.columns).map(
    {"TCGA": "#d62728", "GTEx": "#1f77b4"}
)

# ---------------------------------------------------------
# 7. Generate unsupervised clustering heatmap
# ---------------------------------------------------------

sns.set(style="white")

g = sns.clustermap(
    combined_var,
    method="ward",
    metric="euclidean",
    cmap="vlag",
    col_colors=col_colors,
    xticklabels=False,
    yticklabels=False,
    figsize=(12, 14)
)

plt.show()

# ---------------------------------------------------------
# 8. (Optional) Save figure
# ---------------------------------------------------------

g.savefig(
   r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\TCGA_GTEx_Clustermap.png",
   dpi=300,
   bbox_inches="tight"
)
# =========================================================
# 9. Extract clustering results (from the heatmap)
# =========================================================

from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

# Get the column linkage matrix (sample clustering)
col_linkage = g.dendrogram_col.linkage

# ---------------------------------------------------------
# 10. Define number of clusters (start with 2: tumor vs normal)
# ---------------------------------------------------------

N_CLUSTERS = 2

sample_clusters = fcluster(
    col_linkage,
    t=N_CLUSTERS,
    criterion="maxclust"
)

cluster_df = pd.DataFrame({
    "Sample": combined_var.columns,
    "Cluster": sample_clusters,
    "Type": [
        "TCGA" if col.startswith("TCGA") else "GTEx"
        for col in combined_var.columns
    ]
})

print("\nSample clustering summary:")
print(cluster_df.head())

# ---------------------------------------------------------
# 11. TCGA vs GTEx composition per cluster
# ---------------------------------------------------------

cluster_composition = pd.crosstab(
    cluster_df["Cluster"],
    cluster_df["Type"]
)

print("\nCluster composition (TCGA vs GTEx):")
print(cluster_composition)

# Percentage composition
cluster_percent = cluster_composition.div(
    cluster_composition.sum(axis=1),
    axis=0
) * 100

print("\nCluster composition (%):")
print(cluster_percent.round(2))

# ---------------------------------------------------------
# 12. Silhouette score (quantitative cluster quality)
# ---------------------------------------------------------

sil_score = silhouette_score(
    combined_var.T,
    sample_clusters,
    metric="euclidean"
)

print(f"\nSilhouette score (k={N_CLUSTERS}): {sil_score:.3f}")

# Interpretation helper
if sil_score > 0.5:
    print("→ Strong cluster separation")
elif sil_score > 0.25:
    print("→ Moderate but meaningful separation")
else:
    print("→ Weak separation (possible batch or heterogeneity effects)")

# ---------------------------------------------------------
# 13. PCA confirmation of clustering
# ---------------------------------------------------------

pca = PCA(n_components=2)
pca_coords = pca.fit_transform(combined_var.T)

pca_df = pd.DataFrame(
    pca_coords,
    columns=["PC1", "PC2"],
    index=combined_var.columns
)

pca_df["Cluster"] = sample_clusters
pca_df["Type"] = cluster_df["Type"].values

print("\nPCA variance explained:")
print(f"PC1: {pca.explained_variance_ratio_[0]*100:.2f}%")
print(f"PC2: {pca.explained_variance_ratio_[1]*100:.2f}%")

# ---------------------------------------------------------
# 14. Identify genes driving separation (top contributors)
# ---------------------------------------------------------

pc1_loadings = pd.Series(
    pca.components_[0],
    index=combined_var.index
)

top_pc1_genes = pc1_loadings.abs().sort_values(ascending=False).head(20)

print("\nTop genes contributing to PC1 separation:")
print(top_pc1_genes)

# ---------------------------------------------------------
# 15. Save analysis outputs (optional but recommended)
# ---------------------------------------------------------

cluster_df.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\sample_clusters.csv",
    index=False
)

top_pc1_genes.to_csv(
    r"C:\Users\Ishika\OneDrive\Desktop\Capstone Project\Results\top_PC1_genes.csv"
)
