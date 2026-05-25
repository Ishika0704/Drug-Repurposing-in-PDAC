# Integrated Transcriptomic Drug Repurposing Pipeline for PDAC
This repository contains the computational workflow used for transcriptomics-guided drug repurposing in pancreatic ductal adenocarcinoma (PDAC). The pipeline integrates TCGA, GEO (GSE205154), and GTEx transcriptomic datasets to identify differentially expressed genes, perform connectivity mapping using iLINCS and SigCom LINCS, and validate candidate compounds through molecular docking and network analysis.
Associated manuscript:

Agarwal I, Ganesan S.
Integrated transcriptomic connectivity mapping reveals potential therapeutic candidates for pancreatic ductal adenocarcinoma.
Submitted to Frontiers in Genetics.

## Workflow Overview

1. Data acquisition and preprocessing
2. TCGA, GEO, and GTEx integration
3. Differential expression analysis
4. PDAC signature refinement
5. Connectivity mapping (iLINCS + SigCom LINCS)
6. Drug prioritization
7. Molecular docking validation
8. Gene-drug interaction network analysis
9. Visualization and figure generation

## Datasets Used

- TCGA-PAAD
- GTEx Pancreas
- GEO: GSE205154

## Software Requirements

Python 3.10+
R 4.3+

Python packages:
- pandas
- numpy
- scipy
- statsmodels
- matplotlib
- seaborn

R packages:
- ggplot2
- pheatmap
- limma
- edgeR
- ComplexHeatmap
- igraph
- ggraph

## Pipeline Execution Order
1. Run preprocessing scripts
2. Generate cleaned TCGA, GEO, and GTEx matrices
3. Run differential expression analysis
4. Generate DEG signature
5. Prepare iLINCS/SigCom input files
6. Run connectivity mapping
7. Generate visualizations and networks

## Notes
Some exploratory and developmental scripts used during optimization of the computational workflow are included in the archive/ directory for transparency and reproducibility.

## Acknowledgements
We acknowledge TCGA, GTEx, GEO, iLINCS, SigCom LINCS, STRING, DrugBank, and DGIdb for providing publicly accessible datasets and computational resources.

