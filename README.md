#  Breast Cancer scRNA-seq Analysis

This repository contains the full R pipeline for single-cell RNA sequencing (scRNA-seq) analysis of a breast cancer dataset using **Seurat** and **clusterProfiler**. The goal is to identify distinct cell populations and interpret their biological functions via GO and KEGG pathway enrichment.

---

## Project Structure

breast-cancer-scRNAseq/
│
├── README.md
├── sessionInfo.txt
│
├── analysis/
│ ├── 01_install_packages.R
│ ├── 02_create_seurat_object.R
│ ├── 03_QC_filtering.R
│ ├── 04_normalization_PCA.R
│ ├── 05_clustering_annotation.R
│ ├── 06_marker_identification.R
│ ├── 07_pathway_enrichment.R
│
├── figures/
│ ├── Annotated_UMAP.png
│ ├── GO_Fibroblast.png
│ ├── KEGG_Fibroblast.png
│ ├── GO_Luminal.png
│ ├── KEGG_Luminal.png
│ ├── GO_Macrophage.png
│ ├── KEGG_Macrophage.png
│ ├── GO_Stromal.png
│ ├── KEGG_Stromal.png
│ ├── GO_Endothelial.png
│ ├── KEGG_Endothelial.png
│ ├── GO_Neuroendocrine.png
│ ├── KEGG_Neuroendocrine.png
│
├── data/
│ └── Breast_Cancer_3p_filtered_feature_bc_matrix.h5

---

##  Objective

To characterize breast cancer microenvironment cell populations through single-cell transcriptomics and interpret their biological significance via pathway enrichment.

---

## Dataset

- **Source:** 10x Genomics filtered feature-barcode matrix (`.h5`)
- **Input file:** `Breast_Cancer_3p_filtered_feature_bc_matrix.h5`
- **Platform:** 10x Genomics Chromium 3' single-cell RNA-seq

---

## Methods Summary

1. **Data Loading & QC:** Imported `.h5` matrix, filtered cells by mitochondrial content, UMI count, and feature count.  
2. **Normalization & PCA:** Log-normalization, variable feature selection, and PCA for dimensionality reduction.  
3. **Clustering & UMAP:** Graph-based clustering and visualization via UMAP.  
4. **Annotation:** Cluster identities assigned based on canonical breast cancer cell markers.  
5. **Marker Identification:** Differential gene expression for each cluster.  
6. **Pathway Enrichment:** GO and KEGG analyses using `clusterProfiler`.

---

## Cluster Identities

| Cluster | Identity        | Description |
|----------|-----------------|--------------|
| 0 | **Epithelial** | Core tumor epithelial population |
| 1 | **Stromal** | Stromal or perivascular-like cells |
| 2 | **Luminal** | Luminal-like epithelial cells |
| 3 | **Macrophage** | Tumor-associated macrophages |
| 4 | **Fibroblast** | Cancer-associated fibroblasts (CAFs) |
| 5 | **Neuroendocrine** | Neuronal-like tumor cells |

---

## Biological Interpretation (Pathway Highlights)

- **Fibroblast:** Enriched for ECM organization and collagen binding pathways (GO, KEGG).  
- **Macrophage:** Immune response and cytokine signaling enrichment.  
- **Stromal:** Vascular development and smooth muscle contraction signatures.  
- **Endothelial:** Angiogenesis and cell adhesion pathways.  
- **Neuroendocrine:** Synaptic signaling and neurotransmitter-related GO terms.  
- **Luminal:** Chromosome segregation and mitotic nuclear division processes (highly proliferative).

---

## Key Figures

All UMAPs, GO, and KEGG plots are available in the `figures/` folder.

---

## Dependencies

- R ≥ 4.2  
- Seurat  
- SeuratDisk  
- tidyverse  
- clusterProfiler  
- org.Hs.eg.db  
- GO.db  

Install automatically via:
```r
source("analysis/01_install_packages.R")
# Clone this repo
git clone https://github.com/<your-username>/breast-cancer-scRNAseq.git
cd breast-cancer-scRNAseq

# Open RStudio and run:
source("analysis/02_create_seurat_object.R")
source("analysis/03_QC_filtering.R")
source("analysis/04_normalization_PCA.R")
source("analysis/05_clustering_annotation.R")
source("analysis/06_marker_identification.R")
source("analysis/07_pathway_enrichment.R")
