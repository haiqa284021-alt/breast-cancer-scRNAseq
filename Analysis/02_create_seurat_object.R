library(Seurat)
library(SeuratDisk)
library(tidyverse)
matrix_file <- "data/Breast_Cancer_3p_filtered_feature_bc_matrix.h5"
hdf5_obj <- Read10X_h5(filename = matrix_file)
bc.seurat.obj <- CreateSeuratObject(counts = hdf5_obj, project = "Breast_Cancer")
saveRDS(bc.seurat.obj, file = "data/bc_raw.rds")
