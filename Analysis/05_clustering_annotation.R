# 7. Clustering ------------
bc.seurat.obj <- FindNeighbors(bc.seurat.obj, dims = 1:15)
# understanding resolution -- 
bc.seurat.obj <- FindClusters(bc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(bc.seurat.obj@meta.data)
DimPlot(bc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
# setting identity of clusters
Idents(bc.seurat.obj)
Idents(bc.seurat.obj) <- "RNA_snn_res.0.1"
#setting up identity of bc object to resolution 0.1
#0.1 ident have 6 clusters, so now bc object will have 6 clusters of cells
#so based on that further clustering is done
Idents(bc.seurat.obj)
#UMAP 
bc.seurat.obj <- RunUMAP(bc.seurat.obj, dims = 1:15)
DimPlot(bc.seurat.obj, reduction = "umap")