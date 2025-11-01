#----Normalization Step------
bc.seurat.obj <- NormalizeData(bc.seurat.obj)
#to check which types of command being run on our seurat object
str(bc.seurat.obj)
#Finding highly vairable genes in dataset
bc.seurat.obj <- FindVariableFeatures(
  bc.seurat.obj,
  selection.method = "vst",
  nfeatures = 2000
)
#top 10 variable genes
top_10 <- head(VariableFeatures(bc.seurat.obj), 10)
top_10
#Plot the variable genes
plot1 <- VariableFeaturePlot(bc.seurat.obj)
LabelPoints(plot = plot1, points = top_10, repel = TRUE)
# 5. Scaling -------------
all.genes <- rownames(bc.seurat.obj)
bc.seurat.obj <- ScaleData(bc.seurat.obj, features = all.genes)

str(bc.seurat.obj)
#linear dimension reduction
bc.seurat.obj <- RunPCA(bc.seurat.obj, features = VariableFeatures(object = bc.seurat.obj))
#plot PCA results-- by changing dim number to 1-5 PCA number, you can visualize other PCA too
print(bc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(bc.seurat.obj, dims = 1, cells = 2310, balanced = TRUE)

VizDimLoadings(bc.seurat.obj, dims = 1:2, reduction = "pca")
# determine dimensionality of the data--selection of PCA for downstream analyses
#MUST select PCA having high STD, stop where elbow becomes plateau
ElbowPlot(bc.seurat.obj)
#selected till 15 PCs
