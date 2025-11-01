#now we have to do QC, 1. %of mt.genes per cell-simply proxy for cell death, so higer mt genes = highly damaged cell
View(bc.seurat.obj@meta.data)
#% of MT genes
bc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(
  bc.seurat.obj,
  pattern = "^MT-"
)
head(bc.seurat.obj@meta.data)
View(bc.seurat.obj@meta.data)
#Visualize the metadata, give us idea about the data distribution, and set up QC filters
VlnPlot(bc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)
FeatureScatter(bc.seurat.obj, feature1 = "nFeature_RNA", feature2 = 'nCount_RNA' )
geom_smooth(method = "lm")
#-----filtering cells-----
#check quartiles to see what cut off should be set for QC
quantile(bc.seurat.obj$nFeature_RNA, probs = c(0.95, 0.99))
#major genes lie between 7100 - 8250
quantile(bc.seurat.obj$nCount_RNA, probs = c(0.95, 0.99))
#major UMIS b/w 60k - 91k
#we'll set QC accordingly 
bc.seurat.obj <- subset(
  bc.seurat.obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 8200 &
    nCount_RNA < 92000 &
    percent.mt < 10
)