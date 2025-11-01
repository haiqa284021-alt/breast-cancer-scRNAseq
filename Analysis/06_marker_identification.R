#Looking for well known bc markers in sample cell population-- Panel-based (literature driven)
# marker panels
epi <- c("EPCAM","KRT8","KRT18","KRT19","MUC1","GATA3")
immune <- c("PTPRC","CD3D","CD3E","CD8A","CD4","CD68","CD14","MS4A1")
stroma <- c("THY1","FAP","PDGFRA","PDGFRB","ACTA2","COL1A1","RGS5","PECAM1","VWF")
other <- c("MKI67","TOP2A","VIM","FN1","SNAI1")

# summary dotplot (all panels)
DotPlot(bc.seurat.obj, features = c(epi, immune, stroma, other)) + RotatedAxis()

# per-gene UMAP locations
FeaturePlot(bc.seurat.obj, features = c("EPCAM","PTPRC","THY1","MKI67"), min.cutoff="q10", max.cutoff="q90")

# violin to see expression distributions by cluster
VlnPlot(bc.seurat.obj, features = c("EPCAM","PTPRC","THY1","MKI67"))

# find markers for a cluster (e.g., cluster 0)
# we can find marker genes per cluster and do it manaully one by one
FindMarkers(bc.seurat.obj, ident.1 = 5, min.pct = 0.25, logfc.threshold = 0.25) -> markers_cluster0
head(markers_cluster0)
# run the above program by changing ident to the number of cluster to find respective marker genes
# Cluster 0 is Epithelial
# Cluster 1 likely represents stromal or perivascular-like cells
# cluster 2 is epithelial/stromal transition or luminal-like epithelial cells
# cluster 3 is Tumor-Associated Macrophages (TAMs) or Monocyte-derived macrophages.
# cluster 4 is stromal/fibroblast (cancer-associated fibroblast, CAF)
# cluster 5 is Neuroendocrine-like / Neuronal-like tumor cells

# UMAP annotation based on Cluster types
# 1. Define cluster annotations based on your identified cell types
new.cluster.ids <- c(
  "Epithelial",    # Cluster 0
  "Stromal",       # Cluster 1
  "Luminal",       # Cluster 2
  "Macrophage",    # Cluster 3
  "Fibroblast",    # Cluster 4
  "Neuroendocrine" # Cluster 5
)

# 2. Assign names matching your cluster numbers
names(new.cluster.ids) <- levels(bc.seurat.obj)

# 3. Rename cluster identities
bc.seurat.obj <- RenameIdents(bc.seurat.obj, new.cluster.ids)

# 4. Plot annotated UMAP
DimPlot(bc.seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.8) +
  ggtitle("Annotated UMAP: Breast Cancer Single-Cell Atlas") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Annotated_UMAP.png", width = 6, height = 5, dpi = 300)

# Now find all markers for all 6 clusters
all_markers <- FindAllMarkers(bc.seurat.obj, 
                              only.pos = TRUE, # only genes upregulated in each cluster
                              min.pct = 0.25,  # gene expressed in at least 25% cells
                              logfc.threshold = 0.25 # minimum log fold-change
)
# now let's find top 5 markers per cluster in a single file
library(dplyr)
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
table(top5_markers$cluster)

#Visualize top markers
DotPlot(bc.seurat.obj, features = unique(top5_markers$gene)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 5 Marker Genes per Cluster")
# heatmap
DoHeatmap(bc.seurat.obj, features = unique(top5_markers$gene)) +
  ggtitle("Expression Heatmap of Top 5 Markers per Cluster")