# Pathway enrischment analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db", force = TRUE)
BiocManager::install("clusterProfiler", force = TRUE)
# 1. Install clusterProfiler
BiocManager::install("clusterProfiler", force = TRUE)

# 2. Load both required packages
library(clusterProfiler)
library(org.Hs.eg.db)
# 1. Get top genes from Cluster 5 
levels(bc.seurat.obj)
markers_cluster5 <- FindMarkers(bc.seurat.obj, ident.1 = "Luminal")

genes_cluster5 <- head(markers_cluster5$gene, 50)
genes_cluster5 <- head(rownames(markers_cluster5), 300)
entrez_genes5 <- bitr(genes_cluster5, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)

nrow(entrez_genes5) 
# got mapped 289 genes out of 300
# GO Biological Process enrichment
ego5 <- enrichGO(gene          = entrez_genes5$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)

barplot(ego5, showCategory = 10, title = "GO Enrichment - Cluster 2 (Luminal)")

# KEGG pathway enrichment (optional)
ekegg5 <- enrichKEGG(gene         = entrez_genes5$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)

barplot(ekegg5, showCategory = 10, title = "KEGG Pathways - Cluster 5 (Luminal)")
#Now gonna repeat that code given above for each cluster
# Now I would simple change the ident. 1 by the name of the level (cluster)
# and it will automatically calculate the GO for that
