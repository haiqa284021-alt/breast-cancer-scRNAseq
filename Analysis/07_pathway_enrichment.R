# Pathway enrischment analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db", force = TRUE)
BiocManager::install("clusterProfiler", force = TRUE)
# 1. Install clusterProfiler
BiocManager::install("clusterProfiler", force = TRUE)

# 2. Load both required packages
 
barplot(ekegg5, showCategory = 10, title = "KEGG Pathways - Cluster 5 (Luminal)")
#Now gonna repeat that code given above for each cluster
# Now I would simple change the ident. 1 by the name of the level (cluster)
# and it will automatically calculate the GO for that
