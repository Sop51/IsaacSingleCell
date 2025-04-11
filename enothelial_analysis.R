library(patchwork)

# subset the initial seurat object
endo_apln <- subset(x = zf, idents = c("Endothelial Cell", "apln+"))
# break into timepoints
dpa0 <- subset(x = endo_apln, subset = timepoint == "0dpa")
dpa1 <- subset(x = endo_apln, subset = timepoint == "1dpa")
dpa2 <- subset(x = endo_apln, subset = timepoint == "2dpa")
dpa3 <- subset(x = endo_apln, subset = timepoint == "3dpa")
dpa7 <- subset(x = endo_apln, subset = timepoint == "7dpa")

mock <- subset(x = endo_apln, subset = timepoint == "untreated")
DimPlot(mock, group.by = 'cell.type.12.long') + ggtitle("untreated")
mock 

#look at them individually
# 0dpa first -------------
dpa0dimplot <- DimPlot(dpa0, group.by = 'cell.type.12.long') + ggtitle("0dpa")

# normalize the data
dpa0 <- NormalizeData(dpa0)
# use RNA assat for normalization and selecting variable features
DefaultAssay(dpa0) <- "RNA"
# find variable features
dpa0 <- FindVariableFeatures(dpa0, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(dpa0) <- "integrated"
# scale the data on these top variable features
dpa0 <- ScaleData(dpa0)
# run pca
dpa0 <- RunPCA(dpa0, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(dpa0)
# construct a KNN graph - take dimensionality previously determined
dpa0 <- FindNeighbors(dpa0, dims = 1:20)
# cluster the cells 
dpa0 <- FindClusters(dpa0, resolution = 0.4)
# run umap - only displays LOCAL relationships
dpa0 <- RunUMAP(dpa0, dims = 1:20)

# 1dpa ----------------
dpa1dimplot <- DimPlot(dpa1, group.by = 'cell.type.12.long') + ggtitle("1dpa")

# normalize the data
dpa1 <- NormalizeData(dpa1)
# use RNA assat for normalization and selecting variable features
DefaultAssay(dpa1) <- "RNA"
# find variable features
dpa1 <- FindVariableFeatures(dpa1, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(dpa1) <- "integrated"
# scale the data on these top variable features
dpa1 <- ScaleData(dpa1)
# run pca
dpa1 <- RunPCA(dpa1, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(dpa1)
# construct a KNN graph - take dimensionality previously determined
dpa1 <- FindNeighbors(dpa1, dims = 1:20)
# cluster the cells 
dpa1 <- FindClusters(dpa1, resolution = 0.4)
# run umap - only displays LOCAL relationships
dpa1 <- RunUMAP(dpa1, dims = 1:20)

DimPlot(dpa1) + ggtitle("1dpa")
cluster0.1dpa <- FindMarkers(dpa1, ident.1 = 0, only.pos = TRUE)
cluster1.1dpa <- FindMarkers(dpa1, ident.1 = 1, only.pos = TRUE)

# 2dpa ----------------
dpa2dimplot <- DimPlot(dpa2, group.by = 'cell.type.12.long') + ggtitle("2dpa")

# normalize the data
dpa2 <- NormalizeData(dpa2)
# use RNA assat for normalization and selecting variable features
DefaultAssay(dpa2) <- "RNA"
# find variable features
dpa2 <- FindVariableFeatures(dpa2, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(dpa2) <- "integrated"
# scale the data on these top variable features
dpa2 <- ScaleData(dpa2)
# run pca
dpa2 <- RunPCA(dpa2, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(dpa2)
# construct a KNN graph - take dimensionality previously determined
dpa2 <- FindNeighbors(dpa2, dims = 1:20)
# cluster the cells 
dpa2 <- FindClusters(dpa2, resolution = 0.4)
# run umap - only displays LOCAL relationships
dpa2 <- RunUMAP(dpa2, dims = 1:20)

# 3dpa ----------------
dpa3dimplot <- DimPlot(dpa3, group.by = 'cell.type.12.long') + ggtitle("3dpa")

# normalize the data
dpa3 <- NormalizeData(dpa3)
# use RNA assat for normalization and selecting variable features
DefaultAssay(dpa3) <- "RNA"
# find variable features
dpa3 <- FindVariableFeatures(dpa3, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(dpa3) <- "integrated"
# scale the data on these top variable features
dpa3 <- ScaleData(dpa3)
# run pca
dpa3 <- RunPCA(dpa3, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(dpa3)
# construct a KNN graph - take dimensionality previously determined
dpa3 <- FindNeighbors(dpa3, dims = 1:20)
# cluster the cells 
dpa3 <- FindClusters(dpa3, resolution = 0.4)
# run umap - only displays LOCAL relationships
dpa3 <- RunUMAP(dpa3, dims = 1:20)

# 7dpa ----------------
dpa7dimplot <- DimPlot(dpa7, group.by = 'cell.type.12.long') + ggtitle("7dpa")

# normalize the data
dpa7 <- NormalizeData(dpa7)
# use RNA assat for normalization and selecting variable features
DefaultAssay(dpa7) <- "RNA"
# find variable features
dpa7 <- FindVariableFeatures(dpa7, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(dpa7) <- "integrated"
# scale the data on these top variable features
dpa7 <- ScaleData(dpa7)
# run pca
dpa7 <- RunPCA(dpa7, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(dpa7)
# construct a KNN graph - take dimensionality previously determined
dpa7 <- FindNeighbors(dpa7, dims = 1:20)
# cluster the cells 
dpa7 <- FindClusters(dpa7, resolution = 0.4)
# run umap - only displays LOCAL relationships
dpa7 <- RunUMAP(dpa7, dims = 1:20)

# Arrange in 2 plots per column
dpa0dimplot + dpa1dimplot + dpa2dimplot + dpa3dimplot + dpa7dimplot + 
  plot_layout(ncol = 2)
