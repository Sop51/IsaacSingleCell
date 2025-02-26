# this script is to subcluster within the inital cell clusters

# -------------------- hepatocyte subcluster ----------------------- #
# first pull out the hepatocyte
hep <- subset(x = zf, idents = "Hepatocyte")
# normalize the data
hep <- NormalizeData(hep)
# find variable features
hep <- FindVariableFeatures(hep, selection.method = "vst", nfeatures = 2000)
# scale the data on these top variable features
hep <- ScaleData(hep)
# run pca
hep <- RunPCA(hep, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(hep)
# construct a KNN graph - take dimensionality previously determined
hep <- FindNeighbors(hep, dims = 1:20)
# cluster the cells 
hep <- FindClusters(hep, resolution = 0.4)
# run umap - only displays LOCAL relationships
hep <- RunUMAP(hep, dims = 1:20)

# plot!!
DimPlot(hep, reduction = "umap", features = )

# marker selection
cluster0.hep <- FindMarkers(hep, ident.1 = 0)
cluster1.hep <- FindMarkers(hep, ident.1 = 1)
cluster2.hep <- FindMarkers(hep, ident.1 = 2)
cluster3.hep <- FindMarkers(hep, ident.1 = 3)
cluster4.hep <- FindMarkers(hep, ident.1 = 4)
cluster5.hep <- FindMarkers(hep, ident.1 = 5)
cluster6.hep <- FindMarkers(hep, ident.1 = 6)
cluster7.hep <- FindMarkers(hep, ident.1 = 7)

VlnPlot(hep, features = "anxa2a", group.by = "timepoint")
# ---------------- BECS ----------------------- #
# first pull out the hepatocyte
bec <- subset(x = zf, idents = "Biliary Epithelial Cell")
# normalize the data
bec <- NormalizeData(bec)
# find variable features
bec <- FindVariableFeatures(bec, selection.method = "vst", nfeatures = 2000)
# scale the data on these top variable features
bec <- ScaleData(bec)
# run pca
bec <- RunPCA(bec, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(hep)
# construct a KNN graph - take dimensionality previously determined
bec <- FindNeighbors(bec, dims = 1:20)
# cluster the cells 
bec <- FindClusters(bec, resolution = 0.4)
# run umap - only displays LOCAL relationships
bec <- RunUMAP(bec, dims = 1:20)

# plot!!
DimPlot(bec, reduction = "umap")

# marker selection
cluster0.bec <- FindMarkers(bec, ident.1 = 0)
cluster1.bec <- FindMarkers(bec, ident.1 = 1)
cluster2.bec <- FindMarkers(bec, ident.1 = 2)
cluster3.bec <- FindMarkers(bec, ident.1 = 3)
cluster4.bec <- FindMarkers(bec, ident.1 = 4)
cluster5.bec <- FindMarkers(bec, ident.1 = 5)
cluster6.bec <- FindMarkers(bec, ident.1 = 6)
cluster7.bec <- FindMarkers(bec, ident.1 = 7)

# -------------------- hepatocyte subcluster w integration ----------------------- #
# first pull out the hepatocyte
hep <- subset(x = zf, idents = "Hepatocyte")
# split by time point
timepoint_list <- SplitObject(hep, split.by = "timepoint")
# extract each seurat obj for each time point
untreated <- timepoint_list[["untreated"]]
mock <- timepoint_list[["mock"]]
dpa0 <- timepoint_list[["0dpa"]]
dpa1 <- timepoint_list[["1dpa"]]
dpa2 <- timepoint_list[["2dpa"]]
dpa3 <- timepoint_list[["3dpa"]]
dpa7 <- timepoint_list[["7dpa"]]

# split by sample for each time point
sample_list <- SplitObject(untreated, split.by = "orig.ident")

# normalize and find variable genes for each sample
for (i in 1:length(sample_list)) {
  sample_list[[i]] <- NormalizeData(object = sample_list[[i]])
  sample_list[[i]] <- FindVariableFeatures(object = sample_list[[i]])
}

 # select integration features
features <- SelectIntegrationFeatures(object.list = sample_list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = sample_list,
                                  anchor.features = features,
                                  k.filter = NA)

# integrate the data
integrated <- IntegrateData(anchorset = anchors, k.weight = 30)

# scale the data
integrated <- ScaleData(object = integrated)
integrated <- RunPCA(object = integrated)
integrated <- RunUMAP(object = integrated, dims = 1:50)

DimPlot(integrated, reduction = 'umap', group.by = 'orig.ident')
DimPlot(dpa7, reduction = 'umap', group.by = 'orig.ident')
