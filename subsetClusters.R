# first pull out the hepatocyte
hep <- subset(x = zf, idents = "Hepatocyte")

hep_list <- SplitObject(hep, split.by = "orig.ident")
hep_list <- lapply(hep_list, SCTransform)

# select integration features
features <- SelectIntegrationFeatures(hep_list, nfeatures = 3000)

# prepare for integration (normalization)
hep_list <- PrepSCTIntegration(hep_list, anchor.features = features)

# find anchors for integration
anchors <- FindIntegrationAnchors(hep_list, normalization.method = "SCT", anchor.features = features)

# integrate the data (the actual merging step)
hep_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 46)

# perform PCA on the integrated hepatocyte data
hep_combined <- RunPCA(hep_combined)

# perform UMAP on the PCA reduction
hep_combined <- RunUMAP(hep_combined, reduction = "pca", dims = 1:10)

# find neighbors (based on PCA dimensions)
hep_combined <- FindNeighbors(hep_combined, dims = 1:10)

# find clusters
hep_combined <- FindClusters(hep_combined, dims = 1:10)

