# this script is to subcluster within the inital cell clusters
library(ggpubr)
library(ggplot2)
library(ggsignif)

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
DimPlot(hep, reduction = "umap", group.by = 'timepoint')
UMAPPlot(hep)
FeaturePlot(hep, features = 'anxa2a', group_by = 'timepoint')

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
DimPlot(bec, reduction = "umap", group.by = "timepoint")

# marker selection
cluster0.bec <- FindMarkers(bec, ident.1 = 0)
cluster1.bec <- FindMarkers(bec, ident.1 = 1)
cluster2.bec <- FindMarkers(bec, ident.1 = 2)
cluster3.bec <- FindMarkers(bec, ident.1 = 3)
cluster4.bec <- FindMarkers(bec, ident.1 = 4)
cluster5.bec <- FindMarkers(bec, ident.1 = 5)
cluster6.bec <- FindMarkers(bec, ident.1 = 6)
cluster7.bec <- FindMarkers(bec, ident.1 = 7)

# ---------------- neutrophils ----------------------- #
# first pull out the neutrophils
neutrophil <- subset(x = zf, idents = "Neutrophil")
# normalize the data
neutrophil <- NormalizeData(neutrophil)
# find variable features
neutrophil <- FindVariableFeatures(neutrophil, selection.method = "vst", nfeatures = 2000)
# scale the data on these top variable features
neutrophil <- ScaleData(neutrophil)
# run pca
neutrophil <- RunPCA(neutrophil, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(neutrophil)
# construct a KNN graph - take dimensionality previously determined
neutrophil <- FindNeighbors(neutrophil, dims = 1:20)
# cluster the cells 
neutrophil <- FindClusters(neutrophil, resolution = 0.4)
# run umap - only displays LOCAL relationships
neutrophil <- RunUMAP(neutrophil, dims = 1:20)

# plot!!
DimPlot(neutrophil, reduction = "umap", group.by = "timepoint")
VlnPlot(neutrophil, features = 'myd88', group.by = "timepoint")

# define the feature (gene) of interest
gene <- 'cxcr4b'

# subset the neutrophil data for cxcr4b expression greater than 0
expressing_cells <- neutrophil@meta.data[neutrophil@assays$RNA@data[gene, ] > 0, ]

# function to compute the mean difference between two groups
mean_diff <- function(group1, group2) {
  mean(group1) - mean(group2)
}

# set the number of permutations
n_permutations <- 1000

# observed data
timepoints <- unique(expressing_cells$timepoint)
comparison_pairs <- list()
observed_diff <- numeric(0)
counter <- 1
p_values <- numeric(0)

# loop through timepoints to compute observed differences and p-values
for (i in 1:(length(timepoints)-1)) {
  for (j in (i+1):length(timepoints)) {
    
    group1 <- expressing_cells$timepoint == timepoints[i]
    group2 <- expressing_cells$timepoint == timepoints[j]
    
    # extract values for the two groups using a feature of interest
    data1 <- neutrophil@assays$RNA@data[gene, ][group1]
    data2 <- neutrophil@assays$RNA@data[gene, ][group2]
    
    # compute the observed mean difference
    observed_diff[counter] <- mean_diff(data1, data2)
    
    # perform the permutation test
    perm_diffs <- numeric(n_permutations)
    combined_data <- c(data1, data2)
    
    for (perm in 1:n_permutations) {
      permuted_labels <- sample(combined_data)
      permuted_group1 <- permuted_labels[1:length(data1)]
      permuted_group2 <- permuted_labels[(length(data1)+1):length(combined_data)]
      
      perm_diffs[perm] <- mean_diff(permuted_group1, permuted_group2)
    }
    
    # calculate p-value - the proportion of permuted differences greater than or equal to the observed difference
    p_values[counter] <- mean(abs(perm_diffs) >= abs(observed_diff[counter]))
    
    # store the comparison pairs
    comparison_pairs[[counter]] <- paste(timepoints[i], "vs", timepoints[j])
    
    counter <- counter + 1
  }
}

# create a data frame for plotting
results_df <- data.frame(
  comparison = unlist(comparison_pairs),
  observed_diff = observed_diff,
  p_value = p_values
)

# box plot with permutation test p-values
ggplot(expressing_cells, aes(x = timepoint, y = neutrophil@assays$RNA@data[gene, ][neutrophil@assays$RNA@data[gene, ] > 0], fill = timepoint)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, size = 0.7, color = "black", width = 0.5) +  
  geom_jitter(width = 0.2, alpha = 0.7, size = 2, color = "gray30") +  
  scale_fill_manual(values = cbPalette) +  
  theme_minimal(base_size = 14) +  
  labs(title = paste("Expression of", gene, "in Neutrophils"), 
       x = "Timepoint", 
       y = "Expression") +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.title.x = element_text(size = 16, face = "bold"), 
    axis.title.y = element_text(size = 16, face = "bold"),  
    axis.text = element_text(size = 14),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
    legend.position = "none"  # Remove legend
  ) +
  # add significance bars using geom_signif
  geom_signif(
    comparisons = list(
      c("untreated", "0dpa"),
      c("untreated", "1dpa"),
      c("untreated", "2dpa"),
      c("untreated", "3dpa")
    ),
    annotations = c(
      paste("p =", format(results_df$p_value[results_df$comparison == "untreated vs 0dpa"][1], digits = 3)),
      paste("p =", format(results_df$p_value[results_df$comparison == "untreated vs 1dpa"][1], digits = 3)),
      paste("p =", format(results_df$p_value[results_df$comparison == "untreated vs 2dpa"][1], digits = 3)),
      paste("p =", format(results_df$p_value[results_df$comparison == "untreated vs 3dpa"][1], digits = 3))
    ),
    y_position = c(
      max(neutrophil@assays$RNA@data['cxcr4b', ]) + 0.1,
      max(neutrophil@assays$RNA@data['cxcr4b', ]) + 0.3,
      max(neutrophil@assays$RNA@data['cxcr4b', ]) + 0.5,
      max(neutrophil@assays$RNA@data['cxcr4b', ]) + 0.7
    ),
    tip_length = 0.03
  )

# marker selection
cluster0.neutrophil <- FindMarkers(neutrophil, ident.1 = 0)
cluster1.neutrophil <- FindMarkers(neutrophil, ident.1 = 1)
cluster2.neutrophil <- FindMarkers(neutrophil, ident.1 = 2)
cluster3.neutrophil <- FindMarkers(neutrophil, ident.1 = 3)
cluster4.neutrophil <- FindMarkers(neutrophil, ident.1 = 4)
cluster5.neutrophil <- FindMarkers(neutrophil, ident.1 = 5)
cluster6.neutrophil <- FindMarkers(neutrophil, ident.1 = 6)
cluster7.neutrophil <- FindMarkers(neutrophil, ident.1 = 7)

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
