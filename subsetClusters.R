# this script is to subcluster within the inital cell clusters
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(harmony)
library(Seurat)

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
FeaturePlot(hep, features = 'ezh2')

# marker selection
cluster0.hep <- FindMarkers(hep, ident.1 = 0)
cluster1.hep <- FindMarkers(hep, ident.1 = 1)
cluster2.hep <- FindMarkers(hep, ident.1 = 2)
cluster3.hep <- FindMarkers(hep, ident.1 = 3)
cluster4.hep <- FindMarkers(hep, ident.1 = 4)
cluster5.hep <- FindMarkers(hep, ident.1 = 5)
cluster6.hep <- FindMarkers(hep, ident.1 = 6)
cluster7.hep <- FindMarkers(hep, ident.1 = 7)

VlnPlot(hep, features = "ezh2", group.by = "timepoint")
# ---------------- BECS ----------------------- #
# first pull out the hepatocyte
bec <- subset(x = zf, idents = "Biliary Epithelial Cell")
# use RNA assat for normalization and selecting variable features
DefaultAssay(end) <- "RNA"
# normalize the data
bec <- NormalizeData(bec)
# find variable features
bec <- FindVariableFeatures(bec, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(mac) <- "integrated"
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

# marker selection
cluster0.neutrophil <- FindMarkers(neutrophil, ident.1 = 0)
cluster1.neutrophil <- FindMarkers(neutrophil, ident.1 = 1)
cluster2.neutrophil <- FindMarkers(neutrophil, ident.1 = 2)
cluster3.neutrophil <- FindMarkers(neutrophil, ident.1 = 3)
cluster4.neutrophil <- FindMarkers(neutrophil, ident.1 = 4)
cluster5.neutrophil <- FindMarkers(neutrophil, ident.1 = 5)
cluster6.neutrophil <- FindMarkers(neutrophil, ident.1 = 6)
cluster7.neutrophil <- FindMarkers(neutrophil, ident.1 = 7)

# ----------------- Endothelial Cell Subcluster ------------------ #
# first pull out the endothelial cells
end <- subset(x = zf, idents = "Endothelial Cell")
# use RNA assat for normalization and selecting variable features
DefaultAssay(end) <- "RNA"
# remove genes expressed in small num of cells to avoid simpleLoess warning
gene_counts <- rowSums(as.matrix(GetAssayData(end, slot = "counts")) > 0)
min_cells <- 0.03 * ncol(end)  # 3% of total cells
end <- subset(end, features = names(gene_counts[gene_counts > min_cells]))
# normalize the data
end <- NormalizeData(end)
# find variable features
end <- FindVariableFeatures(end, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(mac) <- "integrated"
# scale the data on these top variable features
end <- ScaleData(end)
# run pca
end <- RunPCA(end, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(end)
# construct a KNN graph - take dimensionality previously determined
end <- FindNeighbors(end, dims = 1:20)
# cluster the cells 
end <- FindClusters(end, resolution = 0.4)
# run umap - only displays LOCAL relationships
end <- RunUMAP(end, dims = 1:20)

# plot!!
x <- DimPlot(end, reduction = "umap")
y <- DimPlot(end, reduction = "umap", group.by = "timepoint")

x | y

# ------------------ Macrophage cell cluster ---------------------- #
# first pull out the macrophage cells
mac <- subset(x = zf, idents = "Macrophage")

# removing ribosomal and heatshock genes
counts <- GetAssayData(mac, assay = "RNA", slot = "counts")
genes_to_remove <- rownames(counts)[grepl("^rps|^hsp|^rpl", rownames(counts))]
counts <- counts[!(rownames(counts) %in% genes_to_remove), ]
mac <- subset(mac, features = rownames(counts))

# use RNA assat for normalization and selecting variable features
DefaultAssay(mac) <- "RNA"
mac <- NormalizeData(mac)
# find variable features
mac <- FindVariableFeatures(mac, selection.method = "vst", nfeatures = 2000)
existing_variable_features <- VariableFeatures(mac)
combined_variable_features <- unique(c(existing_variable_features, genes_of_interest))
VariableFeatures(mac) <- combined_variable_features
# change assay to integrated for the latter steps
DefaultAssay(mac) <- "integrated"
# scale the data on these top variable features
mac <- ScaleData(mac)
# run pca
mac <- RunPCA(mac, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(mac)
# construct a KNN graph - take dimensionality previously determined
mac <- FindNeighbors(mac, dims = 1:20)
# cluster the cells 
mac <- FindClusters(mac, resolution = 0.3)
# run umap - only displays LOCAL relationships
mac <- RunUMAP(mac, dims = 1:20)

# plot!!
x <- FeaturePlot(mac, features = "cd63")
y <- DimPlot(mac, reduction = "umap", group.by = "timepoint")

x | y


DefaultAssay(mac) <- "RNA"
# marker selection
cluster0.mac <- FindMarkers(mac, ident.1 = 0, only.pos = TRUE)
cluster1.mac <- FindMarkers(mac, ident.1 = 1, only.pos = TRUE)
cluster2.mac <- FindMarkers(mac, ident.1 = 2, only.pos = TRUE)
cluster3.mac <- FindMarkers(mac, ident.1 = 3, only.pos = TRUE)
cluster4.mac <- FindMarkers(mac, ident.1 = 4, only.pos = TRUE)
# ------------------------ Lymphocyte subcluster -------------------------- #
# first pull out the lymphocyte cells
lymph <- subset(x = zf, idents = "Lymphocyte")
# use RNA assat for normalization and selecting variable features
DefaultAssay(lymph) <- "RNA"
lymph <- NormalizeData(lymph)
# find variable features
lymph <- FindVariableFeatures(lymph, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(lymph) <- "integrated"
# scale the data on these top variable features
lymph <- ScaleData(lymph)
# run pca
lymph <- RunPCA(lymph, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(lymph)
# construct a KNN graph - take dimensionality previously determined
lymph <- FindNeighbors(lymph, dims = 1:15)
# cluster the cells 
lymph <- FindClusters(lymph, resolution = 0.4)
# run umap - only displays LOCAL relationships
lymph <- RunUMAP(lymph, dims = 1:20)

# plot!!
x <- DimPlot(lymph, reduction = "umap")
y <- DimPlot(lymph, reduction = "umap", group.by = "timepoint")

x | y

# ---------------------- apln+ cell cluster ------------------------------ #
# first pull out the apln+ cells
apln <- subset(x = zf, idents = "apln+")
# use RNA assat for normalization and selecting variable features
DefaultAssay(apln) <- "RNA"
# remove genes expressed in small num of cells to avoid simpleLoess warning
gene_counts <- rowSums(as.matrix(GetAssayData(apln, slot = "counts")) > 0)
min_cells <- 0.01 * ncol(apln)  # 1% of total cells
apln <- subset(apln, features = names(gene_counts[gene_counts > min_cells]))
apln <- NormalizeData(apln)
# find variable features
apln <- FindVariableFeatures(apln, selection.method = "vst", nfeatures = 2000)
# change assay to integrated for the latter steps
DefaultAssay(apln) <- "integrated"
# scale the data on these top variable features
apln <- ScaleData(apln)
# run pca
apln <- RunPCA(apln, npcs = 30)
# create the elbow plot to determine pcs to use going forward
ElbowPlot(apln)
# construct a KNN graph - take dimensionality previously determined
apln <- FindNeighbors(apln, dims = 1:18)
# cluster the cells 
apln <- FindClusters(apln, resolution = 0.2)
# run umap - only displays LOCAL relationships
apln <- RunUMAP(apln, dims = 1:18)

# plot!!
x <- DimPlot(apln, reduction = "umap")
y <- DimPlot(apln, reduction = "umap", group.by = "timepoint")

x | y

DefaultAssay(apln) <- "RNA"
zero <- FindMarkers(apln, ident.1 = 0, only.pos = TRUE)
one <- FindMarkers(apln, ident.1 = 1, only.pos = TRUE)

# Extract and order by p-value
top10_zero_pval <- head(zero, 10)
top10_one_pval <- head(one, 10)

# Extract gene names for the DotPlot
top10_genes_pval <- c(rownames(top10_zero_pval), rownames(top10_one_pval))

# Create a DotPlot of the top 10 genes
# Plot top 10 genes ordered by p-value
# Plot expression of the top 10 genes ordered by p-value
VlnPlot(apln, features = top10_genes_pval, pt.size = 0.1)
# ----------- code to produce a box plot for a cell type across time points ------------- #
celltype <- bec
gene <- 'spint2'

DefaultAssay(celltype) <- 'RNA'
# subset to only include cells expressing the gene
non_zero <- subset(celltype, subset = spint2 > 0)

# convert to a data frame for ggplot
data <- FetchData(non_zero, vars = c(gene, "timepoint"))

# create a box plot
ggplot(data, aes(x = timepoint, y = spint2, fill = timepoint)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black", linewidth = 0.3) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +               
  scale_fill_manual(values = c("#8ECAE6", 
                               "#d6d3cc", "#ffefd3", "#ffc49b", 
                               "#ffb236", "#ff9633","#ff9")) + 
  labs(title = paste("Biliary Epithelial Cell Temporal", gene, "Expression"), 
       x = "Timepoint", 
       y = "Normalized Expression Level") +
  theme_minimal(base_size = 18) +                                                 
  theme(
    legend.position = "none",                                                   
    axis.title = element_text(size = 14, face = "bold"),                         
    axis.text = element_text(size = 16, color = "black"),                        
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),           
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),         
    panel.grid.minor = element_blank(),                                          
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)    
  )

# ---------------------- plot with permutation testing -------------------- #
# define the feature (gene) of interest
gene <- 's1pr4'

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
      c("mock", "0dpa"),
      c("mock", "1dpa"),
      c("mock", "2dpa"),
      c("mock", "3dpa")
    ),
    annotations = c(
      paste("p =", format(results_df$p_value[results_df$comparison == "mock vs 0dpa"][1], digits = 3)),
      paste("p =", format(results_df$p_value[results_df$comparison == "mock vs 1dpa"][1], digits = 3)),
      paste("p =", format(results_df$p_value[results_df$comparison == "mock vs 2dpa"][1], digits = 3)),
      paste("p =", format(results_df$p_value[results_df$comparison == "mock vs 3dpa"][1], digits = 3))
    ),
    y_position = c(
      max(neutrophil@assays$RNA@data[gene, ]) + 0.1,
      max(neutrophil@assays$RNA@data[gene, ]) + 0.3,
      max(neutrophil@assays$RNA@data[gene, ]) + 0.5,
      max(neutrophil@assays$RNA@data[gene, ]) + 0.7
    ),
    tip_length = 0.03
  )