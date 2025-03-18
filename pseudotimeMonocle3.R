library(monocle3)
library(Seurat)
library(dyno)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(dplyr)

# -------------------- getting the data set up ------------------------ #
DefaultAssay(ecm_subset) <- "RNA" # use RNA assay for normalization and selecting variable features

# filter for high quality cells
zf_filtered <- subset(zf, subset = nCount_RNA > 800 &
                    nFeature_RNA > 500)

# remove genes that are only expressed in a small number of cells
expr_mat <- GetAssayData(zf_filtered, layer = "counts")
genes_to_keep <- rowSums(expr_mat > 0) >= 10
zf_filtered <- subset(zf_filtered, features = rownames(expr_mat)[genes_to_keep])

# set cell type as the ident for sub-setting
Idents(zf_filtered) <- zf_filtered$cell.type.12.long

# filter for ECM producing cells
ecm_subset <- subset(zf_filtered, idents = c("Macrophage"), 
                     subset = timepoint %in% c("mock", "untreated") == FALSE)

# pre-processing of data
# normalize the data
ecm_subset <- NormalizeData(ecm_subset)

# find variable features
ecm_subset <- FindVariableFeatures(ecm_subset, selection.method = "vst", nfeatures = 2000)

# change assay to integrated for the latter steps
DefaultAssay(ecm_subset) <- "integrated"

# scale the data on these top variable features
ecm_subset <- ScaleData(ecm_subset)

# run pca
ecm_subset <- RunPCA(ecm_subset, npcs = 30)

# construct a KNN graph - take dimensionality previously determined
ecm_subset <- FindNeighbors(ecm_subset, dims = 1:30)

# cluster the cells 
ecm_subset <- FindClusters(ecm_subset, resolution = 0.2)

# run umap - only displays LOCAL relationships
ecm_subset <- RunUMAP(ecm_subset, dims = 1:30, n.neighbors = 50)

# plot
x <- DimPlot(ecm_subset, reduction = "umap", label = T)
y <- DimPlot(ecm_subset, reduction = "umap", group.by = 'cell.type.12.long', label = T)
x|y

# set default assay back to RNA
DefaultAssay(ecm_subset) <- "RNA"

# ----------------- convert to monacle3
cds <- as.cell_data_set(ecm_subset)

# to get cell metadata
colData(cds)

# to get gene metadata, set new col
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

# retrieve clustering information and store in new object
# assign partitions
reacreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

# assign the cluster information
list_cluster <- ecm_subset@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# assign the umap coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- ecm_subset@reductions$umap@cell.embeddings

# plot
cluster.before.traj <- plot_cells(cds, 
                                  color_cells_by = 'cluster', 
                                  label_groups_by_cluster = FALSE, 
                                  group_label_size = 5) + theme(legend.position = "right")

cluster.names <- plot_cells(cds,
           color_cells_by = "cell.type.12.long",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) + theme(legend.position = "right")

cluster.before.traj|cluster.names

# learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = "timepoint",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5,
           cell_size = 1.5
           ) + theme(legend.position = "right")

# order the cells in pseudotime
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,cds@colData@listData[["timepoint"]] == '0dpa']))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = TRUE,
           label_roots = FALSE,
           label_leaves = TRUE,
           group_label_size = 10)

# cells ordered by monocle3 pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, cell.type.12.long, fill = cell.type.12.long)) +
  geom_boxplot()

# find genes that change as a function ogf pseudotime
deg_ecm <- graph_test(cds, neighbor_graph = 'principal_graph', cores=4)
# filter for DE genes
deg_ecm %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

# showing genes DE across the pathway
FeaturePlot(ecm_subset, features = c('fabp11a', 'si:dkey-102g19.3', 'BX908782.2', 'ctsba', 'ctsd', 'si:ch211-198c19.3'))

# vizualizing pseudotimne in seurat
ecm_subset$pseudotime <- pseudotime(cds)
Idents(ecm_subset) <- ecm_subset$cell.type.12.long
FeaturePlot(ecm_subset, features = 'pseudotime', label = T)

# linear regression to test wether genes change overtime in expression
gene <- c("anxa2a")
cds_subset <- cds[rowData(cds)$gene_short_name %in% gene,]
gene_fits <- fit_models(cds_subset, model_formula_str = "~timepoint")
fit_coefs <- coefficient_table(gene_fits)

plot_genes_hybrid(cds_subset, group_cells_by="timepoint") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
