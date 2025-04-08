library(monocle3)
library(Seurat)
library(dyno)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(dplyr)

# -------------------- getting the data set up ------------------------ #
zf_filtered <- zf

DefaultAssay(ecm_subset) <- "RNA" # use RNA assay for normalization and selecting variable features

# set cell type as the ident for sub-setting
Idents(zf_filtered) <- zf_filtered$cell.type.12.long

# filter for ECM producing cells
bec_subset <- subset(zf_filtered, idents = c("Biliary Epithelial Cell"), 
                     subset = timepoint %in% c("untreated") == FALSE)

hep_subset <- subset(zf_filtered, idents = c("Hepatocyte"), 
                     subset = timepoint %in% c("untreated", "mock", "0dpa", "1dpa") == FALSE)

ecm_subset <- merge(hep_subset, bec_subset)

# After subsetting, drop unused levels from the cell type factor
ecm_subset$cell.type.12.long<- droplevels(ecm_subset$cell.type.12.long)

# pre-processing of data
# normalize the data
ecm_subset <- NormalizeData(ecm_subset)

# find variable features
ecm_subset <- FindVariableFeatures(ecm_subset, selection.method = "vst", nfeatures = 2000, loess.span = 0.5)

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
x <- DimPlot(ecm_subset, reduction = "umap", group.by = 'timepoint')
y <- DimPlot(ecm_subset, reduction = "umap", group.by = 'cell.type.12.long')
x|y

# set default assay back to RNA
DefaultAssay(ecm_subset) <- "RNA"

# ----------------- convert to monacle3 -------------------#
############################### USING MONOCLE3 CLUSTERING ############################
cds <- as.cell_data_set(ecm_subset)
# to get gene metadata, set new col
fData(cds)$gene_short_name <- rownames(fData(cds))
#cds <- preprocess_cds(cds, num_dim = 50)
#cds <- align_cds(cds, alignment_group = "orig.ident")
#plot_pc_variance_explained(cds)
#cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, reduction_method = 'UMAP')
cds <- learn_graph(cds, use_partition = TRUE)

# order the cells in pseudotime
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,cds@colData@listData[["timepoint"]] == 'mock']))


time <- plot_cells(cds, color_cells_by = 'timepoint',
           cell_size = 1,
           label_principal_points = TRUE) + theme(legend.position = "right")

type <- plot_cells(cds, color_cells_by = 'cell.type.12.long',
                   cell_size = 1,
                   label_roots = FALSE,
                   label_branch_points = TRUE,
                   label_leaves = TRUE) + theme(legend.position = "right")

time | type

# find genes that change as a function of pseudotime
deg_ecm <- graph_test(cds, neighbor_graph = 'principal_graph', cores=4)
# filter for DE genes
deg_ecm %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

cds_subset <- choose_cells(cds)

# find genes that change as a function of pseudotime
deg_ecm <- graph_test(cds_subset, neighbor_graph = 'principal_graph', cores=4)
# filter for DE genes
deg_ecm %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

plot_cells(cds, genes=c("hsp90aa1.2", "hspa5", "CABZ01080568.1", "sat1a.1", "cnn2", "mcl1a"),
           label_cell_groups=FALSE,
           label_leaves=FALSE)

# pull out DE genes
pr_deg_ids <- row.names(subset(deg_ecm, q_value < 0.05))


########################### USING SEURATS PRE DONE COMPUTATIONS ###############################
# to get cell metadata
colData(cds)

# to get gene metadata, set new col
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

# retrieve clustering information and store in new object
# assign partitions
recreate.partition <- as.factor(cds@colData$cell.type.12.long)
# Ensure the names match the cell names in the cds object
names(recreate.partition) <- rownames(cds@colData)
# Assign these partitions to the cds clusters
cds@clusters$UMAP$partitions <- recreate.partition

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
cds <- learn_graph(cds, use_partition = TRUE)

plot_cells(cds,
           color_cells_by = "cell.type.12.long",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 7,
           cell_size = 1.5,
           scale_to_range = FALSE
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

# find genes that change as a function of pseudotime
deg_ecm <- graph_test(cds, neighbor_graph = 'principal_graph', cores=4)
# filter for DE genes
deg_ecm %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

deg_ids <- row.names(subset(deg_ecm, q_value < 0.05))

# showing genes DE across the pathway
FeaturePlot(ecm_subset, features = c('cfd', 'krt94', 'mdka', 'anxa1a', 'sparc', 'gstm.3'))

# calculate modules of co-expressed genes
gene_modules <- find_gene_modules(cds[deg_ids,],
                                  resolution=c(10^seq(-6,-1)))
table(gene_modules$module)


# vizualizing pseudotimne in seurat
ecm_subset$pseudotime <- pseudotime(cds)
Idents(ecm_subset) <- ecm_subset$cell.type.12.long
FeaturePlot(ecm_subset, features = 'pseudotime', label = T)

# linear regression to test wether genes change overtime in expression
gene <- c("anxa2a")
cds_subset <- cds[rowData(cds)$gene_short_name %in% gene,]
gene_fits <- fit_models(cds_subset, model_formula_str = "~timepoint")
fit_coefs <- coefficient_table(gene_fits)

# how good is this model at evaluating gene expression?
evaluate_fits(gene_fits)
