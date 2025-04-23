library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)
        
# load in the seurat obj
file_path <- '/Users/sm2949/Desktop/20211217_zf.mtz2.3.cell.types.Robj'
load(file_path)

# okabe-Ito palette extended to 12 colors
okabe_ito_12 <- c(
  "#E69F00", "#56B4E9", "#003434", "#F0E442", "#0072B2", 
  "#D55E00", "#CC79A7", "#999999", "#882255", "#44AA99", 
  "#117733", "#332288"
)

# UMAP plot
umap <- UMAPPlot(zf) + 
  scale_color_manual(values = okabe_ito_12)

# plot gene expression in UMAP
feature <- FeaturePlot(zf, features = c("lima1a"))

VlnPlot(zf, features = c('uhrf1'), group.by = 'cell.type.12.long')

# plot side by side
umap | feature

#--------------------- differential gene expression ------------------------ #
macrophage_table <- FindMarkers(zf, ident.1 = 'Macrophage', only.pos = TRUE)
bec_table <- FindMarkers(zf, ident.1 = 'Biliary Epithelial Cell', only.pos = TRUE)
end_table <- FindMarkers(zf, ident.1 = 'Endothelial Cell', only.pos = TRUE)
hep_table <- FindMarkers(zf, ident.1 = 'Hepatocyte', only.pos = TRUE)
lymph_table <- FindMarkers(zf, ident.1 = 'Lymphocyte', only.pos = TRUE)
neutrophil_table <- FindMarkers(zf, ident.1 = 'Neutrophil', only.pos = TRUE)
ery_table <- FindMarkers(zf, ident.1 = 'Erythrocyte', only.pos = TRUE)
neu_table <- FindMarkers(zf, ident.1 = 'Neuron', only.pos = TRUE)
islet_table <- FindMarkers(zf, ident.1 = 'Islet Cell', only.pos = TRUE)
apln_table <- FindMarkers(zf, ident.1 = 'apln+', only.pos = TRUE)
spink_table <- FindMarkers(zf, ident.1 = 'spink+', only.pos = TRUE)
fthl_table <- FindMarkers(zf, ident.1 = 'fthl+', only.pos = TRUE)

# store tables in a named list
de_tables <- list(
  Macrophage = macrophage_table,
  Biliary_Epithelial_Cell = bec_table,
  Endothelial_Cell = end_table,
  Hepatocyte = hep_table,
  Lymphocyte = lymph_table,
  Neutrophil = neutrophil_table,
  Erythrocyte = ery_table,
  Neuron = neu_table,
  Islet_Cell = islet_table,
  apln = apln_table,
  spink = spink_table,
  fthl = fthl_table
)

# ----------------- small function to identify which cell types genes are sig in ------------------ #
# gene of interest
gene <- "spint2"

# check for significant adjusted p-value in each table
significant_tables <- sapply(names(de_tables), function(ct) {
  table <- de_tables[[ct]]
  if (gene %in% rownames(table)) {
    return(table[gene, "p_val_adj"] < 1.1)
  } else {
    return(FALSE)
  }
})

# print the tables where the gene is significant
significant_cell_types <- names(significant_tables[significant_tables])
if (length(significant_cell_types) > 0) {
  cat("The gene", gene, "has a significant adjusted p-value in:\n")
  print(significant_cell_types)
} else {
  cat("The gene", gene, "is not significantly differentially expressed in any table.\n")
}

# ------------- plotting cell type proportions across timepoints ------------------ #
dim <- DimPlot(zf, cols = custom)
# create a custom color pallete
custom <- c(
  "#011a51",
  "#B6228A",
  "#fab129",
  "#f06c00",
  "#b83700",
  "#143D0C",
  "#94B622",
  "#AAB8D6",
  "#B68222",
  "#27167F",
  "#DFACCC",
  "#FDE335"
)

# extract relevant metadata
metadata <- zf@meta.data

# count the number of cells for each combination of cell type and timepoint
cell_counts <- metadata %>%
  group_by(cell.type.12.long, timepoint) %>%
  summarise(cell_count = n()) %>%
  ungroup()

# stacked bar chart 
stack <- ggplot(cell_counts, aes(x = timepoint, y = cell_count, fill = cell.type.12.long)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked proportions
  labs(
    title = "Cell Type Proportions Across Timepoints",
    x = "Timepoint",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  scale_fill_manual(values = custom) +  # Match colors to custom palette
  theme_minimal() +  # Clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Larger x-axis labels
    axis.text.y = element_text(size = 14),  # Larger y-axis labels
    axis.title.x = element_text(size = 16, face = "bold"),  # Larger x-axis title
    axis.title.y = element_text(size = 16, face = "bold"),  # Larger y-axis title
    legend.text = element_text(size = 14),  # Larger legend text
    legend.title = element_text(size = 16, face = "bold"),  # Larger legend title
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")  # Larger plot title
  )

dim | stack

# ------------------------------- plotting cell type count ------------------------------ #
# count number of cells per cell type
cell_counts <- zf@meta.data %>%
  count(cell.type.12.long) %>%
  arrange(desc(n))

# car plot
ggplot(cell_counts, aes(x = reorder(cell.type.12.long, n), y = n, fill = cell.type.12.long)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 4.5) +  # add cell counts above bars
  scale_fill_viridis_d(option = "D", direction = -1) +  # colorblind-friendly discrete palette
  labs(
    title = "Number of Cells per Cell Type",
    x = "Cell Type",
    y = "Cell Count",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"  # hide legend
  ) +
  ylim(0, max(cell_counts$n) * 1.1)  # add space above bars for labels
