library(Seurat)
library(ggplot2)
library(viridis)
        
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

# plot side by side
umap | feature

# differential gene expression
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

# gene of interest
gene <- "spint2"

# check for significant adjusted p-value in each table
significant_tables <- sapply(names(de_tables), function(ct) {
  table <- de_tables[[ct]]
  if (gene %in% rownames(table)) {
    return(table[gene, "p_val_adj"] < 0.05)
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

