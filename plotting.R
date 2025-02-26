library(Seurat)

# load in the seurat obj
file_path <- '/Users/sophiemarcotte/Desktop/20211217_zf.mtz2.3.cell.types.Robj'
load(file_path)

# UMAP plot
UMAPPlot(zf)

# plot gene expression in UMAP
FeaturePlot(zf, features = c("anxa2a"))

# differential gene expression
cell_type <- 'Macrophage'
diff_table <- FindMarkers(zf, ident.1 = cell_type)