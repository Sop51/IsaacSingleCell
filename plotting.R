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
UMAPPlot(zf) + 
  scale_color_manual(values = okabe_ito_12)

# plot gene expression in UMAP
FeaturePlot(zf, features = c("anxa2a"))

# differential gene expression
cell_type <- 'Macrophage'
diff_table <- FindMarkers(zf, ident.1 = cell_type)