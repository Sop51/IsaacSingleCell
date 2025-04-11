library(ComplexHeatmap)

becdpa1vsdpa2 <- read_csv("/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa1vsdpa2_bil_DE_results.csv")
becdpa2vsdpa3 <- read_csv("/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa2vsdpa3_bil_DE_results.csv")

# filter for significant genes
filter_significant_genes <- function(df, gene_list) {
  df <- df[df$PValue < 0.05 & df$logFC > 1, ]  # Filter for significance and logFC
  return(df)
}

becdpa1vsdpa2_sig <- becdpa1vsdpa2[becdpa1vsdpa2$PValue < 0.05 & abs(becdpa1vsdpa2$logFC) > 1.5, ]
becdpa2vsdpa3_sig <- becdpa2vsdpa3[becdpa2vsdpa3$PValue < 0.05 & abs(becdpa2vsdpa3$logFC) > 1.5, ]

# pull out genes in the matrisome gene list
# load in the matrisome gene set
matrisome_all <- read.csv('/Users/sophiemarcotte/Desktop/Dr_Matrisome_Masterlist_Nauroy et al_2017.xlsx - Dr_Matrisome_Masterlist.csv')
zebrafishMatrisomeGenes <- matrisome_all$Zebrafish.Gene.Symbol

sig.mat.dpa1vsdpa2 <- becdpa1vsdpa2_sig[becdpa1vsdpa2_sig$...1 %in% zebrafishMatrisomeGenes, ]
sig.mat.dpa2vsdpa3 <- becdpa2vsdpa3_sig[becdpa2vsdpa3_sig$...1 %in% zebrafishMatrisomeGenes, ]

# create the matrix
logFC_matrix1vs2 <- as.matrix(sig.mat.dpa1vsdpa2$logFC)
rownames(logFC_matrix1vs2) <- sig.mat.dpa1vsdpa2$...1

# plot
Heatmap(logFC_matrix1vs2, 
        row_names_gp = gpar(fontsize = 8),
        name = "logFC", 
        show_row_names = TRUE, 
        show_column_names = FALSE,
        col = colorRampPalette(c("navy", "white", "orange"))(100), 
        cluster_rows = TRUE,  # Don't cluster rows
        cluster_columns = FALSE,  # Don't cluster columns
        column_title = "logFC Values of DE Matrisome Genes in BECs 2dpa vs dpa" 
)
