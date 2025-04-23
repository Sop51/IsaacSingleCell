library(ComplexHeatmap)
library(readr)
library(Seurat)
library(dplyr)

becdpa1vsdpa2 <- read_csv("/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa1vsdpa2_bil_DE_results.csv")
becdpa2vsdpa3 <- read_csv("/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa2vsdpa3_bil_DE_results.csv")

# filter for significant genes
filter_significant_genes <- function(df, gene_list) {
  df <- df[df$PValue < 0.05 & df$logFC > 1, ]  # Filter for significance and logFC
  return(df)
}

becdpa1vsdpa2_sig <- becdpa1vsdpa2[becdpa1vsdpa2$PValue < 0.05 & becdpa1vsdpa2$logFC < 1.0, ]
becdpa2vsdpa3_sig <- becdpa2vsdpa3[becdpa2vsdpa3$PValue < 0.05 & becdpa2vsdpa3$logFC < 1.0, ]

# pull out genes in the matrisome gene list
# load in the matrisome gene set
matrisome_all <- read.csv('/Users/sophiemarcotte/Desktop/Dr_Matrisome_Masterlist_Nauroy et al_2017.xlsx - Dr_Matrisome_Masterlist.csv')
matrisomeGenes <- matrisome_all$Zebrafish.Gene.Symbol

# filter for only core matrisome genes
core_matrisome <- matrisome_all %>% filter(Matrisome.Division == "Core matrisome")
coreMatrisomeGenes <- core_matrisome$Zebrafish.Gene.Symbol

# filter for only associated matrisome genes
associated_matrisome <- matrisome_all %>% filter(Matrisome.Division == "Matrisome-associated")
associatedMatrisomeGenes <- associated_matrisome$Zebrafish.Gene.Symbol

# pull out the significant core matrisome genes
sig.core.mat.dpa1vsdpa2 <- becdpa1vsdpa2_sig[becdpa1vsdpa2_sig$...1 %in% coreMatrisomeGenes, ]
sig.core.mat.dpa2vsdpa3 <- becdpa2vsdpa3_sig[becdpa2vsdpa3_sig$...1 %in% coreMatrisomeGenes, ]

sig_core_gene_list1 <- sig.core.mat.dpa1vsdpa2$...1
sig_core_gene_list2 <- sig.core.mat.dpa2vsdpa3$...1

# combine the two gene lists
combined_matrisome_list <- unique(c(sig_core_gene_list1, sig_core_gene_list2))

# pull out the significant associated matrisome genes
sig.associated.mat.dpa1vsdpa2 <- becdpa1vsdpa2_sig[becdpa1vsdpa2_sig$...1 %in% associatedMatrisomeGenes, ]
sig.associated.mat.dpa2vsdpa3 <- becdpa2vsdpa3_sig[becdpa2vsdpa3_sig$...1 %in% associatedMatrisomeGenes, ]

sig_associated_gene_list1 <- sig.associated.mat.dpa1vsdpa2$...1
sig_associated_gene_list2 <- sig.associated.mat.dpa2vsdpa3$...1

# combine the two gene lists
combined_associated_matrisome_list <- unique(c(sig_associated_gene_list1, sig_associated_gene_list2))

# pull out the significant ALL matrisome genes
sig.mat.dpa1vsdpa2 <- becdpa1vsdpa2_sig[becdpa1vsdpa2_sig$...1 %in% matrisomeGenes, ]
sig.mat.dpa2vsdpa3 <- becdpa2vsdpa3_sig[becdpa2vsdpa3_sig$...1 %in% matrisomeGenes, ]

sig_gene_list1 <- sig.mat.dpa1vsdpa2$...1
sig_gene_list2 <- sig.mat.dpa2vsdpa3$...1

# combine the two gene lists
all_combined_matrisome_list <- unique(c(sig_gene_list1, sig_gene_list2))

# pull out the top20 up and down regulated genes
#top20_dpa1vsdpa2 <- head(becdpa1vsdpa2_sig[order(becdpa1vsdpa2_sig$PValue), "...1"], 20)
#top20_dpa2vsdpa3 <- head(becdpa2vsdpa3_sig[order(becdpa2vsdpa3_sig$PValue), "...1"], 20)

# pull out genes into lists
#top201 <- top20_dpa1vsdpa2$...1
#top202 <- top20_dpa2vsdpa3$...1

# combine and remove duplicates
combined_top_genes <- unique(c(top201, top202))

# subset the bec object to only include the wanted timepoints
bec_subset <- subset(bec, subset = timepoint %in% c("1dpa", "2dpa", "3dpa"))
bec_subset <- ScaleData(bec_subset, features = all_combined_matrisome_list, assay = "RNA")
bec_subset$timepoint <- droplevels(bec_subset$timepoint)

# plot the matrisome only heatmap
DoHeatmap(bec_subset, 
          features = all_combined_matrisome_list, 
          assay = "RNA", 
          slot = "scale.data", 
          group.by = "timepoint") +
  theme(axis.text.y = element_text(size = 6)) 

# plot the top DE genes heatmap
DoHeatmap(bec_subset, 
          features = combined_top_genes, 
          assay = "RNA", 
          slot = "scale.data", 
          group.by = "timepoint") +
  theme(axis.text.y = element_text(size = 6)) 
