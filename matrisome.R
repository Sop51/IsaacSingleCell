library(readr)
library(fgsea)
library(purrr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GseaVis)
library(Seurat)
library(edgeR)

# load in the matrisome gene set
matrisome_all <- read.csv('/Users/sm2949/Desktop/Dr_Matrisome_Masterlist_Nauroy et al_2017.xlsx - Dr_Matrisome_Masterlist.csv')

# separate into different ECM categories
collagens <- matrisome_all[matrisome_all$Matrisome.Category == "Collagens",]
collagen_genes <- collagens$Zebrafish.Gene.Symbol

ecm_regulators <- matrisome_all[matrisome_all$Matrisome.Category == "ECM Regulators",]
ecm_regulator_genes <- ecm_regulators$Zebrafish.Gene.Symbol

secreted_factors <- matrisome_all[matrisome_all$Matrisome.Category == "Secreted Factors",]
secreted_factors_genes <- secreted_factors$Zebrafish.Gene.Symbol

affiliated_proteins <- matrisome_all[matrisome_all$Matrisome.Category == "ECM-affiliated Proteins",]
affiliated_proteins_genes <- affiliated_proteins$Zebrafish.Gene.Symbol

proteoglycans <- matrisome_all[matrisome_all$Matrisome.Category == "Proteoglycans",]
proteoglycans_genes <- proteoglycans$Zebrafish.Gene.Symbol

ecm_glycoproteins <- matrisome_all[matrisome_all$Matrisome.Category == "ECM Glycoproteins",]
ecm_glycoproteins_genes <- ecm_glycoproteins$Zebrafish.Gene.Symbol

gene_sets <- list(Matrisome = matrisome_all$Zebrafish.Gene.Symbol)

# Create gene sets (gene subsets)
gene_sets <- list(
  Collagens = collagen_genes,
  ECM_Regulators = ecm_regulator_genes,
  Secreted_Factors = secreted_factors_genes,
  ECM_Affiliated_Proteins = affiliated_proteins_genes,
  Proteoglycans = proteoglycans_genes,
  ECM_Glycoproteins = ecm_glycoproteins_genes
)

# helper functions !!!!!!!
# Function to prepare ranked list by filtering p-value and sorting by log2FC
prepare_ranked_list <- function(df) {
    rankings <- sign(df$logFC)*(-log10(df$PValue))
    names(rankings) <- df$...1
    rankings <- sort(rankings, decreasing = TRUE)
    return(rankings)
}

# Function to run GSEA for each timepoint
run_gsea <- function(ranked_list, gene_sets) {
  gsea_results <- fgsea(pathways = gene_sets,
                        stats = ranked_list,
                        scoreType = 'neg',
                        minSize = 5,
                        maxSize = 500)
  return(gsea_results)
}

# Extract NES for each timepoint
get_nes <- function(gsea_results) {
  nes <- gsea_results$NES 
  nes_df <- data.frame(Pathway = gsea_results$pathway, NES = nes)
  return(nes_df)
}


# ------------------- perform for BECs across timepoints ---------------- #
bec0dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_bil_DE_results.csv")
bec1dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_bil_DE_results.csv")
bec2dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_bil_DE_results.csv")
bec3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_bil_DE_results.csv")
bec7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_bil_DE_results.csv")

# Prepare ranked gene lists for each timepoint after filtering by p-value < 0.05 and sorting by log2FC desc
bec0dpa_ranked <- prepare_ranked_list(bec0dpa_df)
bec1dpa_ranked <- prepare_ranked_list(bec1dpa_df)
bec2dpa_ranked <- prepare_ranked_list(bec2dpa_df)
bec3dpa_ranked <- prepare_ranked_list(bec3dpa_df)
bec7dpa_ranked <- prepare_ranked_list(bec7dpa_df)

# Run GSEA for each timepoint
gsea_bec0dpa <- run_gsea(bec0dpa_ranked, gene_sets)
gsea_bec1dpa <- run_gsea(bec1dpa_ranked, gene_sets)
gsea_bec2dpa <- run_gsea(bec2dpa_ranked, gene_sets)
gsea_bec3dpa <- run_gsea(bec3dpa_ranked, gene_sets)
gsea_bec7dpa <- run_gsea(bec7dpa_ranked, gene_sets)

# Extract NES for each timepoint
nes_bec0dpa <- get_nes(gsea_bec0dpa)
nes_bec1dpa <- get_nes(gsea_bec1dpa)
nes_bec2dpa <- get_nes(gsea_bec2dpa)
nes_bec3dpa <- get_nes(gsea_bec3dpa)
nes_bec7dpa <- get_nes(gsea_bec7dpa)

# Merge NES data from all timepoints
merged_nes_bec <- reduce(list(nes_bec0dpa, nes_bec1dpa, nes_bec2dpa, nes_bec3dpa, nes_bec7dpa), full_join, by = "Pathway")

# Rename columns to indicate timepoints
colnames(merged_nes_bec) <- c("Pathway", "bec0dpa", "bec1dpa", "bec2dpa", "bec3dpa", "bec7dpa")

# Filter for pathways that are present in all timepoints 
common_pathways_bec <- merged_nes_bec[complete.cases(merged_nes_bec[, -1]), ]

# Plot NES over time for common pathways
long_nes_bec <- common_pathways_bec %>%
  pivot_longer(cols = starts_with("bec"), names_to = "Timepoint", values_to = "NES")

# remove the cell name from the timepoint column
long_nes_bec$Timepoint <- gsub("^[a-zA-Z]+", "", long_nes_bec$Timepoint)

# plot
ggplot(long_nes_bec, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line(size = 1.0) + 
  geom_point(size = 1.5) +   
  scale_color_brewer(palette = "Set2") +  
  labs(title = "NES Over Liver Regeneration in BECs for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", 
       y = "Normalized Enrichment Score \n (NES)", 
       color = "Category") +
  theme_light(base_size = 12) +  # Cleaner theme with larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), # Rotated x-axis labels
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Centered title
    legend.position = "right", 
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)  # Thicker border
  )

# ------------------- perform for macrophages across timepoints ---------------- #
mac0dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_mac_DE_results.csv")
mac1dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_mac_DE_results.csv")
mac2dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_mac_DE_results.csv")
mac3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_mac_DE_results.csv")
mac7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_mac_DE_results.csv")

# Prepare ranked gene lists for each timepoint after filtering by p-value < 0.05 and sorting by log2FC desc
mac0dpa_ranked <- prepare_ranked_list(mac0dpa_df)
mac1dpa_ranked <- prepare_ranked_list(mac1dpa_df)
mac2dpa_ranked <- prepare_ranked_list(mac2dpa_df)
mac3dpa_ranked <- prepare_ranked_list(mac3dpa_df)
mac7dpa_ranked <- prepare_ranked_list(mac7dpa_df)

# Run GSEA for each timepoint
gsea_mac0dpa <- run_gsea(mac0dpa_ranked, gene_sets)
gsea_mac1dpa <- run_gsea(mac1dpa_ranked, gene_sets)
gsea_mac2dpa <- run_gsea(mac2dpa_ranked, gene_sets)
gsea_mac3dpa <- run_gsea(mac3dpa_ranked, gene_sets)
gsea_mac7dpa <- run_gsea(mac7dpa_ranked, gene_sets)

# Extract NES for each timepoint
nes_mac0dpa <- get_nes(gsea_mac0dpa)
nes_mac1dpa <- get_nes(gsea_mac1dpa)
nes_mac2dpa <- get_nes(gsea_mac2dpa)
nes_mac3dpa <- get_nes(gsea_mac3dpa)
nes_mac7dpa <- get_nes(gsea_mac7dpa)

# Merge NES data from all timepoints
merged_nes_mac <- reduce(list(nes_mac0dpa, nes_mac1dpa, nes_mac2dpa, nes_mac3dpa, nes_mac7dpa), full_join, by = "Pathway")

# Rename columns to indicate timepoints
colnames(merged_nes_mac) <- c("Pathway", "mac0dpa", "mac1dpa", "mac2dpa", "mac3dpa", "mac7dpa")

# Filter for pathways that are present in all timepoints 
common_pathways_mac <- merged_nes_mac[complete.cases(merged_nes_mac[, -1]), ]

# Plot NES over time for common pathways
long_nes_mac <- common_pathways_mac %>%
  pivot_longer(cols = starts_with("mac"), names_to = "Timepoint", values_to = "NES")

# remove the cell name from the timepoint column
long_nes_mac$Timepoint <- gsub("^[a-zA-Z]+", "", long_nes_mac$Timepoint)

# plot
ggplot(long_nes_mac, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line(size = 1.0) + 
  geom_point(size = 1.5) +   
  scale_color_brewer(palette = "Set2") +  
  labs(title = "NES Over Liver Regeneration in Macrophages for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", 
       y = "Normalized Enrichment Score \n (NES)", 
       color = "Category") +
  theme_light(base_size = 12) +  # Cleaner theme with larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), # Rotated x-axis labels
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Centered title
    legend.position = "right", 
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)  # Thicker border
  )

# ------------------- perform for macrophages across timepoints ---------------- #
hep0dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_hep_DE_results.csv")
hep1dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_hep_DE_results.csv")
hep2dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_hep_DE_results.csv")
hep3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_hep_DE_results.csv")
hep7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_hep_DE_results.csv")

# Prepare ranked gene lists for each timepoint after filtering by p-value < 0.05 and sorting by log2FC desc
hep0dpa_ranked <- prepare_ranked_list(hep0dpa_df)
hep1dpa_ranked <- prepare_ranked_list(hep1dpa_df)
hep2dpa_ranked <- prepare_ranked_list(hep2dpa_df)
hep3dpa_ranked <- prepare_ranked_list(hep3dpa_df)
hep7dpa_ranked <- prepare_ranked_list(hep7dpa_df)

# Run GSEA for each timepoint
gsea_hep0dpa <- run_gsea(hep0dpa_ranked, gene_sets)
gsea_hep1dpa <- run_gsea(hep1dpa_ranked, gene_sets)
gsea_hep2dpa <- run_gsea(hep2dpa_ranked, gene_sets)
gsea_hep3dpa <- run_gsea(hep3dpa_ranked, gene_sets)
gsea_hep7dpa <- run_gsea(hep7dpa_ranked, gene_sets)

# Extract NES for each timepoint
nes_hep0dpa <- get_nes(gsea_hep0dpa)
nes_hep1dpa <- get_nes(gsea_hep1dpa)
nes_hep2dpa <- get_nes(gsea_hep2dpa)
nes_hep3dpa <- get_nes(gsea_hep3dpa)
nes_hep7dpa <- get_nes(gsea_hep7dpa)

# Merge NES data from all timepoints
merged_nes_hep <- reduce(list(nes_hep0dpa, nes_hep1dpa, nes_hep2dpa, nes_hep3dpa, nes_hep7dpa), full_join, by = "Pathway")

# Rename columns to indicate timepoints
colnames(merged_nes_hep) <- c("Pathway", "hep0dpa", "hep1dpa", "hep2dpa", "hep3dpa", "hep7dpa")

# Filter for pathways that are present in all timepoints 
common_pathways_hep <- merged_nes_hep[complete.cases(merged_nes_hep[, -1]), ]

# Plot NES over time for common pathways
long_nes_hep <- common_pathways_hep %>%
  pivot_longer(cols = starts_with("hep"), names_to = "Timepoint", values_to = "NES")

# remove the cell name from the timepoint column
long_nes_hep$Timepoint <- gsub("^[a-zA-Z]+", "", long_nes_hep$Timepoint)

# plot
ggplot(long_nes_hep, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line(size = 1.0) + 
  geom_point(size = 1.5) +   
  scale_color_brewer(palette = "Set2") +  
  labs(title = "NES Over Liver Regeneration in Hepatocytes for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", 
       y = "Normalized Enrichment Score \n (NES)", 
       color = "Category") +
  theme_light(base_size = 12) +  # Cleaner theme with larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), # Rotated x-axis labels
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Centered title
    legend.position = "right", 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)  # Thicker border
  )
# ---------------------- create all plots in a row ---------------------------- #
# Add a column to indicate the cell type
long_nes_bec$CellType <- "BECs"
long_nes_mac$CellType <- "Macrophages"
long_nes_hep$CellType <- "Hepatocytes"

# Combine all datasets
combined_nes <- bind_rows(long_nes_bec, long_nes_mac, long_nes_hep)

# Plot with faceting
ggplot(combined_nes, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line(size = 1.0) + 
  geom_point(size = 1.5) +   
  scale_color_brewer(palette = "Set2") +  
  labs(title = "NES Over Liver Regeneration for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", 
       y = "Normalized Enrichment Score (NES)", 
       color = "Category") +
  theme_light(base_size = 12) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    legend.position = "right", 
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),  
    strip.text = element_text(size = 12, face = "bold", color = "black")  # Change facet label text color to black
  ) +
  facet_wrap(~CellType, nrow = 1)  # Facet in a single row
# ---------------------- creating a plot of ALL cell types ------------------------- #
long_nes_bec$CellType <- "BEC"
long_nes_mac$CellType <- "Macrophage"
long_nes_hep$CellType <- "Hepatocyte"

# combine into one df
long_nes_all <- bind_rows(long_nes_bec, long_nes_mac, long_nes_hep)

# remove the cell name from the timepoint column
long_nes_all$Timepoint <- gsub("^[a-zA-Z]+", "", long_nes_all$Timepoint)

# plot
ggplot(long_nes_all, aes(x = Timepoint, y = NES, color = CellType, group = CellType)) +
  geom_line(size = 1.2) +  
  geom_point(size = 2) +   
  scale_color_manual(values = c("Hepatocyte" = "#011a51", 
                                "BEC" = "#B6228A", 
                                "Macrophage" = "#f06c00")) + 
  labs(title = "NES Over Liver Regeneration for Matrisome Gene Sets Across Cell Types", 
       x = "Timepoint (dpa)", 
       y = "Normalized Enrichment Score \n (NES)", 
       color = "Cell Type") +
  theme_minimal(base_size = 12) +  # Increase base font size for readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"), # Rotate x-axis labels
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.title = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5), # Centered title
    legend.position = "right", 
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a box around the plot
  )
# -------------------------- create the heatmap --------------------------- #
# NOTE : THIS IS FOR HEPATOCYTES
# read in the file containing leading gene edge info for each cell type
leading_edges <- read_csv('/Users/sm2949/Desktop/v2SingleCellLeadingEdge.csv')

# reformat the categories 
leading_edges <- leading_edges %>%
  mutate(Category = case_when(
    Category == "ecm affiliated proteins" ~ "ECM-affiliated proteins",
    Category == "ecm glycoproteins" ~ "ECM glycoproteins",
    Category == "ecm regulators" ~ "ECM regulators",
    Category == "secreted factors" ~ "Secreted factors",
    TRUE ~ Category # Keep other values unchanged
  )) 

# NOTE: THIS SECTION IS FOR HEPATOCYTES 
# extract the genes for each timepoint
leading_0dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Hepatocyte", '0dpa'])
leading_1dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Hepatocyte", '1dpa'])
leading_2dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Hepatocyte", '2dpa'])
leading_3dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Hepatocyte", '3dpa'])
leading_7dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Hepatocyte", '7dpa'])

# combine into one list 
all_leading_genes_hep <- c(leading_0dpa_genes, 
                           leading_1dpa_genes, 
                           leading_2dpa_genes, 
                           leading_3dpa_genes, 
                           leading_7dpa_genes)

# subset the count matrix for only the wanted genes
leading_edge_hep <- normalized_hep[rownames(normalized_hep) %in% all_leading_genes_hep, ]
leading_edge_hep <- as.data.frame(leading_edge_hep)
leading_edge_hep <- leading_edge_hep[, !colnames(leading_edge_hep) %in% c("mock")]
leading_edge_hep <- as.matrix(leading_edge_hep)

# Create a metadata frame to hold timepoints and KEGG categories for the genes
gene_metadata <- leading_edges %>% 
  filter(`cell type` == "Hepatocyte") %>%
  select('0dpa', '1dpa', '2dpa', '3dpa', '7dpa', 'Category') %>%
  gather(key = "timepoint", value = "gene", -Category) %>%
  filter(gene %in% rownames(leading_edge_hep)) %>%
  distinct(gene, Category, .keep_all = TRUE)

# Create a new data frame that matches the order of genes in `leading_edge_hep`
ordered_gene_metadata <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_hep)) %>%
  arrange(match(gene, rownames(leading_edge_hep))) %>%
  select(gene, Category)

# Set the gene column as the row names
ordered_gene_metadata <- as.data.frame(ordered_gene_metadata)
rownames(ordered_gene_metadata) <- ordered_gene_metadata$gene
ordered_gene_metadata <- ordered_gene_metadata %>% select(-gene)

annotation_colors <- list(
  Category = c("ECM-affiliated proteins" = "#eab69f",
               "ECM glycoproteins" = "#e07a5f", 
               "ECM regulators" = "#3d405b",
               "Secreted factors" = "#f2cc8f")
)

ha <- HeatmapAnnotation(Category = ordered_gene_metadata$Category,
                        which = 'row',
                        col = annotation_colors)

timepoint_split <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_hep)) %>%
  arrange(match(gene, rownames(leading_edge_hep))) %>%
  pull(timepoint)

# Z-score normalize the matrix by row (genes)
leading_edge_hep_zscored <- t(scale(t(leading_edge_hep)))

# remove columns that start with "mock" or "untreated"
leading_edge_hep_zscored <- leading_edge_hep_zscored %>%
  as.data.frame() %>%
  select(-starts_with("mock"), -starts_with("untreated")) 

# clean up the column names to remove the library
clean_names <- gsub("LIB11[[:alpha:]]$", "", colnames(leading_edge_hep_zscored))

# ensure uniqueness by appending .rep1, .rep2, etc.
colnames(leading_edge_hep_zscored) <- ave(clean_names, clean_names, FUN = function(x) {
  if (length(x) > 1) paste0(x, ".rep", seq_along(x)) else x
})

leading_edge_hep_zscored <- as.matrix(leading_edge_hep_zscored)

# generate heatmap with annotations
Heatmap(leading_edge_hep_zscored, right_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2),
        column_title = "Leading Edge Matrisome Genes Per Pathway In Hepatocytes",
        heatmap_legend_param = list(
          title = gt_render("<span style='color:black'>**Expression**</span>"), 
          at = c(-2, 0, 2), 
          labels = gt_render(c("Low Expression", "No Change", "High Expression"))
        ))

# NOTE : THIS IS FOR MACROPHAGES
# extract the genes for each timepoint
leading_0dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Macrophage", '0dpa'])
leading_1dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Macrophage", '1dpa'])
leading_2dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Macrophage", '2dpa'])
leading_3dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Macrophage", '3dpa'])
leading_7dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "Macrophage", '7dpa'])

# combine into one list 
all_leading_genes_mac <- c(leading_0dpa_genes, 
                           leading_1dpa_genes, 
                           leading_2dpa_genes, 
                           leading_3dpa_genes, 
                           leading_7dpa_genes)

# subset the count matrix for only the wanted genes
leading_edge_mac <- normalized_mac[rownames(normalized_mac) %in% all_leading_genes_mac, ]
leading_edge_mac <- as.data.frame(leading_edge_mac)
leading_edge_mac <- leading_edge_mac[, !colnames(leading_edge_mac) %in% c("mock")]
leading_edge_mac <- as.matrix(leading_edge_mac)

# Create a metadata frame to hold timepoints and KEGG categories for the genes
gene_metadata <- leading_edges %>% 
  filter(`cell type` == "Macrophage") %>%
  select('0dpa', '1dpa', '2dpa', '3dpa', '7dpa', 'Category') %>%
  gather(key = "timepoint", value = "gene", -Category) %>%
  filter(gene %in% rownames(leading_edge_mac)) %>%
  distinct(gene, Category, .keep_all = TRUE)

# Create a new data frame that matches the order of genes in `leading_edge_mac`
ordered_gene_metadata <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_mac)) %>%
  arrange(match(gene, rownames(leading_edge_mac))) %>%
  select(gene, Category)

# Set the gene column as the row names
ordered_gene_metadata <- as.data.frame(ordered_gene_metadata)
rownames(ordered_gene_metadata) <- ordered_gene_metadata$gene
ordered_gene_metadata <- ordered_gene_metadata %>% select(-gene)


ha <- HeatmapAnnotation(Category = ordered_gene_metadata$Category,
                        which = 'row',
                        col = annotation_colors)

timepoint_split <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_mac)) %>%
  arrange(match(gene, rownames(leading_edge_mac))) %>%
  pull(timepoint)

# Z-score normalize the matrix by row (genes)
leading_edge_mac_zscored <- t(scale(t(leading_edge_mac)))

# remove columns that start with "mock" or "untreated"
leading_edge_mac_zscored <- leading_edge_mac_zscored %>%
  as.data.frame() %>%
  select(-starts_with("mock"), -starts_with("untreated")) 

# clean up the column names to remove the library
clean_names <- gsub("LIB11[[:alpha:]]$", "", colnames(leading_edge_mac_zscored))

# ensure uniqueness by appending .rep1, .rep2, etc.
colnames(leading_edge_mac_zscored) <- ave(clean_names, clean_names, FUN = function(x) {
  if (length(x) > 1) paste0(x, ".rep", seq_along(x)) else x
})

leading_edge_mac_zscored <- as.matrix(leading_edge_mac_zscored)

# generate heatmap with annotations
Heatmap(leading_edge_mac_zscored, right_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2),
        column_title = "Leading Edge Matrisome Genes Per Pathway In Macrophages",
        heatmap_legend_param = list(
          title = gt_render("<span style='color:black'>**Expression**</span>"), 
          at = c(-2, 0, 2), 
          labels = gt_render(c("Low Expression", "No Change", "High Expression"))
        ))

# NOTE : THIS IS FOR BECS
# extract the genes for each timepoint
leading_0dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "BEC", '0dpa'])
leading_1dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "BEC", '1dpa'])
leading_2dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "BEC", '2dpa'])
leading_3dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "BEC", '3dpa'])
leading_7dpa_genes <- unlist(leading_edges[leading_edges$`cell type` == "BEC", '7dpa'])

# combine into one list 
all_leading_genes_bec <- c(leading_0dpa_genes, 
                           leading_1dpa_genes, 
                           leading_2dpa_genes, 
                           leading_3dpa_genes, 
                           leading_7dpa_genes)

# subset the count matrix for only the wanted genes
leading_edge_bec <- normalized_bec[rownames(normalized_bec) %in% all_leading_genes_bec, ]
leading_edge_bec <- as.data.frame(leading_edge_bec)
leading_edge_bec <- leading_edge_bec[, !colnames(leading_edge_bec) %in% c("mock")]
leading_edge_bec <- as.matrix(leading_edge_bec)

# Create a metadata frame to hold timepoints and KEGG categories for the genes
gene_metadata <- leading_edges %>% 
  filter(`cell type` == "BEC") %>%
  select('0dpa', '1dpa', '2dpa', '3dpa', '7dpa', 'Category') %>%
  gather(key = "timepoint", value = "gene", -Category) %>%
  filter(gene %in% rownames(leading_edge_bec)) %>%
  distinct(gene, Category, .keep_all = TRUE)

# Create a new data frame that matches the order of genes in `leading_edge_bec`
ordered_gene_metadata <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_bec)) %>%
  arrange(match(gene, rownames(leading_edge_bec))) %>%
  select(gene, Category)

# Set the gene column as the row names
ordered_gene_metadata <- as.data.frame(ordered_gene_metadata)
rownames(ordered_gene_metadata) <- ordered_gene_metadata$gene
ordered_gene_metadata <- ordered_gene_metadata %>% select(-gene)


ha <- HeatmapAnnotation(Category = ordered_gene_metadata$Category,
                        which = 'row',
                        col = annotation_colors)

timepoint_split <- gene_metadata %>%
  filter(gene %in% rownames(leading_edge_bec)) %>%
  arrange(match(gene, rownames(leading_edge_bec))) %>%
  pull(timepoint)

# Z-score normalize the matrix by row (genes)
leading_edge_bec_zscored <- t(scale(t(leading_edge_bec)))

# remove columns that start with "mock" or "untreated"
leading_edge_bec_zscored <- leading_edge_bec_zscored %>%
  as.data.frame() %>%
  select(-starts_with("mock"), -starts_with("untreated")) 

# clean up the column names to remove the library
clean_names <- gsub("LIB11[[:alpha:]]$", "", colnames(leading_edge_bec_zscored))

# ensure uniqueness by appending .rep1, .rep2, etc.
colnames(leading_edge_bec_zscored) <- ave(clean_names, clean_names, FUN = function(x) {
  if (length(x) > 1) paste0(x, ".rep", seq_along(x)) else x
})

leading_edge_bec_zscored <- as.matrix(leading_edge_bec_zscored)

# generate heatmap with annotations
Heatmap(leading_edge_bec_zscored, right_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2),
        column_title = "Leading Edge Matrisome Genes Per Pathway In BECs",
        heatmap_legend_param = list(
          title = gt_render("<span style='color:black'>**Expression**</span>"), 
          at = c(-2, 0, 2), 
          labels = gt_render(c("Low Expression", "No Change", "High Expression"))
        ))

