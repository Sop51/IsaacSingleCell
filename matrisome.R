library(readr)
library(fgsea)
library(purrr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GseaVis)

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

background_genes <- rownames(zf)

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
                        scoreType = 'std',
                        minSize = 10,
                        maxSize = 500,
                        nproc = 1,
                        geneL)
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

ggplot(long_nes_bec, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In BECs for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

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

ggplot(long_nes_mac, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In Macrophages for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

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

plotEnrichment(gene_sets[["Secreted_Factors"]],
               hep1dpa_ranked) + labs(title="Matrisome")


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

ggplot(long_nes_hep, aes(x = Timepoint, y = NES, color = Pathway, group = Pathway)) +
  geom_line() +
  geom_point() +
  labs(title = "NES Over Regeneration In Hepatocytes for Matrisome Gene Sets", 
       x = "Timepoint (dpa)", y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()

