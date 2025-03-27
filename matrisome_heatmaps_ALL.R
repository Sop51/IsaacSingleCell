library(readr)
library(pheatmap)

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

all_matrisome_genes <- list(Matrisome = matrisome_all$Zebrafish.Gene.Symbol)

# helper functions !!!!!!!
filter_matrisome_genes <- function(df) {
  matrisome_genes <- all_matrisome_genes$Matrisome
  df_filtered <- df[df$logFC > 1 & df$PValue < 0.05, ]  # Filter for significance
  df_filtered <- df_filtered[df_filtered$...1 %in% matrisome_genes , ]  # Keep only matrisome genes
  return(df_filtered)
}

generate_heatmap_matrix <- function(timepoint_dfs){
  # extract gene names
  all_genes <- unique(unlist(lapply(timepoint_dfs, function(df) df$...1)))
  # create a matrix to store log2FC values
  heatmap_matrix <- data.frame(Gene = all_genes)
  # fill in log2FC values for each timepoint
  for (timepoint in names(timepoint_dfs)) {
    df <- timepoint_dfs[[timepoint]]
    heatmap_matrix[[timepoint]] <- df$logFC[match(heatmap_matrix$Gene, df$...1)]
  }
  # replace NAs with 0 
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  # convert to matrix format
  rownames(heatmap_matrix) <- heatmap_matrix$Gene
  heatmap_matrix <- as.matrix(heatmap_matrix[, -1])
  return(heatmap_matrix)
}

# ------------------- perform for BECs across timepoints ---------------- #
bec0dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_bil_DE_results.csv")
bec1dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_bil_DE_results.csv")
bec2dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_bil_DE_results.csv")
bec3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_bil_DE_results.csv")
bec7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_bil_DE_results.csv")

# prepare ranked gene lists for each timepoint after filtering by p-value < 0.05 and sorting by log2FC desc
bec0dpa_ranked <- filter_matrisome_genes(bec0dpa_df)
bec1dpa_ranked <- filter_matrisome_genes(bec1dpa_df)
bec2dpa_ranked <- filter_matrisome_genes(bec2dpa_df)
bec3dpa_ranked <- filter_matrisome_genes(bec3dpa_df)
bec7dpa_ranked <- filter_matrisome_genes(bec7dpa_df)

# list of ranked timepoint dataframes
timepoint_bec_dfs <- list(
  "0dpa" = bec0dpa_ranked,
  "1dpa" = bec1dpa_ranked,
  "2dpa" = bec2dpa_ranked,
  "3dpa" = bec3dpa_ranked,
  "7dpa" = bec7dpa_ranked
)

# generate heatmap
pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  cluster_rows = TRUE, 
  cluster_cols = FALSE,  
  scale = "row",  
  fontsize_row = 7,  
  fontsize_col = 10,
  main = "Differentially Expressed Matrisome Genes in BECs Across Regeneration"
)

# ------------------- perform for BECs across timepoints ---------------- #
mac0dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_mac_DE_results.csv")
mac1dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_mac_DE_results.csv")
mac2dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_mac_DE_results.csv")
mac3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_mac_DE_results.csv")
mac7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_mac_DE_results.csv")

# prepare ranked gene lists for each timepoint after filtering by p-value < 0.05 and sorting by log2FC desc
mac0dpa_ranked <- filter_matrisome_genes(mac0dpa_df)
mac1dpa_ranked <- filter_matrisome_genes(mac1dpa_df)
mac2dpa_ranked <- filter_matrisome_genes(mac2dpa_df)
mac3dpa_ranked <- filter_matrisome_genes(mac3dpa_df)
mac7dpa_ranked <- filter_matrisome_genes(mac7dpa_df)

# list of ranked timepoint dataframes
timepoint_mac_dfs <- list(
  "0dpa" = mac0dpa_ranked,
  "1dpa" = mac1dpa_ranked,
  "2dpa" = mac2dpa_ranked,
  "3dpa" = mac3dpa_ranked,
  "7dpa" = mac7dpa_ranked
)

mac_heatmap_matrix <- generate_heatmap_matrix(timepoint_mac_dfs)

# generate heatmap
pheatmap(
  mac_heatmap_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),  
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  scale = "row",  
  fontsize_row = 7,  
  fontsize_col = 10,
  main = "Differentially Expressed Matrisome Genes in Macrophages Across Regeneration"
)

# ------------------- perform for hpeatocytes across timepoints ---------------- #
hep0dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_hep_DE_results.csv")
hep1dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_hep_DE_results.csv")
hep2dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_hep_DE_results.csv")
hep3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_hep_DE_results.csv")
hep7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_hep_DE_results.csv")

# prepare ranked gene lists for each timepoint after filtering by p-value < 0.05 and sorting by log2FC desc
hep0dpa_ranked <- filter_matrisome_genes(hep0dpa_df)
hep1dpa_ranked <- filter_matrisome_genes(hep1dpa_df)
hep2dpa_ranked <- filter_matrisome_genes(hep2dpa_df)
hep3dpa_ranked <- filter_matrisome_genes(hep3dpa_df)
hep7dpa_ranked <- filter_matrisome_genes(hep7dpa_df)

# list of ranked timepoint dataframes
timepoint_hep_dfs <- list(
  "0dpa" = hep0dpa_ranked,
  "1dpa" = hep1dpa_ranked,
  "2dpa" = hep2dpa_ranked,
  "3dpa" = hep3dpa_ranked,
  "7dpa" = hep7dpa_ranked
)

hep_heatmap_matrix <- generate_heatmap_matrix(timepoint_hep_dfs)

# generate heatmap
pheatmap(
  hep_heatmap_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),  
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  scale = "row",  
  fontsize_row = 7,  
  fontsize_col = 10,
  main = "Differentially Expressed Matrisome Genes in Hepatocytes Across Regeneration"
)
