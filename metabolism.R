library(KEGGREST)
library(dplyr)
library(pheatmap)

# read in the differential expression files
hep3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_hep_DE_results.csv")
hep7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_hep_DE_results.csv")
bec3dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_bil_DE_results.csv")
bec7dpa_df <- read_csv("/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_bil_DE_results.csv")

# store dataframes in a named list for easy access
de_dfs <- list(
  "hep3dpa" = hep3dpa_df,
  "hep7dpa" = hep7dpa_df,
  "bec3dpa" = bec3dpa_df,
  "bec7dpa" = bec7dpa_df
) 

# define pathways of interest
pathways <- c("dre00010", "dre00071", "dre00061", "dre00020", "dre00120", 
              "dre00561", "dre00564", "dre00980")

pathway_names <- c("Glycolysis / Gluconeogenesis", "Fatty acid degradation", "Fatty acid biosynthesis",
                   "Citrate cycle (TCA cycle)", "Primary bile acid biosynthesis", "Glycerolipid metabolism",
                   "Glycerophospholipid metabolism", "Metabolism of xenobiotics by cytochrome P450")

# function to get genes from a KEGG pathway
get_kegg_genes <- function(pathway) {
  genes <- keggGet(pathway)[[1]]$GENE
  if (is.null(genes)) return(NULL)
  genes <- genes[seq(1, length(genes), by = 2)]  
  return(genes)
}

# retrieve genes for all pathways and store as a named list
kegg_genes_list <- setNames(lapply(pathways, get_kegg_genes), pathways)

# Function to convert Entrez IDs to Gene Symbols
convert_to_symbols <- function(entrez_ids) {
  symbols <- mapIds(org.Dr.eg.db, 
                    keys = entrez_ids, 
                    column = "SYMBOL", 
                    keytype = "ENTREZID", 
                    multiVals = "first")  # Get the first mapped value
  return(symbols)
}

# Apply the function to each pathway
gene_symbols_list <- lapply(kegg_genes_list, convert_to_symbols)

# pull out the genes from the de tables for each pathway 
dre00010_genes <- gene_symbols_list$dre00010
bec3dpa_dre00010 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00010_genes)
bec7dpa_dre00010 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00010_genes)
hep3dpa_dre00010 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00010_genes)
hep7dpa_dre00010 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00010_genes)

dre00071_genes <- gene_symbols_list$dre00071
bec3dpa_dre00071 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00071_genes)
bec7dpa_dre00071 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00071_genes)
hep3dpa_dre00071 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00071_genes)
hep7dpa_dre00071 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00071_genes)

dre000611_genes <- gene_symbols_list$dre00061
bec3dpa_dre000611 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre000611_genes)
bec7dpa_dre000611 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre000611_genes)
hep3dpa_dre000611 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre000611_genes)
hep7dpa_dre000611 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre000611_genes)

dre00020_genes <- gene_symbols_list$dre00020
bec3dpa_dre00020 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00020_genes)
bec7dpa_dre00020 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00020_genes)
hep3dpa_dre00020 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00020_genes)
hep7dpa_dre00020 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00020_genes)

dre00120_genes <- gene_symbols_list$dre00120
bec3dpa_dre00120 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00120_genes)
bec7dpa_dre00120 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00120_genes)
hep3dpa_dre00120 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00120_genes)
hep7dpa_dre00120 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00120_genes)

dre00561_genes <- gene_symbols_list$dre00561
bec3dpa_dre00561 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00561_genes)
bec7dpa_dre00561 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00561_genes)
hep3dpa_dre00561 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00561_genes)
hep7dpa_dre00561 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00561_genes)

dre00564_genes <- gene_symbols_list$dre00564
bec3dpa_dre00564 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00564_genes)
bec7dpa_dre00564 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00564_genes)
hep3dpa_dre00564 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00564_genes)
hep7dpa_dre00564 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00564_genes)

dre00980_genes <- gene_symbols_list$dre00980
bec3dpa_dre00980 <- bec3dpa_df %>% filter(bec3dpa_df[[1]] %in% dre00980_genes)
bec7dpa_dre00980 <- bec7dpa_df %>% filter(bec7dpa_df[[1]] %in% dre00980_genes)
hep3dpa_dre00980 <- hep3dpa_df %>% filter(hep3dpa_df[[1]] %in% dre00980_genes)
hep7dpa_dre00980 <- hep7dpa_df %>% filter(hep7dpa_df[[1]] %in% dre00980_genes)

# define the function to filter for significant genes
filter_significant_genes <- function(data, gene_list, p_value_threshold = 0.05, fold_change_threshold = 1) {
  # filter the data for the genes in the gene list and for significant results
  filtered_data <- data %>%
    filter(data[[1]] %in% gene_list, PValue < p_value_threshold, abs(logFC) > fold_change_threshold)
  
  return(filtered_data)
}

# apply the function 
bec3dpa_dre00010_significant <- filter_significant_genes(bec3dpa_dre00010, dre00010_genes)
bec7dpa_dre00010_significant <- filter_significant_genes(bec7dpa_dre00010, dre00010_genes)
hep3dpa_dre00010_significant <- filter_significant_genes(hep3dpa_dre00010, dre00010_genes)
hep7dpa_dre00010_significant <- filter_significant_genes(hep7dpa_dre00010, dre00010_genes)

# Define the function to extract rows and columns based on significant genes and timepoint
extract_significant_genes_from_normalized <- function(significant_data, normalized_data, timepoint) {
  # Get the gene names from the first column of the significant data
  significant_genes <- significant_data[[1]]
  
  # Extract the rows from normalized_hep where row names match significant gene names
  matched_data <- normalized_data[rownames(normalized_data) %in% significant_genes, ]
  
  # Subset columns based on timepoint (3dpa or 7dpa)
  if (timepoint == "3dpa") {
    # For 3dpa, take columns 10-12 and 16-18
    matched_data <- matched_data[, c(10:12, 16:18)]
  } else if (timepoint == "7dpa") {
    # For 7dpa, take columns 13-15 and 16-18
    matched_data <- matched_data[, c(13:15, 16:18)]
  }
  
  return(as.data.frame(matched_data))
}

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")


