library(KEGGREST)
library(dplyr)
library(pheatmap)
library(readr)
library(AnnotationDbi)
library(pheatmap)
library(clusterProfiler)
library(org.Dr.eg.db)

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

# define the function to filter for significant genes
filter_significant_genes <- function(data, p_value_threshold = 0.05, fold_change_threshold = 1) {
  # filter the data for the genes in the gene list and for significant results
  filtered_data <- data %>%
    filter(PValue < p_value_threshold, abs(logFC) > fold_change_threshold)
  
  return(filtered_data)
}

hep3dpa_df <- filter_significant_genes(hep3dpa_df)
hep7dpa_df <- filter_significant_genes(hep7dpa_df)
bec3dpa_df <- filter_significant_genes(bec3dpa_df)
bec7dpa_df <- filter_significant_genes(bec7dpa_df)

# -------------- KEGG / GO --------------------- #
# separate into up and down regulated genes
up_hep3dpa <- hep3dpa_df$...1[hep3dpa_df$logFC > 0]
down_hep3dpa <- hep3dpa_df$...1[hep3dpa_df$logFC < 0]

up_hep7dpa <- hep7dpa_df$...1[hep7dpa_df$logFC > 0]
down_hep7dpa <- hep7dpa_df$...1[hep7dpa_df$logFC < 0]

up_bec3dpa <- bec3dpa_df$...1[bec3dpa_df$logFC > 0]
down_bec3dpa <- bec3dpa_df$...1[bec3dpa_df$logFC < 0]

up_bec7dpa <- bec7dpa_df$...1[bec7dpa_df$logFC > 0]
down_bec7dpa <- bec7dpa_df$...1[bec7dpa_df$logFC < 0]

# GO enrichment for upregulated and downregulated genes for each dataset

# Hep3dpa
go_up_hep3dpa <- enrichGO(gene = up_hep3dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
go_down_hep3dpa <- enrichGO(gene = down_hep3dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

dotplot(go_up_hep3dpa) + ggtitle("GO Enrichment for Upregulated Genes in Hep3dpa")
dotplot(go_down_hep3dpa) + ggtitle("GO Enrichment for Downregulated Genes in Hep3dpa")

# Hep7dpa
go_up_hep7dpa <- enrichGO(gene = up_hep7dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
go_down_hep7dpa <- enrichGO(gene = down_hep7dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

dotplot(go_up_hep7dpa) + ggtitle("GO Enrichment for Upregulated Genes in Hep7dpa")
dotplot(go_down_hep7dpa) + ggtitle("GO Enrichment for Downregulated Genes in Hep7dpa")

# Bec3dpa
go_up_bec3dpa <- enrichGO(gene = up_bec3dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
go_down_bec3dpa <- enrichGO(gene = down_bec3dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

dotplot(go_up_bec3dpa) + ggtitle("GO Enrichment for Upregulated Genes in Bec3dpa")
dotplot(go_down_bec3dpa) + ggtitle("GO Enrichment for Downregulated Genes in Bec3dpa")

# Bec7dpa
go_up_bec7dpa <- enrichGO(gene = up_bec7dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
go_down_bec7dpa <- enrichGO(gene = down_bec7dpa, OrgDb = org.Dr.eg.db, keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)

dotplot(go_up_bec7dpa) + ggtitle("GO Enrichment for Upregulated Genes in Bec7dpa")
dotplot(go_down_bec7dpa) + ggtitle("GO Enrichment for Downregulated Genes in Bec7dpa")

# KEGG enrichment for upregulated and downregulated genes for each dataset

# Function to convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
  return(entrez_ids$ENTREZID)
}

# Example for Hep3dpa upregulated genes
up_hep3dpa_entrez <- convert_to_entrez(up_hep3dpa)

# Example for Hep3dpa downregulated genes
down_hep3dpa_entrez <- convert_to_entrez(down_hep3dpa)

# Repeat for Hep7dpa, Bec3dpa, Bec7dpa
up_hep7dpa_entrez <- convert_to_entrez(up_hep7dpa)
down_hep7dpa_entrez <- convert_to_entrez(down_hep7dpa)

up_bec3dpa_entrez <- convert_to_entrez(up_bec3dpa)
down_bec3dpa_entrez <- convert_to_entrez(down_bec3dpa)

up_bec7dpa_entrez <- convert_to_entrez(up_bec7dpa)
down_bec7dpa_entrez <- convert_to_entrez(down_bec7dpa)

# Hep3dpa KEGG enrichment
kegg_up_hep3dpa <- enrichKEGG(gene = up_hep3dpa_entrez, organism = "dre", pvalueCutoff = 0.05)
kegg_down_hep3dpa <- enrichKEGG(gene = down_hep3dpa_entrez, organism = "dre", pvalueCutoff = 0.05)

dotplot(kegg_up_hep3dpa) + ggtitle("KEGG Enrichment for Upregulated Genes in Hep3dpa")
dotplot(kegg_down_hep3dpa) + ggtitle("KEGG Enrichment for Downregulated Genes in Hep3dpa")

# Hep7dpa KEGG enrichment
kegg_up_hep7dpa <- enrichKEGG(gene = up_hep7dpa_entrez, organism = "dre", pvalueCutoff = 0.05)
kegg_down_hep7dpa <- enrichKEGG(gene = down_hep7dpa_entrez, organism = "dre", pvalueCutoff = 0.05)

dotplot(kegg_up_hep7dpa) + ggtitle("KEGG Enrichment for Upregulated Genes in Hep7dpa")
dotplot(kegg_down_hep7dpa) + ggtitle("KEGG Enrichment for Downregulated Genes in Hep7dpa")

# Bec3dpa KEGG enrichment
kegg_up_bec3dpa <- enrichKEGG(gene = up_bec3dpa_entrez, organism = "dre", pvalueCutoff = 0.05)
kegg_down_bec3dpa <- enrichKEGG(gene = down_bec3dpa_entrez, organism = "dre", pvalueCutoff = 0.05)

dotplot(kegg_up_bec3dpa) + ggtitle("KEGG Enrichment for Upregulated Genes in Bec3dpa")
dotplot(kegg_down_bec3dpa) + ggtitle("KEGG Enrichment for Downregulated Genes in Bec3dpa")

# Bec7dpa KEGG enrichment
kegg_up_bec7dpa <- enrichKEGG(gene = up_bec7dpa_entrez, organism = "dre", pvalueCutoff = 0.05)
kegg_down_bec7dpa <- enrichKEGG(gene = down_bec7dpa_entrez, organism = "dre", pvalueCutoff = 0.05)

dotplot(kegg_up_bec7dpa) + ggtitle("KEGG Enrichment for Upregulated Genes in Bec7dpa")
dotplot(kegg_down_bec7dpa) + ggtitle("KEGG Enrichment for Downregulated Genes in Bec7dpa")


# ------------ Metabolism Genes data manipulation ----------------- #
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

<<<<<<< HEAD
# define the function to filter for significant genes
filter_significant_genes <- function(data, gene_list, p_value_threshold = 0.05, fold_change_threshold = 1) {
  # filter the data for the genes in the gene list and for significant results
  filtered_data <- data %>%
    filter(data[[1]] %in% gene_list, PValue < p_value_threshold, abs(logFC) > fold_change_threshold)
  
  return(filtered_data)
}

# Apply the function for each pathway ID
bec3dpa_dre00010_significant <- filter_significant_genes(bec3dpa_dre00010, dre00010_genes)
bec7dpa_dre00010_significant <- filter_significant_genes(bec7dpa_dre00010, dre00010_genes)
hep3dpa_dre00010_significant <- filter_significant_genes(hep3dpa_dre00010, dre00010_genes)
hep7dpa_dre00010_significant <- filter_significant_genes(hep7dpa_dre00010, dre00010_genes)

bec3dpa_dre00071_significant <- filter_significant_genes(bec3dpa_dre00071, dre00071_genes)
bec7dpa_dre00071_significant <- filter_significant_genes(bec7dpa_dre00071, dre00071_genes)
hep3dpa_dre00071_significant <- filter_significant_genes(hep3dpa_dre00071, dre00071_genes)
hep7dpa_dre00071_significant <- filter_significant_genes(hep7dpa_dre00071, dre00071_genes)

bec3dpa_dre00061_significant <- filter_significant_genes(bec3dpa_dre000611, dre00061_genes)
bec7dpa_dre00061_significant <- filter_significant_genes(bec7dpa_dre00061, dre00061_genes)
hep3dpa_dre00061_significant <- filter_significant_genes(hep3dpa_dre00061, dre00061_genes)
hep7dpa_dre00061_significant <- filter_significant_genes(hep7dpa_dre00061, dre00061_genes)

bec3dpa_dre00020_significant <- filter_significant_genes(bec3dpa_dre00020, dre00020_genes)
bec7dpa_dre00020_significant <- filter_significant_genes(bec7dpa_dre00020, dre00020_genes)
hep3dpa_dre00020_significant <- filter_significant_genes(hep3dpa_dre00020, dre00020_genes)
hep7dpa_dre00020_significant <- filter_significant_genes(hep7dpa_dre00020, dre00020_genes)

bec3dpa_dre00120_significant <- filter_significant_genes(bec3dpa_dre00120, dre00120_genes)
bec7dpa_dre00120_significant <- filter_significant_genes(bec7dpa_dre00120, dre00120_genes)
hep3dpa_dre00120_significant <- filter_significant_genes(hep3dpa_dre00120, dre00120_genes)
hep7dpa_dre00120_significant <- filter_significant_genes(hep7dpa_dre00120, dre00120_genes)

bec3dpa_dre00561_significant <- filter_significant_genes(bec3dpa_dre00561, dre00561_genes)
bec7dpa_dre00561_significant <- filter_significant_genes(bec7dpa_dre00561, dre00561_genes)
hep3dpa_dre00561_significant <- filter_significant_genes(hep3dpa_dre00561, dre00561_genes)
hep7dpa_dre00561_significant <- filter_significant_genes(hep7dpa_dre00561, dre00561_genes)

bec3dpa_dre00564_significant <- filter_significant_genes(bec3dpa_dre00564, dre00564_genes)
bec7dpa_dre00564_significant <- filter_significant_genes(bec7dpa_dre00564, dre00564_genes)
hep3dpa_dre00564_significant <- filter_significant_genes(hep3dpa_dre00564, dre00564_genes)
hep7dpa_dre00564_significant <- filter_significant_genes(hep7dpa_dre00564, dre00564_genes)

bec3dpa_dre00980_significant <- filter_significant_genes(bec3dpa_dre00980, dre00980_genes)
bec7dpa_dre00980_significant <- filter_significant_genes(bec7dpa_dre00980, dre00980_genes)
hep3dpa_dre00980_significant <- filter_significant_genes(hep3dpa_dre00980, dre00980_genes)
hep7dpa_dre00980_significant <- filter_significant_genes(hep7dpa_dre00980, dre00980_genes)


=======
>>>>>>> 39b2c21076be93d90d183adb0669351a7afb3d81
# Define the function to extract rows and columns based on significant genes and timepoint
extract_significant_genes_from_normalized <- function(significant_data, normalized_data, timepoint) {
  # Get the gene names from the first column of the significant data
  significant_genes <- significant_data[[1]]
  
  # Extract the rows from normalized_hep where row names match significant gene names
  matched_data <- normalized_data[rownames(normalized_data) %in% significant_genes, , drop = FALSE]
  
  # Subset columns based on timepoint (3dpa or 7dpa)
  if (timepoint == "3dpa") {
    # For 3dpa, take columns 10-12 and 16-18
    matched_data <- matched_data[, c(10:12, 16:18), drop = FALSE]
  } else if (timepoint == "7dpa") {
    # For 7dpa, take columns 13-15 and 16-18
    matched_data <- matched_data[, c(13:15, 16:18), drop = FALSE]
  }
  
  return(as.data.frame(matched_data))
}

# define pathways of interest
pathways <- c("dre00010", "dre00071", "dre00061", "dre00020", "dre00120", 
              "dre00561", "dre00564", "dre00980")

<<<<<<< HEAD
bec3dpa_dre00071_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00071_significant, normalized_hep, "3dpa")
bec7dpa_dre00071_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00071_significant, normalized_hep, "7dpa")
hep3dpa_dre00071_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00071_significant, normalized_hep, "3dpa")
hep7dpa_dre00071_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00071_significant, normalized_hep, "7dpa")

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")

bec3dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec3dpa_dre00010_significant, normalized_hep, "3dpa")
bec7dpa_dre00010_matched <- extract_significant_genes_from_normalized(bec7dpa_dre00010_significant, normalized_hep, "7dpa")
hep3dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep3dpa_dre00010_significant, normalized_hep, "3dpa")
hep7dpa_dre00010_matched <- extract_significant_genes_from_normalized(hep7dpa_dre00010_significant, normalized_hep, "7dpa")
=======
pathway_names <- c("Glycolysis / Gluconeogenesis", "Fatty acid degradation", "Fatty acid biosynthesis",
                   "Citrate cycle (TCA cycle)", "Primary bile acid biosynthesis", "Glycerolipid metabolism",
                   "Glycerophospholipid metabolism", "Metabolism of xenobiotics by cytochrome P450")

# ---------------------- run for bec's at 3dpa for each pathway -------------------#
bec3dpa_dre00010_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00010, normalized_bec, "3dpa")
bec3dpa_dre00071_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00071, normalized_bec, "3dpa")
bec3dpa_dre00061_counts <- extract_significant_genes_from_normalized(bec3dpa_dre000611, normalized_bec, "3dpa")
bec3dpa_dre00020_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00020, normalized_bec, "3dpa")
bec3dpa_dre00120_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00120, normalized_bec, "3dpa")
bec3dpa_dre00561_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00561, normalized_bec, "3dpa")
bec3dpa_dre00564_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00564, normalized_bec, "3dpa")
bec3dpa_dre00980_counts <- extract_significant_genes_from_normalized(bec3dpa_dre00980, normalized_bec, "3dpa")

# Create a mapping from pathway ID to pathway name
pathway_mapping <- setNames(pathway_names, pathways)

# Put your count data frames into a named list
counts_list <- list(
  "dre00010" = bec3dpa_dre00010_counts,
  "dre00071" = bec3dpa_dre00071_counts,
  "dre00061" = bec3dpa_dre00061_counts,  
  "dre00020" = bec3dpa_dre00020_counts,
  "dre00120" = bec3dpa_dre00120_counts,
  "dre00561" = bec3dpa_dre00561_counts,
  "dre00564" = bec3dpa_dre00564_counts,
  "dre00980" = bec3dpa_dre00980_counts
)

# loop through each data frame in the list and add Gene and Pathway columns
for(pid in names(counts_list)) {
  df <- counts_list[[pid]]
  # Add a column for gene names, assuming rownames contain the gene names
  df$Gene <- rownames(df)
  # Add a column for the pathway name (using our mapping)
  df$Pathway <- pathway_mapping[pid]
  # Save back to the list
  counts_list[[pid]] <- df
}

# combine all the data frames into one
combined_df <- do.call(rbind, counts_list)

# remove duplicates: keep only the first occurrence of each gene
combined_df_unique <- combined_df %>% distinct(Gene, .keep_all = TRUE)

# set the row names to the Gene column
rownames(combined_df_unique) <- combined_df_unique$Gene

# Define the new names manually based on the indices
new_colnames <- c("3dpa.1", "3dpa.2", "3dpa.3", "mock.1", "mock.2", "mock.3", "Gene", "Pathway") 

# Assign these names to your dataframe
colnames(combined_df_unique) <- new_colnames

# Remove non-numeric columns and convert to matrix for heatmap
heatmap_data <- combined_df_unique %>% dplyr::select(-Gene, -Pathway)
heatmap_matrix <- as.matrix(heatmap_data)

# Create row annotation for pathways
row_annotation <- data.frame(Pathway = combined_df_unique$Pathway)
rownames(row_annotation) <- combined_df_unique$Gene

# Define colors for pathways
pathway_colors <- RColorBrewer::brewer.pal(length(unique(row_annotation$Pathway)), "Set3")
names(pathway_colors) <- unique(row_annotation$Pathway)
annotation_colors <- list(Pathway = pathway_colors)

# Generate heatmap
pheatmap(
  heatmap_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  scale = "row",
  fontsize_row = 6,  
  fontsize_col = 10,
  main = "Pathway-Specific Differential Expression of \n Liver Metabolism and Functional Genes in BECs"
)

# ------------------ run for bec's at 7dpa for each pathway ------------------------ #
bec7dpa_dre00010_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00010, normalized_bec, "7dpa")
bec7dpa_dre00071_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00071, normalized_bec, "7dpa")
bec7dpa_dre00061_counts <- extract_significant_genes_from_normalized(bec7dpa_dre000611, normalized_bec, "7dpa")
bec7dpa_dre00020_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00020, normalized_bec, "7dpa")
bec7dpa_dre00120_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00120, normalized_bec, "7dpa")
bec7dpa_dre00561_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00561, normalized_bec, "7dpa")
bec7dpa_dre00564_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00564, normalized_bec, "7dpa")
bec7dpa_dre00980_counts <- extract_significant_genes_from_normalized(bec7dpa_dre00980, normalized_bec, "7dpa")

# Create a mapping from pathway ID to pathway name
pathway_mapping <- setNames(pathway_names, pathways)

# Put your count data frames into a named list
counts_list <- list(
  "dre00010" = bec7dpa_dre00010_counts,
  "dre00071" = bec7dpa_dre00071_counts,
  #"dre00061" = bec7dpa_dre00061_counts,  
  "dre00020" = bec7dpa_dre00020_counts,
  "dre00120" = bec7dpa_dre00120_counts,
  "dre00561" = bec7dpa_dre00561_counts,
  "dre00564" = bec7dpa_dre00564_counts,
  "dre00980" = bec7dpa_dre00980_counts
)

# loop through each data frame in the list and add Gene and Pathway columns
for(pid in names(counts_list)) {
  df <- counts_list[[pid]]
  # Add a column for gene names, assuming rownames contain the gene names
  df$Gene <- rownames(df)
  # Add a column for the pathway name (using our mapping)
  df$Pathway <- pathway_mapping[pid]
  # Save back to the list
  counts_list[[pid]] <- df
}

# combine all the data frames into one
combined_df <- do.call(rbind, counts_list)

# remove duplicates: keep only the first occurrence of each gene
combined_df_unique <- combined_df %>% distinct(Gene, .keep_all = TRUE)

# set the row names to the Gene column
rownames(combined_df_unique) <- combined_df_unique$Gene

# Define the new names manually based on the indices
new_colnames <- c("7dpa.1", "7dpa.2", "7dpa.3", "mock.1", "mock.2", "mock.3", "Gene", "Pathway") 

# Assign these names to your dataframe
colnames(combined_df_unique) <- new_colnames

# Remove non-numeric columns and convert to matrix for heatmap
heatmap_data <- combined_df_unique %>% dplyr::select(-Gene, -Pathway)
heatmap_matrix <- as.matrix(heatmap_data)

# Create row annotation for pathways
row_annotation <- data.frame(Pathway = combined_df_unique$Pathway)
rownames(row_annotation) <- combined_df_unique$Gene

# Define colors for pathways
pathway_colors <- RColorBrewer::brewer.pal(length(unique(row_annotation$Pathway)), "Set3")
names(pathway_colors) <- unique(row_annotation$Pathway)
annotation_colors <- list(Pathway = pathway_colors)

# Generate heatmap
pheatmap(
  heatmap_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  scale = "row",
  fontsize_row = 6,  
  fontsize_col = 10,
  main = "Pathway-Specific Differential Expression of \n Liver Metabolism and Functional Genes in BECs"
)

# ---------------------- run for heps at 3dpa for each pathway -------------------#
hep3dpa_dre00010_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00010, normalized_hep, "3dpa")
hep3dpa_dre00071_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00071, normalized_hep, "3dpa")
hep3dpa_dre00061_counts <- extract_significant_genes_from_normalized(hep3dpa_dre000611, normalized_hep, "3dpa")
hep3dpa_dre00020_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00020, normalized_hep, "3dpa")
hep3dpa_dre00120_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00120, normalized_hep, "3dpa")
hep3dpa_dre00561_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00561, normalized_hep, "3dpa")
hep3dpa_dre00564_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00564, normalized_hep, "3dpa")
hep3dpa_dre00980_counts <- extract_significant_genes_from_normalized(hep3dpa_dre00980, normalized_hep, "3dpa")

# Create a mapping from pathway ID to pathway name
pathway_mapping <- setNames(pathway_names, pathways)

# Put your count data frames into a named list
counts_list <- list(
  "dre00010" = hep3dpa_dre00010_counts,
  "dre00071" = hep3dpa_dre00071_counts,
  "dre00061" = hep3dpa_dre00061_counts,  
  "dre00020" = hep3dpa_dre00020_counts,
  "dre00120" = hep3dpa_dre00120_counts,
  "dre00561" = hep3dpa_dre00561_counts,
  "dre00564" = hep3dpa_dre00564_counts,
  "dre00980" = hep3dpa_dre00980_counts
)

# loop through each data frame in the list and add Gene and Pathway columns
for(pid in names(counts_list)) {
  df <- counts_list[[pid]]
  # Add a column for gene names, assuming rownames contain the gene names
  df$Gene <- rownames(df)
  # Add a column for the pathway name (using our mapping)
  df$Pathway <- pathway_mapping[pid]
  # Save back to the list
  counts_list[[pid]] <- df
}

# combine all the data frames into one
combined_df <- do.call(rbind, counts_list)

# remove duplicates: keep only the first occurrence of each gene
combined_df_unique <- combined_df %>% distinct(Gene, .keep_all = TRUE)

# set the row names to the Gene column
rownames(combined_df_unique) <- combined_df_unique$Gene

# Define the new names manually based on the indices
new_colnames <- c("3dpa.1", "3dpa.2", "3dpa.3", "mock.1", "mock.2", "mock.3", "Gene", "Pathway") 

# Assign these names to your dataframe
colnames(combined_df_unique) <- new_colnames

# Remove non-numeric columns and convert to matrix for heatmap
heatmap_data <- combined_df_unique %>% dplyr::select(-Gene, -Pathway)
heatmap_matrix <- as.matrix(heatmap_data)

# Create row annotation for pathways
row_annotation <- data.frame(Pathway = combined_df_unique$Pathway)
rownames(row_annotation) <- combined_df_unique$Gene

# Define colors for pathways
pathway_colors <- RColorBrewer::brewer.pal(length(unique(row_annotation$Pathway)), "Set3")
names(pathway_colors) <- unique(row_annotation$Pathway)
annotation_colors <- list(Pathway = pathway_colors)

# Generate heatmap
pheatmap(
  heatmap_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  scale = "row",
  fontsize_row = 6,  
  fontsize_col = 10,
  main = "Pathway-Specific Differential Expression of Liver \n Metabolism and Functional Genes in Hepatocytes"
)

# ------------------ run for hep at 7dpa for each pathway ------------------------ #
hep7dpa_dre00010_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00010, normalized_hep, "7dpa")
hep7dpa_dre00071_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00071, normalized_hep, "7dpa")
hep7dpa_dre00061_counts <- extract_significant_genes_from_normalized(hep7dpa_dre000611, normalized_hep, "7dpa")
hep7dpa_dre00020_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00020, normalized_hep, "7dpa")
hep7dpa_dre00120_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00120, normalized_hep, "7dpa")
hep7dpa_dre00561_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00561, normalized_hep, "7dpa")
hep7dpa_dre00564_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00564, normalized_hep, "7dpa")
hep7dpa_dre00980_counts <- extract_significant_genes_from_normalized(hep7dpa_dre00980, normalized_hep, "7dpa")

# Create a mapping from pathway ID to pathway name
pathway_mapping <- setNames(pathway_names, pathways)

# Put your count data frames into a named list
counts_list <- list(
  "dre00010" = hep7dpa_dre00010_counts,
  "dre00071" = hep7dpa_dre00071_counts,
  #"dre00061" = hep7dpa_dre00061_counts,  
  "dre00020" = hep7dpa_dre00020_counts,
  #"dre00120" = hep7dpa_dre00120_counts,
  "dre00561" = hep7dpa_dre00561_counts,
  "dre00564" = hep7dpa_dre00564_counts,
  "dre00980" = hep7dpa_dre00980_counts
)

# loop through each data frame in the list and add Gene and Pathway columns
for(pid in names(counts_list)) {
  df <- counts_list[[pid]]
  # Add a column for gene names, assuming rownames contain the gene names
  df$Gene <- rownames(df)
  # Add a column for the pathway name (using our mapping)
  df$Pathway <- pathway_mapping[pid]
  # Save back to the list
  counts_list[[pid]] <- df
}

# combine all the data frames into one
combined_df <- do.call(rbind, counts_list)

# remove duplicates: keep only the first occurrence of each gene
combined_df_unique <- combined_df %>% distinct(Gene, .keep_all = TRUE)

# set the row names to the Gene column
rownames(combined_df_unique) <- combined_df_unique$Gene

# Define the new names manually based on the indices
new_colnames <- c("7dpa.1", "7dpa.2", "7dpa.3", "mock.1", "mock.2", "mock.3", "Gene", "Pathway") 

# Assign these names to your dataframe
colnames(combined_df_unique) <- new_colnames

# Remove non-numeric columns and convert to matrix for heatmap
heatmap_data <- combined_df_unique %>% dplyr::select(-Gene, -Pathway)
heatmap_matrix <- as.matrix(heatmap_data)

# Create row annotation for pathways
row_annotation <- data.frame(Pathway = combined_df_unique$Pathway)
rownames(row_annotation) <- combined_df_unique$Gene

# Define colors for pathways
pathway_colors <- RColorBrewer::brewer.pal(length(unique(row_annotation$Pathway)), "Set3")
names(pathway_colors) <- unique(row_annotation$Pathway)
annotation_colors <- list(Pathway = pathway_colors)

# Generate heatmap
pheatmap(
  heatmap_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  scale = "row",
  fontsize_row = 6,  
  fontsize_col = 10,
  main = "Pathway-Specific Differential Expression of \n Liver Metabolism and Functional Genes in Hepatocytes"
)
>>>>>>> 39b2c21076be93d90d183adb0669351a7afb3d81
