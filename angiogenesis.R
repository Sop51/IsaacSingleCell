library(KEGGREST)
library(biomaRt)
library(viridis)
library(dplyr)

# helper functions
# filter for significant genes
filter_significant_genes <- function(df, gene_list) {
  df <- df[df$PValue < 0.05 & df$logFC > 1.5, ] 
  df <- df[rownames(df) %in% unique_angio_genes, ]
}

# create a copy of the seurat object
zf_angio <- zf

# get the angiogenesis genes
# kegg
mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")
kegg_angio_genes <- keggGet("dre04370")[[1]]$GENE
kegg_angio_gene_ids <- unique(kegg_angio_genes[seq(1, length(kegg_angio_genes), by = 2)])
angio_gene_mapping <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                                 filters = "entrezgene_id",
                                 values = kegg_angio_gene_ids,
                                 mart = mart)

angiogenesis_genes <- angio_gene_mapping$external_gene_name

# go
go_angio_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006"),
  filters = "go",
  values = "GO:0001525",
  mart = mart
)

go_angio_genes <- unique(go_angio_genes$external_gene_name)

# get final list
unique_angio_genes <- unique(c(angiogenesis_genes, go_angio_genes))

# create a grouping variable
zf_angio$group <- paste0(zf_angio$cell.type.12.long, "_", zf_angio$timepoint)

# calculate average expression per cell type per timepoint
avg_expr <- AverageExpression(
  zf_angio,
  assays = 'RNA',
  features = unique_angio_genes,
  group.by = "group",
  return.seurat = FALSE,
  layer = 'data'
)

# pull into a dataframe and set gene as a column
avg_expr_df <- avg_expr[[1]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene")

# reshape into long format
avg_expr_long <- avg_expr_df %>%
  pivot_longer(-gene, names_to = "cell_condition", values_to = "expression")

# create cell type and condition (timepoint) column
avg_expr_long <- avg_expr_long %>%
  separate(cell_condition, into = c("cell_type", "condition"), sep = "-(?=[^+-]*$)", remove = FALSE)

# remove the untreated group and any introduced NA values
avg_expr_long <- avg_expr_long %>%
  filter(condition != "untreated", !is.na(expression))

# separate mock for joining
mock_expr <- avg_expr_long %>%
  filter(condition == "mock") %>%
  rename(mock_expression = expression) %>%
  dplyr::select(gene, cell_type, mock_expression)

# join mock expression to each gene/timepoint
logfc_df <- avg_expr_long %>%
  filter(condition != "mock") %>%
  left_join(mock_expr, by = c("gene", "cell_type")) %>%
  mutate(logFC = expression - mock_expression)

# count genes with logFC > 1.5
upregulated_counts <- logfc_df %>%
  filter(logFC > 1.5) %>%
  group_by(cell_type, condition) %>%
  summarise(num_up_genes = n(), .groups = "drop")

# create the plot
ggplot(upregulated_counts, aes(x = cell_type, y = num_up_genes, fill = condition)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  scale_fill_viridis(discrete = TRUE, option = "D", name = "Timepoint") +
  labs(
    title = "Upregulated Angiogenesis Genes (logFC > 1.5)",
    x = "Cell Type",
    y = "Number of Genes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
