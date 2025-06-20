library(biomaRt)
library(CellChat)
library(Seurat)

# initally need to convert to human orthologs
mapping <- read_csv('/Users/sophiemarcotte/Desktop/mart_export (1).txt')
colnames(mapping) <- c("human_gene", "zebrafish_gene")

# remove duplicates or NAs
mapping <- mapping[!duplicated(mapping$zebrafish_gene), ]
mapping <- mapping[!is.na(mapping$human_gene), ]

# only keep genes that are in the seuratobj
mapping <- mapping[!is.na(mapping$human_gene) & mapping$zebrafish_gene %in% rownames(zf), ]

# get the expression matrix
expr_mat <- GetAssayData(zf, assay = 'RNA', slot = 'data')

# create a list of expression values for each gene
# allow z genes with more than one ortholog to be there twice
expr_list <- lapply(seq_len(nrow(mapping)), function(i) {
  zf_gene <- mapping$zebrafish_gene[i]
  hs_gene <- mapping$human_gene[i]
  
  vec <- expr_mat[zf_gene, , drop = FALSE]
  rownames(vec) <- hs_gene
  return(vec)
})

# put the lists back into a matrix
new_expr_mat <- do.call(rbind, expr_list)

# collapse duplicate human genes by summing expression - just be aware downstream
collapsed_expr <- rowsum(as.matrix(new_expr_mat), group = rownames(new_expr_mat))

zf[["RNA_human"]] <- CreateAssayObject(counts = collapsed_expr)
DefaultAssay(zf) <- "RNA_human"
Idents(zf) <- zf@meta.data$cell.type.12.long

zf <- NormalizeData(zf,scale.factor = 10000)

# PREPARE DATA TO RUN 
# pull out counts matrix
labels <- Idents(zf)
data.input <- zf[["RNA_human"]]@data
meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
