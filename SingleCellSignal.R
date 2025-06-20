library(biomaRt)
library(CellChat)
library(Seurat)

# ----------------- MAP SEURAT OBJ TO HUMAN ORTHOLOGS ---------------- #
# initally need to convert to human orthologs ----
mapping <- read_csv('/Users/sm2949/Desktop/mart_export.txt')
colnames(mapping) <- c("zebrafish_gene", "human_gene")

# remove duplicates or NAs ----
mapping <- mapping[!is.na(mapping$human_gene), ]

# only keep genes that are in the seurat obj ----
mapping <- mapping[mapping$zebrafish_gene %in% rownames(zf), ]

# get the expression matrix ----
expr_mat <- GetAssayData(zf, assay = 'RNA', slot = 'data')

# create a list of expression values for each gene ----
# allow z genes with more than one ortholog to be there twice
expr_list <- lapply(seq_len(nrow(mapping)), function(i) {
  zf_gene <- mapping$zebrafish_gene[i]
  hs_gene <- mapping$human_gene[i]
  
  vec <- expr_mat[zf_gene, , drop = FALSE]
  rownames(vec) <- hs_gene
  return(vec)
})

# put the lists back into a matrix ----
new_expr_mat <- do.call(rbind, expr_list)

# collapse duplicate human genes by summing expression - just be aware downstream ----
collapsed_expr <- rowsum(as.matrix(new_expr_mat), group = rownames(new_expr_mat))

zf[["RNA_human"]] <- CreateAssayObject(counts = collapsed_expr)
DefaultAssay(zf) <- "RNA_human"
Idents(zf) <- zf@meta.data$cell.type.12.long

# ------------------- PREPARE DATA TO RUN ------------------------ #
# pull out counts matrix and prepare the metadata ----
zf <- NormalizeData(zf)
labels <- Idents(zf)
data.input <- zf[["RNA"]]@data
meta <- data.frame(labels = labels, row.names = names(labels))

# create a cell chat object ----
cellchat <- createCellChat(object = data.input, group.by =
                             "labels", meta = meta)

# select the L-R interaction database ----
CellChatDB <- CellChatDB.zebrafish 

# select the whole database ----
CellChatDB.use <- CellChatDB

# set the selected database in the object then subset the exp matrix to relevant genes ----
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

# identify over expressed L and R from above in data ----
cellchat <- identifyOverExpressedGenes(cellchat)

# for each over expressed L or R, see if its associated L or R is over expressed as well ----
cellchat <- identifyOverExpressedInteractions(cellchat)

# smooth the data since our sequencing depth is lower ----
# cellchat <- smoothData(cellchat, adj = PPI.zebrafish) no zebrafish :(

# infer cell-cell communication at a L-R pair level ----
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05,
                              raw.use = TRUE) 

# filter the cell-cell communication based on the number of cells in each group (10 min) ----
cellchat <- filterCommunication(cellchat, min.cells = 10)

# infer cell-cell communication at the pathway level ----
cellchat <- computeCommunProbPathway(cellchat)

# calculate the aggregated cell-cell communication level (two cell types here) ----
#cellchat <- aggregateNet(cellchat)
sources.use = c("Macrophage")
targets.use = c("Biliary Epithelial Cell")
cellchat <- aggregateNet(cellchat, sources.use = sources.use, targets.use = targets.use)

# save results as an rds 
saveRDS(cellchat, file = "cellchat_mac_bec_z.rds")

# look at the signaling pathways 
pathways.show.all <- cellchat@netP$pathways
pathways.show <- c("MIF")

netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

pairLR.MIF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MIF[1,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
