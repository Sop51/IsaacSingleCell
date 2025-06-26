library(biomaRt)
library(CellChat)
library(Seurat)
library(anndata)
library(zellkonverter)

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
# online data
ad <- read_h5ad("/Users/sm2949/Downloads/GSE272484_CellCousin_final.h5ad")
counts <- t(as.matrix(ad$X))
data.input <- normalizeData(counts)
meta <- ad$obs
meta$labels <- meta[["anno6"]]

# pull out counts matrix and prepare the metadata ----
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
cellchat <- computeCommunProb(cellchat, type = "triMean",
                              raw.use = TRUE) 

# filter the cell-cell communication based on the number of cells in each group (10 min) ----
cellchat <- filterCommunication(cellchat, min.cells = 10)

# infer cell-cell communication at the pathway level ----
cellchat <- computeCommunProbPathway(cellchat)

# calculate the aggregated cell-cell communication level (two cell types here) ----
cellchat <- aggregateNet(cellchat)
sources.use = c("Macrophages")
targets.use = c("Cholangiocytes-Control", "Endothelial", "HSCs")
cellchat <- aggregateNet(cellchat, sources.use = sources.use, targets.use = targets.use)

# save results as an rds ----
saveRDS(cellchat, file = "cellchat_GSE272484.rds")

# do some visualization ----
groupSize <- as.numeric(table(cellchat@idents))
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# look at the signaling pathways  ----
pathways.show.all <- cellchat@netP$pathways
pathways.show <- c("JAM")

# macrophage ----
netVisual_bubble(cellchat, sources.use = 21, targets.use = c(1:21), remove.isolate = FALSE)

# endothelial cells ----
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:11), remove.isolate = FALSE)

# hepatocyte cells ----
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:11), remove.isolate = FALSE)

# biliary epithelial cell ----
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:11), remove.isolate = FALSE)

# looking at a specific pathway
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL)

pairLR.NOTCH <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.NOTCH[2,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

genes.use <- extractEnrichedLR(cellchat, signaling = "JAM", geneLR.return = TRUE)$geneLR
Seurat::VlnPlot(zf, features = genes.use)

# identify signaling roles (dominant senders and recievers) ----
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pathways.show <- c("TIGIT")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# global communication patterns ----
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

# try some manual lookups
ligands <- read_excel('/Users/sm2949/Desktop/sigmeta.xlsx')
ligands_z <- left_join(ligands, mapping, by = c("Ligand" = "human_gene"))
mac_1_ligands <- sig_genes_t1[rownames(sig_genes_t1) %in% ligands_z$zebrafish_gene, ]
mac_2_ligands <- sig_genes_t2[rownames(sig_genes_t2) %in% ligands_z$zebrafish_gene, ]

# pull out receptors expressed by BECs @ 1+2pda
receptors <- read.csv('/Users/sm2949/Desktop/receptors.txt', sep='\t')
receptors_z <- mapping[mapping$human_gene %in% receptors$Hgnc.Symbol, ]
bil_1_receptors <- sig_genes_t1_bil[rownames(sig_genes_t1_bil) %in% receptors_z$zebrafish_gene, ]
bil_2_receptors <- sig_genes_t2_bil[rownames(sig_genes_t2_bil) %in% receptors_z$zebrafish_gene, ]

library(clusterProfiler)
library(org.Dr.eg.db)
library(enrichplot)
entrez_ids <- bitr(rownames(sig_genes_t2), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dr.eg.db")
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID,
                          organism = "dre", # zebrafish
                          pvalueCutoff = 0.05)
go_enrich <- enrichGO(gene = entrez_ids$ENTREZID,
                      OrgDb = org.Dr.eg.db,
                      ont = "BP", # Biological Process
                      pvalueCutoff = 0.05)
go_results <- as.data.frame(go_enrich@result)
