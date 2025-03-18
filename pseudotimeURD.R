library(rdl)
library(Seurat)
library(URD)
library(devtools)

# subset the seurat obj to only ecm producing
cell_types_of_interest <- c("Macrophage")  
zf_subset <- subset(zf, subset = cell.type.12.long %in% cell_types_of_interest & timepoint != "untreated" & timepoint != "mock")

counts <- as.matrix(zf_subset@assays[['RNA']]@counts)
meta.data <- zf_subset@meta.data

# create an URD object, filters the data, normalize and log transform
urd_obj <- URD::createURD(count.data = counts, meta = meta.data, min.cells=3, min.counts=3)

# copy metadata information over
urd_obj@group.ids$timepoint <- as.character(urd_obj@meta[rownames(urd_obj@group.ids),"timepoint"])
urd_obj@group.ids$cell.type <- as.character(urd_obj@meta[rownames(urd_obj@group.ids),"cell.type.12.long"])

# extract the timepoints
timepoints <- sort(unique(urd_obj@group.ids$timepoint))

# loop over each timepoint to get variable genes
var.by.timepoint <- lapply(timepoints, function(timepoint) {
  # get cells for the specific timepoint
  cells_in_timepoint <- URD::cellsInCluster(urd_obj, "timepoint", timepoint)
  # find variable genes for those cells
  URD::findVariableGenes(
    urd_obj, 
    cells.fit = cells_in_timepoint, 
    set.object.var.genes = FALSE, 
    diffCV.cutoff = 0.3, 
    mean.min = 0.005, 
    mean.max = 100, 
    main.use = paste0("Timepoint ", timepoint), 
    do.plot = TRUE
  )
})

# combine the results from each timepoint into a single list of variable genes and load into the URD object
var.genes <- sort(unique(unlist(var.by.timepoint)))
urd_obj@var.genes <- var.genes

# calculate PCA and consider those PCs that with standard deviation 2x expected by noise as significant
urd_obj <- URD::calcPCA(urd_obj, mp.factor = 2)
URD::pcSDPlot(urd_obj)

# calculate tSNE
set.seed(19)
urd_obj <- URD::calcTsne(object = urd_obj)
URD::plotDim(urd_obj, "timepoint", plot.title = "tSNE: Timepoint")

# knn = sqrt(n.cells)) = sqrt(6276) = 79.22121, round since few cells
urd_obj <- URD::calcDM(urd_obj, knn = 100, sigma=20)
URD::plotDimArray(urd_obj, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma local, 119 NNs): Timepoint", label="timepoint", plot.title="", legend=F)
URD::plotDim(urd_obj, "timepoint", transitions.plot = 10000, plot.title="Timepoint (with transitions)")

# Use all cells from the first stage as the root
root.cells <- URD::cellsInCluster(urd_obj, "timepoint", "7dpa")

# Then we run 'flood' simulations
obj.floods <- URD::floodPseudotime(urd_obj, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

# The we process the simulations into a pseudotime
urd_obj <- URD::floodPseudotimeProcess(urd_obj, obj.floods, floods.name="pseudotime")

URD::pseudotimePlotStabilityOverall(urd_obj)
URD::plotDim(urd_obj, "pseudotime")

URD::plotDists(urd_obj, "pseudotime", "timepoint", plot.title="Pseudotime by timepoint")

# Create a subsetted object of just those cells from the final timepoint
urd_obj.0dpa <- URD::urdSubset(urd_obj, cells.keep=URD::cellsInCluster(urd_obj, "timepoint", "7dpa"))

# Use the variable genes that were calculated only on 7dpa
urd_obj.0dpa@var.genes <- var.by.timepoint[[1]]

# Calculate PCA and tSNE
urd_obj.0dpa <- URD::calcPCA(urd_obj.0dpa, mp.factor = 1.5)
URD::pcSDPlot(urd_obj.0dpa)

urd_obj.0dpa <- URD::calcTsne(urd_obj.0dpa)

# Calculate graph clustering of these cells
urd_obj.0dpa <- URD::graphClustering(urd_obj.0dpa, num.nn = 80, do.jaccard=T, method="Louvain")
URD::plotDim(urd_obj.0dpa, "Louvain-80", plot.title = "Louvain (50 NN) graph clustering", point.size=3)

# Copy cluster identities from axial.6somite object to a new clustering ("tip.clusters") in the full axial object.
urd_obj@group.ids[rownames(urd_obj.0dpa@group.ids), "tip.clusters"] <- urd_obj.0dpa@group.ids$`Louvain-80`

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
urd.ptlogistic <- URD::pseudotimeDetermineLogistic(urd_obj, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)

# Bias the transition matrix acording to pseudotime
urd.biased.tm <- as.matrix(URD::pseudotimeWeightTransitionMatrix(urd_obj, "pseudotime", logistic.params=urd.ptlogistic))

# Simulate the biased random walks from each tip
urd.walks <- URD::simulateRandomWalksFromTips(urd_obj, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = urd.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)

# Process the biased random walks into visitation frequencies
urd_obj <- URD::processRandomWalksFromTips(urd_obj, urd.walks, verbose = F)
URD::plotDim(urd_obj, "tip.clusters", plot.title="Cells in each tip")

# Load the cells used for each tip into the URD object
urd_obj.tree <- URD::loadTipCells(urd_obj, "tip.clusters")

# Build the tree
urd_obj.tree <- URD::buildTree(urd_obj.tree, pseudotime = "pseudotime", tips.use=1:2, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

URD::plotTree(urd_obj.tree, "timepoint", title="Timepoint")

# Generate the force-directed layout
urd_obj.tree <- URD::treeForceDirectedLayout(urd_obj.tree, num.nn=100, cut.unconnected.segments=2, verbose=T)

URD::plotTreeForce(urd_obj.tree, "anxa2a", title = "GSC", title.cex = 2, title.line=2.5)

