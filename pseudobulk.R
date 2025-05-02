library(Seurat)
library(glmGamPoi)
library(tidyverse)
library(edgeR)
library(dplyr)
library(tidyverse)

# ------------------- pseudobulk set up ----------------------- #
# make a copy fo the original Seurat obj
zf_pseudo <- zf

# create a new meta data col with the time point and library
zf_pseudo$samples <- paste0(zf_pseudo$timepoint, zf_pseudo$orig.ident)

# aggregate expression per sample per cell type
cts <- AggregateExpression(zf_pseudo, 
                    group.by = c("cell.type.12.short", "samples"),
                    assays = 'RNA',
                    slot = 'counts',
                    return.seurat = FALSE)

# convert to a data frame
cts <- as.data.frame(cts$RNA)

# transpose the data frame
cts.t <- t(cts)

# get values where to split
splitrows <- gsub('_.*', '', rownames(cts.t))

# split the data frame
cts.split <- split.data.frame(cts.t,
                 f = factor(splitrows))

# fix the columns and transpose
cts.split.modified <- lapply(cts.split, function(x){
  # for each element remove the cell type name
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

# ------------------- running deseq analysis for bec ------------ #
# get the counts matrix
counts_bil <- cts.split.modified$Bil

# generate sample level metadata
colData <- data.frame(samples = colnames(counts_bil))

# create the metadata
colData <- colData %>%
  mutate(condition = ifelse(grepl("untreated", samples), "untreated",
                            ifelse(grepl("0dpa", samples), "0dpa",
                            ifelse(grepl("1dpa", samples), "1dpa",
                            ifelse(grepl("2dpa", samples), "2dpa",
                            ifelse(grepl("3dpa", samples), "3dpa",
                            ifelse(grepl("mock", samples), "mock",
                            ifelse(grepl("7dpa", samples), "7dpa", NA)))))))) %>% column_to_rownames(var = 'samples')

# set the condition to a factor vairable
groups <- as.factor(colData$condition)

# create the DGElist object
y <- DGEList(counts=counts_bil,group=groups)
# filter by expression
keep <- filterByExpr(y)
# filter by library size
y <- y[keep,,keep.lib.sizes=FALSE]
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_bec <- cpm(y)
# create the model design
design <- model.matrix(~0+groups)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

dpa0_bil_qlf <- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0))
dpa0_bil <- as.data.frame(dpa0_bil_qlf)

dpa1_bil_qlf <- glmQLFTest(fit, contrast = c(0,1,0,0,0,-1,0))
dpa1_bil <- as.data.frame(dpa1_bil_qlf)

dpa2_bil_qlf <- glmQLFTest(fit, contrast = c(0,0,1,0,0,-1,0))
dpa2_bil <- as.data.frame(dpa2_bil_qlf)

dpa3_bil_qlf <- glmQLFTest(fit, contrast = c(0,0,0,1,0,-1,0))
dpa3_bil <- as.data.frame(dpa3_bil_qlf)

dpa7_bil_qlf <- glmQLFTest(fit, contrast = c(0,0,0,0,1,-1,0))
dpa7_bil <- as.data.frame(dpa7_bil_qlf)

# run special comparisons for progenitor BEC detection
dpa1vsdpa2 <- glmQLFTest(fit, contrast = c(0,-1,1,0,0,0,0))
dpa1vsdpa2 <- as.data.frame(dpa1vsdpa2)
dpa2vsdpa3 <- glmQLFTest(fit, contrast = c(0,0,1,-1,0,0,0))
dpa2vsdpa3 <- as.data.frame(dpa2vsdpa3)

# save to a csv
write.csv(dpa1vsdpa2, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa1vsdpa2_bil_DE_results.csv", row.names = TRUE)
write.csv(dpa2vsdpa3, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa2vsdpa3_bil_DE_results.csv", row.names = TRUE)

# save each DE results data frame to CSV
write.csv(dpa0_bil, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa0_bil_DE_results.csv", row.names = TRUE)
write.csv(dpa1_bil, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa1_bil_DE_results.csv", row.names = TRUE)
write.csv(dpa2_bil, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa2_bil_DE_results.csv", row.names = TRUE)
write.csv(dpa3_bil, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa3_bil_DE_results.csv", row.names = TRUE)
write.csv(dpa7_bil, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa7_bil_DE_results.csv", row.names = TRUE)
# ------------------- pseudobulk for hepatocytes ----------------------- #
# get the counts matrix
counts_hep <- cts.split.modified$Hep

# generate sample level metadata
colData <- data.frame(samples = colnames(counts_hep))

# create the metadata
colData <- colData %>%
  mutate(condition = ifelse(grepl("untreated", samples), "untreated",
                     ifelse(grepl("0dpa", samples), "0dpa",
                     ifelse(grepl("1dpa", samples), "1dpa",
                     ifelse(grepl("2dpa", samples), "2dpa",
                     ifelse(grepl("3dpa", samples), "3dpa",
                     ifelse(grepl("mock", samples), "mock",
                     ifelse(grepl("7dpa", samples), "7dpa", NA)))))))) %>% column_to_rownames(var = 'samples')

# set the condition to a factor vairable
groups <- as.factor(colData$condition)

# create the DGElist object
y <- DGEList(counts=counts_hep,group=groups)
# filter by expression
keep <- filterByExpr(y)
# filter by library size
y <- y[keep,,keep.lib.sizes=FALSE]
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_hep <- cpm(y)
# create the model design
design <- model.matrix(~0+groups)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

dpa0_hep_qlf <- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0))
dpa0_hep <- as.data.frame(dpa0_hep_qlf)

dpa1_hep_qlf <- glmQLFTest(fit, contrast = c(0,1,0,0,0,-1,0))
dpa1_hep <- as.data.frame(dpa1_hep_qlf)

dpa2_hep_qlf <- glmQLFTest(fit, contrast = c(0,0,1,0,0,-1,0))
dpa2_hep <- as.data.frame(dpa2_hep_qlf)

dpa3_hep_qlf <- glmQLFTest(fit, contrast = c(0,0,0,1,0,-1,0))
dpa3_hep <- as.data.frame(dpa3_hep_qlf)

dpa7_hep_qlf <- glmQLFTest(fit, contrast = c(0,0,0,0,1,-1,0))
dpa7_hep <- as.data.frame(dpa7_hep_qlf)

# save each DE results data frame to CSV
write.csv(dpa0_hep, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa0_hep_DE_results.csv", row.names = TRUE)
write.csv(dpa1_hep, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa1_hep_DE_results.csv", row.names = TRUE)
write.csv(dpa2_hep, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa2_hep_DE_results.csv", row.names = TRUE)
write.csv(dpa3_hep, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa3_hep_DE_results.csv", row.names = TRUE)
write.csv(dpa7_hep, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa7_hep_DE_results.csv", row.names = TRUE)

# run special comparisons for newly formed hepatocytes
dpa2vsdpa3 <- glmQLFTest(fit, contrast = c(0,0,1,-1,0,0,0))
dpa2vsdpa3 <- as.data.frame(dpa2vsdpa3)
dpa3vsdpa7 <- glmQLFTest(fit, contrast = c(0,0,0,1,-1,0,0))
dpa3vsdpa7 <- as.data.frame(dpa3vsdpa7)

# save to a csv
write.csv(dpa2vsdpa3, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa2vsdpa3_hep_DE_results.csv", row.names = TRUE)
write.csv(dpa3vsdpa7, "/Users/sophiemarcotte/Desktop/SingleCellV2WithinClusterDE/dpa3vsdpa7_hep_DE_results.csv", row.names = TRUE)

# ------------------- pseudobulk for macrophages ----------------------- #
# get the counts matrix
counts_mac <- cts.split.modified$Mac

# generate sample level metadata
colData <- data.frame(samples = colnames(counts_mac))

# create the metadata
colData <- colData %>%
  mutate(condition = ifelse(grepl("untreated", samples), "untreated",
                     ifelse(grepl("0dpa", samples), "0dpa",
                     ifelse(grepl("1dpa", samples), "1dpa",
                     ifelse(grepl("2dpa", samples), "2dpa",
                     ifelse(grepl("3dpa", samples), "3dpa",
                     ifelse(grepl("mock", samples), "mock",
                     ifelse(grepl("7dpa", samples), "7dpa", NA)))))))) %>% column_to_rownames(var = 'samples')

# set the condition to a factor vairable
groups <- as.factor(colData$condition)

# create the DGElist object
y <- DGEList(counts=counts_mac,group=groups)
# filter by expression
keep <- filterByExpr(y)
# filter by library size
y <- y[keep,,keep.lib.sizes=FALSE]
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_mac <- cpm(y)
# create the model design
design <- model.matrix(~0+groups)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

dpa0_mac_qlf <- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0))
dpa0_mac <- as.data.frame(dpa0_mac_qlf)

dpa1_mac_qlf <- glmQLFTest(fit, contrast = c(0,1,0,0,0,-1,0))
dpa1_mac <- as.data.frame(dpa1_mac_qlf)

dpa2_mac_qlf <- glmQLFTest(fit, contrast = c(0,0,1,0,0,-1,0))
dpa2_mac <- as.data.frame(dpa2_mac_qlf)

dpa3_mac_qlf <- glmQLFTest(fit, contrast = c(0,0,0,1,0,-1,0))
dpa3_mac <- as.data.frame(dpa3_mac_qlf)

dpa7_mac_qlf <- glmQLFTest(fit, contrast = c(0,0,0,0,1,-1,0))
dpa7_mac <- as.data.frame(dpa7_mac_qlf)

# save each DE results data frame to CSV
write.csv(dpa0_mac, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_mac_DE_results.csv", row.names = TRUE)
write.csv(dpa1_mac, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_mac_DE_results.csv", row.names = TRUE)
write.csv(dpa2_mac, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_mac_DE_results.csv", row.names = TRUE)
write.csv(dpa3_mac, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_mac_DE_results.csv", row.names = TRUE)
write.csv(dpa7_mac, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_mac_DE_results.csv", row.names = TRUE)

# ------------------- pseudobulk for apln+ ----------------------- #
# get the counts matrix
counts_apln <- cts.split.modified$`apln+`

# generate sample level metadata
colData <- data.frame(samples = colnames(counts_apln))

# create the metadata
colData <- colData %>%
  mutate(condition = ifelse(grepl("untreated", samples), "untreated",
                      ifelse(grepl("0dpa", samples), "0dpa",
                      ifelse(grepl("1dpa", samples), "1dpa",
                      ifelse(grepl("2dpa", samples), "2dpa",
                      ifelse(grepl("3dpa", samples), "3dpa",
                      ifelse(grepl("mock", samples), "mock",
                      ifelse(grepl("7dpa", samples), "7dpa", NA)))))))) %>% column_to_rownames(var = 'samples')

# set the condition to a factor vairable
groups <- as.factor(colData$condition)

# create the DGElist object
y <- DGEList(counts=counts_apln,group=groups)
# filter by expression
keep <- filterByExpr(y)
# filter by library size
y <- y[keep,,keep.lib.sizes=FALSE]
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_apln <- cpm(y)
# create the model design
design <- model.matrix(~0+groups)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

dpa0_apln_qlf <- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0))
dpa0_apln <- as.data.frame(dpa0_apln_qlf)

dpa1_apln_qlf <- glmQLFTest(fit, contrast = c(0,1,0,0,0,-1,0))
dpa1_apln <- as.data.frame(dpa1_apln_qlf)

dpa2_apln_qlf <- glmQLFTest(fit, contrast = c(0,0,1,0,0,-1,0))
dpa2_apln <- as.data.frame(dpa2_apln_qlf)

dpa3_apln_qlf <- glmQLFTest(fit, contrast = c(0,0,0,1,0,-1,0))
dpa3_apln <- as.data.frame(dpa3_apln_qlf)

dpa7_apln_qlf <- glmQLFTest(fit, contrast = c(0,0,0,0,1,-1,0))
dpa7_apln <- as.data.frame(dpa7_apln_qlf)

# save each DE results data frame to CSV
write.csv(dpa0_apln, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_apln_DE_results.csv", row.names = TRUE)
write.csv(dpa1_apln, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_apln_DE_results.csv", row.names = TRUE)
write.csv(dpa2_apln, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_apln_DE_results.csv", row.names = TRUE)
write.csv(dpa3_apln, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_apln_DE_results.csv", row.names = TRUE)
write.csv(dpa7_apln, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_apln_DE_results.csv", row.names = TRUE)

# ------------------- pseudobulk for end ----------------------- #
# get the counts matrix
counts_end <- cts.split.modified$End

# generate sample level metadata
colData <- data.frame(samples = colnames(counts_end))

# create the metadata
colData <- colData %>%
  mutate(condition = ifelse(grepl("untreated", samples), "untreated",
                            ifelse(grepl("0dpa", samples), "0dpa",
                                   ifelse(grepl("1dpa", samples), "1dpa",
                                          ifelse(grepl("2dpa", samples), "2dpa",
                                                 ifelse(grepl("3dpa", samples), "3dpa",
                                                        ifelse(grepl("mock", samples), "mock",
                                                               ifelse(grepl("7dpa", samples), "7dpa", NA)))))))) %>% column_to_rownames(var = 'samples')

# set the condition to a factor vairable
groups <- as.factor(colData$condition)

# create the DGElist object
y <- DGEList(counts=counts_end,group=groups)
# filter by expression
keep <- filterByExpr(y)
# filter by library size
y <- y[keep,,keep.lib.sizes=FALSE]
# normalize by library size
y <- normLibSizes(y)
# calculate the normalization factors
y <- calcNormFactors(y)
normalized_end <- cpm(y)
# create the model design
design <- model.matrix(~0+groups)
# estimate dispersion
y <- estimateDisp(y, design)
colnames(design) <- levels(groups)
# fit the model
fit <- glmQLFit(y, design, robust=TRUE)

dpa0_end_qlf <- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0))
dpa0_end <- as.data.frame(dpa0_end_qlf)

dpa1_end_qlf <- glmQLFTest(fit, contrast = c(0,1,0,0,0,-1,0))
dpa1_end <- as.data.frame(dpa1_end_qlf)

dpa2_end_qlf <- glmQLFTest(fit, contrast = c(0,0,1,0,0,-1,0))
dpa2_end <- as.data.frame(dpa2_end_qlf)

dpa3_end_qlf <- glmQLFTest(fit, contrast = c(0,0,0,1,0,-1,0))
dpa3_end <- as.data.frame(dpa3_end_qlf)

dpa7_end_qlf <- glmQLFTest(fit, contrast = c(0,0,0,0,1,-1,0))
dpa7_end <- as.data.frame(dpa7_end_qlf)

# save each DE results data frame to CSV
write.csv(dpa0_end, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa0_end_DE_results.csv", row.names = TRUE)
write.csv(dpa1_end, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa1_end_DE_results.csv", row.names = TRUE)
write.csv(dpa2_end, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa2_end_DE_results.csv", row.names = TRUE)
write.csv(dpa3_end, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa3_end_DE_results.csv", row.names = TRUE)
write.csv(dpa7_end, "/Users/sm2949/Desktop/SingleCellV2WithinClusterDE/dpa7_end_DE_results.csv", row.names = TRUE)
