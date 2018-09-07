## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ascend)
data(datasets)
RawSet <- datasets$RawSet
EMSet <- datasets$CompleteSet

## ----cell_info_knitr-----------------------------------------------------
knitr::kable(data.table::data.table(
  cell_barcode = c("Cell1-1", "Cell2-1", "...", "...", "Cellx-x"), 
  batch = c("1", "1", "...", "...", "x")))

## ----gene_info_knitr-----------------------------------------------------
knitr::kable(data.table::data.table(
  gene_identifier1 = c("GENE1", "GENE2", "...", "...", "GENEX"), 
  gene_identifier2 = c("IDENTIFIER1", "IDENTIFIER2", "...", "...", 
                                      "IDENTIFIERX")))

## ----set_controls, eval = FALSE------------------------------------------
#  # Example controls
#  controls <- list(Mt = grep("^Mt-", rownames(expression_matrix), ignore.case = TRUE)
#                   Rb = grep("^Rps|^Rpl", rownames(expression_matrix), ignore.case = TRUE)

## ----load_data, eval = FALSE---------------------------------------------
#  # Creating an EMSet from scratch
#  em_set <- newEMSet(assays = list(counts = expression_matrix),
#                     colInfo = col_info_data_frame,
#                     rowInfo = row_info_data_frame,
#                     controls = control_list)
#  
#  # Loading an EMSet from a SingleCellExperiment
#  em_set <- newEMSet(SingleCellExperiment,
#                     colInfo = col_info_data_frame,
#                     rowInfo = row_info_data_frame,
#                     controls = control_list)
#  

## ----get_set_examples----------------------------------------------------
# Get assays
count_matrix <- counts(EMSet)
norm_matrix <- normcounts(EMSet)
logcounts_matrix <- logcounts(EMSet)

# Set assays
counts(EMSet) <- count_matrix
normcounts(EMSet) <- norm_matrix
logcounts(EMSet) <- logcounts_matrix

# Get gene and cell information
col_info <- colInfo(EMSet)
col_data <- colData(EMSet)
row_info <- rowInfo(EMSet)
row_data <- rowData(EMSet)

# Set gene and cell information
colInfo(EMSet) <- col_info
colData(EMSet) <- col_data
rowInfo(EMSet) <- row_info
rowData(EMSet) <- row_data

# Get reduced dimensionality data
tsne_matrix <- reducedDim(EMSet, "TSNE")
pca_matrix <- reducedDim(EMSet, "PCA")

# Get cluster analysis
clusterAnalysis <- clusterAnalysis(EMSet)

# Get progress log
progressLog <- progressLog(EMSet)

## ----dataframe_accesors--------------------------------------------------
# Reduce EMSet to first ten cells and genes
tiny_EMSet <- EMSet[1:10,1:10]

# Review content in smaller dataset 
print(counts(tiny_EMSet))
print(colInfo(tiny_EMSet))
print(rowInfo(tiny_EMSet))

## ----subsetCondition-----------------------------------------------------
# Subset batch 1 from the combined EMSet
Batch1_EMSet <- subsetCondition(EMSet, by = "batch", conditions = list(batch = c(1)))

## ----update_object, eval = FALSE-----------------------------------------
#  # Read in old EMSet stored as an RDS file
#  legacy_EMSet <- readRDS("LegacyEMSet.rds")
#  
#  # Update to new object
#  # Please ensure your new object has the same name as the old object
#  legacy_EMSet <- updateObject(legacy_EMSet)
#  

## ----runTSNE-------------------------------------------------------------
EMSet <- runTSNE(EMSet, seed = 1)

## ----plotTSNE, fig.width = 6, fig.height = 5-----------------------------
tsne_plot <- plotTSNE(EMSet, group = "cluster")
tsne_plot

## ----colour_tsne, fig.width = 6, fig.height = 5--------------------------
library(ggplot2)
tsne_plot <- tsne_plot + scale_color_manual(values=c("#bb5f4c", 
                                                     "#8e5db0", 
                                                     "#729b57"))
tsne_plot

## ----convert2SCE---------------------------------------------------------
# Convert to SingleCellExperiment
SingleCellExperiment <- EMSet2SCE(EMSet)

# Revert to SingleCellExperiment
EMSet <- SCE2EMSet(SingleCellExperiment)


## ----scran_normalisation-------------------------------------------------
scran_normalised <- scranNormalise(RawSet, quickCluster = FALSE, min.mean = 1e-05)

## ----scran_cellcycle, eval = FALSE---------------------------------------
#  # Convert identifiers to ENSEMBL gene identifiers
#  ensembl_set <- convertGeneID(RawSet, new.annotation = "ensembl_gene_id")
#  
#  # Load training data from scran
#  training_data <- readRDS(system.file("exdata", "human_cycle_markers.rds",
#                                       package = "scran"))
#  
#  # Run cell cycle
#  scran_cellcycle <- scranCellCycle(ensembl_set, training_set = training_data)
#  
#  # Show cell cycle results
#  colInfo(scran_cellcycle)

## ----DESeq, eval = FALSE-------------------------------------------------
#  cluster2_vs_others <- runDESeq(EMSet, group = "cluster",
#                                 condition.a = 2, condition.b = c(1, 3),
#                                 ngenes = 1500)

## ----DESeq2, eval = FALSE------------------------------------------------
#  cluster1_vs_others <- runDESeq2(EMSet, group = "cluster",
#                                 condition.a = 1, condition.b = c(2, 3),
#                                 ngenes = 1500)

