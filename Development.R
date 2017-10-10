# Rebuilding AEMSet object
# file.dir <- "/Users/a.senabouth/Data/IPSCRetina_scRNA_Aggr_V2/outs/filtered_gene_bc_matrices_mex/GRCh38p7"
devtools::load_all("/Users/a.senabouth/CodeRepositories/ASCEND")
#
library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar = TRUE), default = TRUE)
#
# aem.set <- CellRangerToASCEND("/Users/a.senabouth/Data/IPSCRetina_scRNA_Aggr_V2", "GRCh38p7")
# aem.set <- ConvertGeneAnnotation(aem.set, "gene_symbol", "ensembl_id")
#
# # Load reference cell cycle genes from scran
# hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
#
# aem.set <- scranCellCycle(aem.set, hs.pairs)
# aem.set <- ConvertGeneAnnotation(aem.set, "ensembl_id", "gene_symbol")
#
# filtered.set <- FilterByOutliers(aem.set, CellThreshold = 3, ControlThreshold = 3)
# filtered.set <- FilterByCustomControl("Mt", 20, filtered.set)
# filtered.set <- FilterByCustomControl("Rb", 50, filtered.set)
# filtered.set <- FilterByExpressedGenesPerCell(filtered.set, 1)
#
# # Run Scran Normalisation
# scran.obj <- scranNormalise(filtered.set, quickCluster = FALSE)
#
# # PCA
# pca.obj <- RunPCA(scran.obj)
# pca.obj <- ReduceDimensions(pca.obj, 20)
#
# # Clusters
# clustered.obj <- FindOptimalClusters(pca.obj)
# rgc.obj <- clustered.obj
# devtools::use_data(rgc.obj)
data(rgc.obj)
mds.plot <- PlotMDS(rgc.obj, PCA = TRUE, condition = "batch")
pca.plot <- PlotPCA(rgc.obj, dim1 = 1, dim2 = 2, condition = "batch")
tsne.plot <- PlotTSNE(rgc.obj, PCA = TRUE, condition = "cluster")
# Confounding Factors
# RegressConfoundingFactors function
#candidate.genes <- c("CDK4","CCND1","NOC2L","ATAD3C", "CCNL2", "RP5-902P8.12")
#regressed.obj <- RegressConfoundingFactors(scran.obj, candidate.genes)

#expression.matrix <- GetExpressionMatrix(regressed.obj, format = "data.frame")

