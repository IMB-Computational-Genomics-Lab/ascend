# Rebuilding AEMSet object
file.dir <- "/Users/a.senabouth/Data/IPSCRetina_scRNA_Aggr_V2/outs/filtered_gene_bc_matrices_mex/GRCh38p7"
devtools::load_all("/Users/a.senabouth/CodeRepositories/ASCEND")
aem.set <- CellRangerToASCEND("/Users/a.senabouth/Data/IPSCRetina_scRNA_Aggr_V2", "GRCh38p7")

filtered.set <- FilterByOutliers(aem.set, CellThreshold = 3, ControlThreshold = 3)
filtered.set <- FilterByCustomControl("Mt", 20, filtered.set)
filtered.set <- FilterByCustomControl("Rb", 50, filtered.set)
filtered.set <- FilterByExpressedGenesPerCell(filtered.set, 0.01)

# Run Scran Normalisation
scran.obj <- scranNormalise(filtered.set, quickCluster = FALSE)
rle.obj <- NormaliseByRLE(filtered.set)
rle.obj <- ExcludeControl(rle.obj, "Mt")
rle.obj <- ExcludeControl(rle.obj, "Rb")

# PCA
pca.obj <- RunPCA(scran.obj)
PlotPCAVariance(pca.obj, n = 100)
pca.obj <- ReduceDimensions(pca.obj, n = 20)
#
# # Cluster
# clustered.obj <- FindOptimalClusters(pca.obj)
#
# # Run Differential Expression
# ### Cluster differential expression function
# cell.information <- clustered.obj@CellInformation
# batch.info <- unlist(cell.information$batch)
# batch.info[which(batch.info == 1)] <- "Positive"
# batch.info[which(batch.info == 2)] <- "Negative"
# cell.information$THY1 <- batch.info
# clustered.obj@CellInformation <- cell.information
#
# clean.obj <- SubsetCluster(clustered.obj, c("1"))
# clean.pca <- RunPCA(clean.obj)
# clean.pca <- ReduceDimensions(clean.pca, 20)
# clean.cluster <- FindOptimalClusters(clean.pca)
#
# clean.de.expression <- RunDiffExpression(clean.cluster, column = "cluster", conditions = c("1", "2", "3"))
#
# clean.pca <- PlotPCA(clean.cluster, column = "cluster")
# clean.tsne <- PlotTSNE(clean.cluster, PCA = TRUE, column = "cluster")
# clean.mds <- PlotMDS(clean.cluster, PCA = TRUE, column = "cluster")
