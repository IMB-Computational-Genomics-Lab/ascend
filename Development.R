# Rebuilding AEMSet object
file.dir <- "/Users/a.senabouth/Data/IPSCRetina_scRNA_Aggr_V2/outs/filtered_gene_bc_matrices_mex/GRCh38p7"
devtools::load_all("/Users/a.senabouth/CodeRepositories/ASCEND")

library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

aem.set <- CellRangerToASCEND("/Users/a.senabouth/Data/IPSCRetina_scRNA_Aggr_V2", "GRCh38p7")

normalise.batches <- NormaliseBatches(aem.set)
filtered.set <- FilterByOutliers(aem.set, CellThreshold = 3, ControlThreshold = 3)
filtered.set <- FilterByCustomControl("Mt", 20, filtered.set)
filtered.set <- FilterByCustomControl("Rb", 50, filtered.set)
filtered.set <- FilterByExpressedGenesPerCell(filtered.set, 1)

# Run Scran Normalisation
scran.obj <- scranNormalise(filtered.set, quickCluster = FALSE)

# PCA
pca.obj <- RunPCA(scran.obj)
#PlotPCAVariance(pca.obj, n = 100)
pca.obj <- ReduceDimensions(pca.obj, n = 20)

# Cluster
clustered.obj <- FindOptimalClusters(pca.obj)
clean.obj <- SubsetCluster(clustered.obj, "1")

# Run PCA
pca.obj <- RunPCA(clean.obj)
pca.obj <- ReduceDimensions(pca.obj)
clean.clusters <- FindOptimalClusters(pca.obj)

# Troubleshooting Plots
cell.df <- GetCellInfo(clean.clusters)

# Label THY1 status
THY1 <- cell.df$batch
THY1[which(THY1 == "1")] <- TRUE
THY1[which(THY1 == "2")] <- FALSE
cell.df$THY1 <- THY1

# Label cells based on BRN3B/POU4F2 expression - gene involved in differentiation of RGCs
expression.matrix <- GetExpressionMatrix(clean.clusters, format = "data.frame")
POU4F2 <- rep(FALSE, ncol(expression.matrix))
POU4F2[which(expression.matrix["POU4F2", ] > 0)] <- TRUE
cell.df$POU4F2 <- POU4F2
cluster_colours <- cell.df$cluster
cluster_colours[which(cluster_colours == 1)] <- "#6c45cb"
cluster_colours[which(cluster_colours == 2)] <- "#d66d00"
cluster_colours[which(cluster_colours == 3)] <- "#007047"
cell.df$cluster_colours <- cluster_colours

object <- clean.clusters
object <- ReplaceCellInfo(object, cell.df)

