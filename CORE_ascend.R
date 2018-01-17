library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

devtools::load_all("~/CodeRepositories/ascend-dev")
object <- readRDS("/Volumes/Anne's External HD/mmMeninges_scRNA/mmMeninges_scRNA_Clustered_V2.rds")
conservative = TRUE
ngenes = NULL
windows = seq(0.025:1, by=0.025)
remove_outlier = TRUE

# User inputs a EMSet
if (class(object) == "EMSet"){
  # Making sure user has run PCA and reduced dimensions
  if ( !object@Log$PCA ){
    stop("Please run RunPCA followed by ReduceDimensions on this object before running this function.")
  }
  if ( is.null(object@PCA$PCA) ){
    stop("Please reduce this dataset with RunPCA.")
  }
} else {
  stop("Please supply an EMSet.")
}
  
  # Generate distance matrix and hclust objects as normal
print("Performing unsupervised clustering...")
# Reduce PCA to 20 dimensions
if (ncol(object@PCA$PCA) > 20){
  pca.matrix <- object@PCA$PCA[,1:20]
} else{
  pca.matrix <- object@PCA$PCA
}

# Generate a distance matrix from the PCA object
distance.matrix <- stats::dist(pca.matrix)
original.tree <- stats::hclust(distance.matrix, method="ward.D2")
original.clusters <- unname(dynamicTreeCut::cutreeDynamic(original.tree,
                                                          distM=as.matrix(distance.matrix), 
                                                          verbose=0))
# Check if outlier clusters have been generated
if (0 %in% original.clusters){
  outlier_barcode_list <- c()
  if (remove_outlier){
    limit <- 0
    while ((0 %in% original.clusters) && (limit <= 5)){
      # If this is the next time the loop has been run, add the old barcodes to 
      # the outlier_barcode_list to store in the object 
      if (exists("outlier_barcodes")){
        outlier_barcode_list <- c(outlier_barcode_list, outlier_barcodes)
      }
      
      # Get indices of outliers
      outlier_idx <- which(original.clusters == 0)
      print(sprintf("%i outliers detected in dataset. Removing cells and repeating PCA.", length(outlier_idx)))
      outlier_barcodes <- colnames(object@ExpressionMatrix)[outlier_idx]
      keep_barcodes <- colnames(object@ExpressionMatrix)[-outlier_idx]
      object <- SubsetCells(object, keep_barcodes)
      object <- RunPCA(object)
      pca.matrix <- object@PCA$PCA[,1:20]
      distance.matrix <- stats::dist(pca.matrix)
      original.tree <- stats::hclust(distance.matrix, method="ward.D2")
      original.clusters <- unname(dynamicTreeCut::cutreeDynamic(original.tree,
                                                                distM=as.matrix(distance.matrix), verbose=0))
      if (!(0 %in% original.clusters)){
        outlier_barcode_list <- c(outlier_barcode_list, outlier_barcodes)
        break
      } else{
        limit <- limit + 1
      }
    }
  } else{
    stop("Your dataset may contain cells that are too distinct from the main 
         population of cells. We recommend you run this function with 
         'remove_outlier = TRUE' or check the cell-cell normalisation of your
         dataset.")
  }
}

# Generate clustering matrix
# Features number of clusters produced at varying dendrogram cut heights
# Perform 40 dynamic tree cuts at varying heights
print("Generating clusters by running dynamicTreeCut at different heights...")

# Need to optimise this steps
print(system.time(cluster.list <- BiocParallel::bplapply(windows, RetrieveCluster,
                        hclust.obj = original.tree,
                        distance.matrix = distance.matrix)))
print(system.time(height.list <- lapply(windows, function(x) x)))
  
cluster.matrix <- GenerateClusteringMatrix(original.clusters, cluster.list)
  
print("Calculating rand indices...")
# Generate RAND index
rand.idx.matrix <- GenerateRandIndexMatrix(cluster.matrix, original.clusters, cluster.list)
  
print("Calculating stability values...")
# Generate Stability Values
rand.idx.matrix <- GenerateStabilityValues(rand.idx.matrix)
  
print("Aggregating data...")
# Add new cluster counts to rand matrix
cluster.counts <- BiocParallel::bplapply(cluster.list, function(x) length(unique(x)) )
rand.idx.matrix$cluster_count <- as.vector(as.numeric(cluster.counts))
rand.idx.matrix$stability_count <- rand.idx.matrix$stability_count/40
  
print("Finding optimal number of clusters...")
# Find Optimal Values
key.stats <- BuildKeyStat(rand.idx.matrix)
optimal.idx <- FindOptimalResult(key.stats, conservative = conservative)
optimal.cluster.list <- cluster.list[[optimal.idx]]
optimal.tree.height <- height.list[[optimal.idx]]
optimal.cluster.number <- cluster.counts[[optimal.idx]]
  
print("Optimal number of clusters found! Returning output...")
cell.labels <- original.tree$labels
  
# Add information to the EMSet object
output.list <- list(
  DistanceMatrix = distance.matrix,
  Hclust = original.tree,
  PutativeClusters = setNames(original.clusters, cell.labels),
  ClusteringMatrix = cluster.matrix,
  Clusters = setNames(optimal.cluster.list, cell.labels),
  NumberOfClusters = optimal.cluster.number,
  OptimalTreeHeight = optimal.tree.height,
  KeyStats = key.stats,
  RandMatrix = rand.idx.matrix
)

if (exists("outlier_barcode_list")){
  output.list$UnclusteredCells <- outlier_barcode_list
  object@Log$UnclusteredCells <- outlier_barcode_list
}

# Append all of this information to the EMSet object
object@Clusters <- output.list
object@CellInformation$cluster <- unlist(optimal.cluster.list)
saveRDS(object, "ClusteredEMSet.rds")
