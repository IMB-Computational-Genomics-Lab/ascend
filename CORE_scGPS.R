devtools::load_all("~/CodeRepositories/ascend-dev")

#@param remove_outlier a vector containing IDs for clusters to be removed the default vector contains 0, as 0 is the cluster with singletons

Cluster_with_OutlierRemovalOption <- function(object = NULL, ngenes= 1500, windows = seq(0.025:1, by=0.025),
                                              remove_outlier = c(0)){
  
  print("Calculating distance matrix")
  #function for the highest resolution clustering (i.e. no window applied)
  firstRoundClustering <- function(object = NULL){
    exprs_mat <- assay(object)
    exprs_mat_topVar <- topvar_scGPS(exprs_mat, ngenes = ngenes)
    #Take the top variable genes
    #make a transpose
    exprs_mat_t <-t(exprs_mat_topVar)
    dist_mat <- rcpp_parallel_distance(exprs_mat_t)
    print("Performing hierarchical clustering")
    original.tree <- fastcluster::hclust(as.dist(dist_mat), method="ward.D2")
    #The original clusters to be used as the reference
    print("Finding clustering information")
    original.clusters <- unname(cutreeDynamic(original.tree, distM=as.matrix(dist_mat), verbose=0))
    original.tree$labels <- original.clusters
    return(list("tree" = original.tree, "cluster_ref" = original.clusters, "dist_mat" = dist_mat))
  }
  
  removeOutlierCluster <-function(object = NULL,object_rmOutlier = object_rmOutlier){
    #check for singletons
    firstRound_out <- firstRoundClustering(object)
    firstRound_cluster <- as.data.frame(table(firstRound_out$cluster_ref))
    cluster_toRemove <- which(firstRound_out$cluster_ref %in% object_rmOutlier)
    
    if(length(cluster_toRemove) > 0){
      print(paste0("Removing ", length(cluster_toRemove)," cells as outliers..."))
      if(length(cluster_toRemove) == ncol(object)){stop("All cells were removed. Check your outliers.")}
      object_rmOutlier <- object[,-cluster_toRemove]
      firstRound_out <- firstRoundClustering(object_rmOutlier)
    }
    return(list("firstRound_out" = firstRound_out,
                "cellsRemoved" = colnames(object[,cluster_toRemove]),
                "cellsForClustering" = colnames(object[,-cluster_toRemove]))
    )
  }
  
  firstRoundPostRemoval <- removeOutlierCluster(object = object, object_rmOutlier = remove_outlier)
  firstRound_out <-firstRoundPostRemoval$firstRound_out
  #return variables for the next step
  original.tree <- firstRound_out$tree
  original.clusters <-firstRound_out$cluster_ref
  dist_mat <-firstRound_out$dist_mat
  
  clustering_param <-list()
  for (i in 1:length(windows)){
    
    namelist =paste0("window",windows[i])
    toadd <-as.vector(cutreeDynamic(original.tree, distM=as.matrix(dist_mat),
                                    minSplitHeight=windows[i], verbose=0))
    
    print(paste0("writing clustering result for run ", i))
    clustering_param[[i]] <-  list(toadd)
    names(clustering_param[[i]]) <- namelist
  }
  
  names(clustering_param[[i]]) <- "cluster_ref"
  print("Done clustering, moving to stability calculation...")
  return(list("list_clusters" = clustering_param, "tree" = original.tree,
              "cluster_ref" = original.clusters,
              "cellsRemoved"= firstRoundPostRemoval$cellsRemoved,
              "cellsForClustering"= firstRoundPostRemoval$cellsForClustering))
  
}