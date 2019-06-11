################################################################################
#
# ascend_clustering.R
# description: Functions related to the clustering of data.
#
################################################################################

getConsecSeq <- function(x, direction = c("forward", "reverse")){
  consecutive <- c()
  # This is for when we are looking for consecutive numbers in the FORWARD direction
  if (direction == "forward"){
    for (i in 1:length(x)){
      if (i != 1){
        current_idx <- i
        previous_idx <- i - 1
        previous_val <- x[previous_idx]
        current_val <- x[current_idx]
        step_val <- current_val - previous_val
      } else{
        current_idx <- i
        next_idx <- i + 1
        current_val <- x[current_idx]
        next_val <- x[next_idx]
        step_val <- next_val - current_val
      }
      
      if (step_val == 1){
        if (length(consecutive) > 0){
          if (current_val - consecutive[length(consecutive)] == 1){
            consecutive <- c(consecutive, current_val)
          }
        } else{
          consecutive <- c(consecutive, current_val)
        }
      }
    }
  }
  
  # For when we are looking for consecutive numbers in the reverse direction
  if (direction == "reverse"){
    # Do it backwards
    start_point <- length(x)
    end_point <- 1
    
    for(i in start_point:end_point){
      if (i != start_point){
        current_idx <- i
        previous_idx <- i + 1
        previous_val <- x[previous_idx]
        current_val <- x[current_idx]
        step_val <- previous_val - current_val
      } else{
        current_idx <- i
        next_idx <- i - 1
        current_val <- x[current_idx]
        next_val <- x[next_idx]
        step_val <- current_val - next_val
      }
      
      # If this number is consecutive with its neighbour, see if it is
      # consecutive with the current list...
      if (step_val == 1){
        if (length(consecutive) > 0){
          if (consecutive[length(consecutive)] - current_val == 1){
            consecutive <- c(consecutive, current_val)
          }
        } else{
          # If it's the first consecutive value, then just add it to the list
          consecutive <- c(consecutive, current_val)
        }
      }
    }
    consecutive <- rev(consecutive)
  }
  return(consecutive)
}

findOptimalResult <- function(key_stats, conservative = TRUE, nres = 40){
  
  # Get stability values for most clusters and least clusters
  # Which one is greater than the other one?
  max_min_idx <- c(1, nres)
  stability <- key_stats$Stability
  minmax_stabilities <- stability[max_min_idx]
  
  # 1. Optimal result is lowest resolution if plateau accounts for at least 50% 
  # of results
  # 2. Otherwise, check if high resolution plateau is more stable than the
  # cumulative stability of the transition region.
  if ((minmax_stabilities[2]) >= 0.5){
    optimal_stability_idx <- nres
  } else{
    # Find where transition region begins by finding
    # where low resolution plateau ends
    # and where high resolution plateu begins
    if (minmax_stabilities[1] != minmax_stabilities[2]){
      left_plateau <- getConsecSeq(which(stability == stability[1]), direction = "forward")
      right_plateau <- getConsecSeq(which(stability == stability[nres]), direction = "reverse")
    } else{
      left_plateau <- getConsecSeq(which(stability == stability[1]), direction = "forward")
      right_plateau <- getConsecSeq(which(stability == stability[nres]), direction = "reverse")
      right_plateau <- right_plateau[which(!(right_plateau %in% left_plateau))]
    }
    
    # Get stabilities of transition region
    transition_stabilities <- stability[-(c(left_plateau,right_plateau))]
    plateaus <- c(left_plateau, right_plateau)
    
    # If left plateau is longer than transition region, it is optimal
    if (length(left_plateau) > length(transition_stabilities)){
      optimal_stability_idx <- right_plateau[1] 
    } else{
      # Otherwise get the most stable result in the transition region
      max_stability_indices <- which(stability == max(transition_stabilities))
      # Remove indices that are in either plateaus
      max_stability_indices <- max_stability_indices[which(!max_stability_indices %in% plateaus)]
      # Determine if equally stable regions have different number of clusters.
      putative_clusters <- unique(key_stats[max_stability_indices, "ClusterCount"])
      if (length(putative_clusters) > 1){
        # If conservative flag is true, use the lower cluster count
        if (conservative){
          putative_cluster <- min(putative_clusters) 
        } else{
          putative_cluster <- max(putative_clusters)          
        }
      } else{
        putative_cluster <- putative_clusters
      }
      # Get optimal index
      optimal_stability_idx <- min(which(key_stats$ClusterCount == putative_cluster & key_stats$ConsecutiveRI == 1))
    }
  }
  return(optimal_stability_idx)
}

buildKeyStat <- function(rand_matrix = NULL, nres = 40){
  step <- signif(1/nres, digits = 3)
  # Set up a key stat dataframe
  rand.matrix <- as.data.frame(rand_matrix)
  key_stats_df <- cbind(as.numeric(rand_matrix$order)*step,
                        rand_matrix$stability_count,
                        rand_matrix$cluster_index_ref,
                        rand_matrix$cluster_index_consec,
                        rand_matrix$cluster_count)
  colnames(key_stats_df) <- c('Height',
                              'Stability',
                              'RandIndex',
                              'ConsecutiveRI',
                              "ClusterCount")
  key_stats_df <- as.data.frame(key_stats_df)
  key_stats_df$Height <-as.character(key_stats_df$Height)
  return(key_stats_df)
}

calcStability <- function(rand_idx_matrix = NULL, nres = 40){
  # WITHIN GENERATE STABILITY VALUES
  stability_values <- rand_idx_matrix$cluster_index_consec
  
  # RATIONALE BEHIND COUNTER STEPS
  # Basically, we are looking for which value remains constant the longest
  # Set up counter to keep track for each column
  general_counter <- rep(0, length(stability_values))
  
  # For first result...
  general_counter[1] <- 1
  
  # Loop over the rest
  for (i in 1:length(stability_values)){
    if (i < nres){
      # Check forward
      if (stability_values[i] == stability_values[i+1]){
        general_counter[i+1] <- general_counter[i]+1
      } else{
        general_counter[i+1] <- 1
      }
    } else{
      if (stability_values[i] == stability_values[i-1]){
        general_counter[i] <- general_counter[i-1]+1
      } else{
        general_counter[i] <- 1
      }
    }
  }
  
  # Reset the counter to find where there is no change
  flat_counter <- general_counter
  for (i in 1:length(general_counter)){
    if (i != nres){
      if (flat_counter[i] == 1 && flat_counter[i+1]==1){
        flat_counter[i] <-0
      }
    } else{
      if (flat_counter[i] == 1 && flat_counter[i-1]==1){
        flat_counter[i] <-0
      }
    }
  }
  
  if (0 %in% flat_counter){
    flat_idx <- which(flat_counter == 0)
  } else{
    flat_idx <- c(1)
  }
  
  # Reset the counter to the count values
  adjusted_counter <- general_counter
  
  # Setup beginning of the counter from the left-most side of the flat point
  if (flat_idx[1] == 1){
    adjusted_counter[1] <- 1
  } else{
    adjusted_counter[1:flat_idx[1] - 1] <- general_counter[flat_idx[1] - 1]
  }
  
  # Set up the end of the counter from the right-most side of the flat point
  # Situation 1 - Where the flat point goes right to the end of the counter
  # Situation 2 - Where the flat point does not end
  if (flat_idx[length(flat_idx)] == length(general_counter)){
    adjusted_counter[flat_idx[length(flat_idx)]] <- 1
  } else {
    adjusted_counter[(flat_idx[length(flat_idx)]+1):length(general_counter)] <- general_counter[length(general_counter)]
    adjusted_counter[flat_idx[length(flat_idx)]] <- 1
  }
  
  # Setup the counter in the middle - if required
  if (length(flat_idx) > 1){
    for (i in 2:length(flat_idx)-1){
      adjusted_counter[(flat_idx[i]+1):(flat_idx[i+1]-1)] <- general_counter[flat_idx[i+1]-1]
      adjusted_counter[flat_idx[i]] <- 1
    }
  }
  
  rand_idx_matrix$stability_count <- adjusted_counter
  return(rand_idx_matrix)
}


chooseNew <- function(n = NULL, k = NULL) {
  n <- c(n); out1 <- rep(0,length(n));
  for (i in c(1:length(n)) ){
    if ( n[i]<k ) {out1[i] <- 0}
    else {out1[i] <- choose(n[i],k) }
  }
  return(out1)
}

calcRandIndex <- function (x, adjust = TRUE)
{
  # Setting up variables for calculations
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  nn <- 0
  m <- nrow(x)
  n <- ncol(x)
  
  # Loop iterating over ranges
  for (i in 1:m) {
    c <- 0
    for (j in 1:n) {
      a <- a + chooseNew(x[i, j], 2)
      nj <- sum(x[, j])
      c <- c + chooseNew(nj, 2)
    }
    ni <- sum(x[i, ])
    b <- b + chooseNew(ni, 2)
    nn <- nn + ni
  }
  if (adjust) {
    d <- chooseNew(nn, 2)
    adrand <- (a - (b * c)/d)/(0.5 * (b + c) - (b * c)/d)
    return(adrand)
  }
  else {
    b <- b - a
    c <- c - a
    d <- chooseNew(nn, 2) - a - b - c
    rand <- (a + d)/(a + b + c + d)
    return(rand)
  }
}

calcRandIndexMatrix <- function(cluster_matrix, original_clusters, cluster_list){
  # Values to output to
  cluster_index_consec <-list()
  cluster_index_ref <-list()
  
  # Grab column one to seed values for the new table
  input_table <- table(cluster_matrix[,1], original_clusters)
  cluster_index_consec[[1]] <- 1
  cluster_index_ref[[1]] <- calcRandIndex(input_table)
  
  # Populate the rest of the matrix
  for (i in 2:(ncol(cluster_matrix) - 1)){
    cluster_index_consec[[i]] <- calcRandIndex(table(cluster_matrix[, i], cluster_matrix[ ,i-1]))
    cluster_index_ref[[i]] <- calcRandIndex(table(cluster_matrix[, i], original_clusters))
  }
  
  # Organise rand indexes into a table, assign order/result numbers
  cluster_index_consec <-unlist(cluster_index_consec)
  cluster_index_ref <-unlist(cluster_index_ref)
  rand_idx_matrix <-as.data.frame(cbind(cluster_index_consec, cluster_index_ref))
  rand_idx_matrix$order <-row.names(rand_idx_matrix)
  return(rand_idx_matrix)
}

retrieveCluster <- function(x, hclust_obj = NULL, distance_matrix = NULL){
  print(sprintf("Calculating clusters from cut height %1.2f", x))
  clusters <- dynamicTreeCut::cutreeDynamic(hclust_obj,
                                            distM=distance_matrix,
                                            minSplitHeight=x, verbose=0)
  output_list <- list() 
  output_list[[as.character(x)]] <- clusters
  return(output_list)
}

#' runCORE
#'
#' This function determines the optimal number of clusters for a dataset.
#' This function first generates a distance matrix and a hclust object, and then
#' cuts the tree at different heights.
#'
#' This will return an EMSet with the following objects:
#' \describe{
#'     \item{DistanceMatrix}{A distance matrix.}
#'     \item{Hclust}{A hclust object.}
#'     \item{PutativeClusters}{Cluster identities generated by dynamicTreeCut.}
#'     \item{ClusteringMatrix}{A matrix containing a cluster identities from cutting at 40 different heights.}
#'     \item{Clusters}{Optimum cluster identities for each cell.}
#'     \item{NumberOfClusters}{Number of clusters.}
#'     \item{OptimalTreeHeight}{Optimal tree height used to generate cluster identities.}
#'     \item{KeyStats}{A dataframe containing information on each generated clustering result, that is used to determine the optimal cluster number.}
#' }
#' @param object An EMSet object that has undergone PCA reduction.
#' @param conservative Use conservative (more stable) clustering result
#' (TRUE or FALSE). Default: TRUE.
#' @param dims Number of PC components to use in distance matrix generation. Default: 20
#' @param nres Number of resolutions to test, ranging from 20 to 100. Default: 40.
#' @param remove.outliers Remove cells that weren't assigned a cluster with
#' dynamicTreeCut. This is indicative of outlier cells within the sample.
#' Default: FALSE.
#' @param ... ...
#' 
#' @examples 
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Run CORE with default parameters
#' em_set <- runCORE(em_set, conservative = TRUE,
#' dims = 20, nres = 40, remove.outliers = TRUE)
#' 
#' @return An \code{\linkS4class{EMSet}} with cluster information loaded into 
#' the clusterAnalysis and colInfo slots
#' 
#' @include ascend_objects.R
#' @include ascend_dimreduction.R
#' @importFrom stats dist setNames
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom BiocParallel bplapply
#' @export
#'
setGeneric("runCORE", def = function(object, 
                                     ..., 
                                     conservative,
                                     nres,
                                     remove.outliers) {
  standardGeneric("runCORE")  
})

#' @rdname runCORE
setMethod("runCORE", signature("EMSet"), function(object, 
                                                  conservative = TRUE,
                                                  dims = 20,
                                                  nres = 40, 
                                                  remove.outliers = FALSE){
  
  # These defaults aren't getting passed on...
  if (missing(conservative)){
    conservative = TRUE
  }
  
  if (missing(dims)){
    dims = 20
  }
  
  if (missing(nres)){
    nres = 40
  }
  
  if (missing(remove.outliers)){
    remove.outliers = FALSE
  }
  
  if (!("PCA" %in% SingleCellExperiment::reducedDimNames(object))){
    stop("Please calculate PCA values using runPCA before using this function.")
  }
  
  pca_matrix <- SingleCellExperiment::reducedDim(object, "PCA")
  
  if (ncol(pca_matrix) >= dims){
    pca_matrix <- pca_matrix[ , 1:dims]
  }
  
  if (nres != 40){
    minres <- 20
    if (nres < minres | nres > 100){
      stop(sprintf("nres should be set to a value between %i and 100. Please try again.", minres))
    }
  }
  
  # Generate sliding windows
  windows <- seq((1/nres):1, by = 1/nres)
  
  generateClusters <- function(pca_matrix = NULL){
    print("Initiating clustering...")
    print("Calculating distance matrix...")
    distance_matrix <- stats::dist(pca_matrix)
    print("Generating hclust object via fastcluster...")
    loadNamespace("fastcluster")
    original_tree <- fastcluster::hclust(distance_matrix, method="ward.D2")
    print("Using dynamicTreeCut to generate reference set of clusters...")
    clusters <- unname(dynamicTreeCut::cutreeDynamic(original_tree,
                                                     distM = as.matrix(distance_matrix),
                                                     verbose=0))
    return(list(distanceMatrix = distance_matrix, hClust = original_tree, putativeClusters = clusters))
  }
  
  cluster_results <- FALSE
  outlier_cells <- c()
  clustering_result <- generateClusters(pca_matrix = pca_matrix)
  counter <- 0
  
  while (cluster_results == FALSE){
    if (0 %in% clustering_result$putativeClusters){
      if (counter < 5){
        # Add to outlier cells 
        outlier_cells <- c(outlier_cells, clustering_result$hClust$labels[which(clustering_result$putativeClusters == 0)])
        
        # Remove outlier cells from dataset
        object <- object[, !(colnames(object) %in% outlier_cells)]
        
        # Run PCA again
        object <- runPCA(object, scaling = TRUE, ngenes = 1500)
        
        pca_matrix <- SingleCellExperiment::reducedDim(object, "PCA")
        if (ncol(pca_matrix) > dims){
          pca_matrix <- pca_matrix[ , 1:dims]
        }
        
        clustering_result <- generateClusters(pca_matrix = pca_matrix)
        counter <- counter + 1
      }
      else{
        stop("Many outliers have been detected. Please review your filtering and normalisation methods before trying to cluster again.")
      }
    } else{
      cluster_results <- TRUE
    }
  }
  
  cluster_list <- sapply(windows, retrieveCluster, 
                         hclust_obj = clustering_result$hClust,
                         distance_matrix = as.matrix(clustering_result$distanceMatrix))
  
  cluster_matrix <- as.data.frame(matrix(unlist(cluster_list), ncol = length(cluster_list), byrow = FALSE))
  colnames(cluster_matrix) <- names(cluster_list)
  rownames(cluster_matrix) <- clustering_result$hClust$labels
  cluster_matrix$REF <- clustering_result$putativeClusters
  
  print("Calculating rand indices...")
  rand_idx_matrix <- calcRandIndexMatrix(cluster_matrix, clustering_result$putativeClusters)
  
  print("Calculating stability values...")
  rand_idx_matrix <- calcStability(rand_idx_matrix, nres)
  
  print("Aggregating data...")
  # Add new cluster counts to rand matrix
  cluster_counts <- apply(cluster_matrix[1:nres], 2, function(x) length(unique(x)))
  rand_idx_matrix$cluster_count <- as.vector(as.numeric(cluster_counts))
  rand_idx_matrix$stability_count <- rand_idx_matrix$stability_count/nres
  
  print("Finding optimal number of clusters...")
  key_stats <- buildKeyStat(rand_idx_matrix, nres)
  optimal_idx <- findOptimalResult(key_stats, conservative = conservative, nres)
  optimal_cluster_list <- cluster_matrix[, optimal_idx]
  optimal_tree_height <- colnames(cluster_matrix)[[optimal_idx]]
  optimal_cluster_number <- cluster_counts[[optimal_tree_height]]
  
  print("Optimal number of clusters found! Returning output...")
  cell_labels <- clustering_result$hClust$labels
  
  # Add information to the EMSet object
  clustering_result$putativeClusters = stats::setNames(clustering_result$putativeClusters, cell_labels)
  clustering_result$clusteringMatrix = cluster_matrix
  clustering_result$nClusters = as.numeric(optimal_cluster_number)
  clustering_result$clusters = optimal_cluster_list
  clustering_result$optimalTreeHeight = as.numeric(optimal_tree_height)
  clustering_result$keyStats = key_stats
  
  log <- progressLog(object)
  log$clustering <- list(clustering = TRUE,
                         nres = nres,
                         remove.outliers = remove.outliers,
                         conservative = conservative,
                         dims = dims)
  
  
  if (length(outlier_cells) > 0){
    clustering_result$unclustered <- outlier_cells
    log$clustering$unclustered_cells <- outlier_cells
  }
  
  # Append all of this information to the EMSet object
  clusterAnalysis(object) <- clustering_result
  colInfo(object)$cluster <- optimal_cluster_list
  progressLog(object) <- log
  return(object)
})
