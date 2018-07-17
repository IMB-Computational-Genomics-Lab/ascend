################################################################################
#
# ascend_clustering.R
# description: Functions related to the clustering of data.
#
################################################################################

#' @export
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

#' @export
findOptimalResult <- function(key_stats, conservative = TRUE, nres = 40){
  # Get stability values for most clusters and least clusters
  # Which one is greater than the other one?
  # # FindOptimalIndex
  max_min_idx <- c(1,nres)
  stability <- key_stats$Stability
  endpoint_stability_idx <- max_min_idx[which(stability[max_min_idx] ==
                                                max(stability[max_min_idx]))]
  if (length(endpoint_stability_idx) > 1){
    endpoint_stability_idx <- 1
  }
  endpoint_stability_value <- stability[endpoint_stability_idx]
  
  if (endpoint_stability_value >= 0.5){
    optimal_idx <- endpoint_stability_idx
  } else{
    if (stability[1] != stability[nres]){
      left_plateau <- getConsecSeq(which(stability == stability[1]), direction = "forward")
      right_plateau <- getConsecSeq(which(stability == stability[nres]), direction = "reverse")
    } else{
      # Remove indices in right plateau that are in left plateau
      left_plateau <- getConsecSeq(which(stability == stability[1]), direction = "forward")
      right_plateau <- getConsecSeq(which(stability == stability[nres]), direction = "reverse")
      right_plateau <- right_plateau[-(which(right_plateau %in% left_plateau))]
    }
    
    # Get rand indexes
    rand_idx_values <- key_stats$RandIndex[-c(left_plateau, right_plateau)]
    unique_rand_idx <- unique(rand_idx_values)
    
    consecutive_vals <- list()
    
    # Search for stability between the plateau
    for (idx in unique_rand_idx){
      # Get the indices where the index is present
      indexes <- which(key_stats$RandIndex == idx)
      
      # Remove indices that cross over into the plateau
      indexes <- setdiff(indexes, left_plateau)
      indexes <- setdiff(indexes, right_plateau)
      
      # If the indices plateau, get the consecutive values
      if (length(indexes) > 1){
        consecutive_indexes <- getConsecSeq(indexes, direction = "forward")
        consecutive_vals[[as.character(idx)]] <- consecutive_indexes
      }
    }
    
    # List of consecutive indexes for each unique rand index
    unique_rand_length <- lapply(consecutive_vals, length)
    max_length <- max(unlist(unique_rand_length))
    max_idx <- names(unique_rand_length)[which(unique_rand_length == max_length)]
    
    if (length(max_idx) > 1){
      # Conservative Check
      if (conservative == FALSE){
        max_idx <- max_idx[1]
      } else{
        max_idx <- max_idx[length(max_idx)]
      }
      optimal_idx <- min(which(key_stats$RandIndex == max_idx))
    } else{
      # Get Minimum
      min_idx <- lapply(max_idx, function(x) return(min(consecutive_vals[[x]])))
      optimal_idx <- min(unlist(min_idx))
    }
  }
  return(optimal_idx)
}

#' @export
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

#' @export
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

#' @export
chooseNew <- function(n = NULL, k = NULL) {
  n <- c(n); out1 <- rep(0,length(n));
  for (i in c(1:length(n)) ){
    if ( n[i]<k ) {out1[i] <- 0}
    else {out1[i] <- choose(n[i],k) }
  }
  return(out1)
}

#' @export
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

#' @export
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
  rand_idx_matrix <-as_data_frame(cbind(cluster_index_consec, cluster_index_ref))
  rand_idx_matrix$order <-row.names(rand_idx_matrix)
  return(rand_idx_matrix)
}

#' @export
retrieveCluster <- function(height_list, hclust_obj = NULL, distance_matrix = NULL){
  clusters <- sapply(height_list, function(height) dynamicTreeCut::cutreeDynamic(hclust_obj,
                                                                                 distM=as.matrix(distance_matrix),
                                                                                 minSplitHeight=height, verbose=0))
  colnames(clusters) <- height_list
  return(clusters)
}

#' @include ascend_objects.R
#' @export
setGeneric("runCORE", def = function(object, 
                                     ..., 
                                     conservative,
                                     nres,
                                     remove.outlier) {
  standardGeneric("runCORE")  
})

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
#' @param nres Number of resolutions to test, ranging from 20 to 100. Default: 40.
#' @param remove.outlier Remove cells that weren't assigned a cluster with
#' dynamicTreeCut. This is indicative of outlier cells within the sample.
#' Default: FALSE.
#' @return An \code{\linkS4class{EMSet}} with cluster information loaded into the
#' Clusters slot.
#' @include ascend_objects.R
#' @include ascend_dimreduction.R
#' @importFrom stats dist hclust setNames
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom BiocParallel bplapply
#' @export
#'
setMethod("runCORE", signature("EMSet"), function(object, 
                                                  conservative = TRUE,
                                                  nres = 40, 
                                                  remove.outlier = FALSE){
  
  # Fill in missing arguments
  if (missing(conservative)){
    conservative <- TRUE
  }
  if (missing(nres)){
    nres <- 40
  }
  if (missing(remove.outlier)){
    remove.outlier <- TRUE
  }
                                 
  # Checks
  ## 1. Make sure it is PCA reduced
  if(is.null(SingleCellExperiment::reducedDim(object, "PCA"))){
    stop("Please run runPCA on this object before using this function.")
  }
  
  # Prepare arguments
  ## Use top 20 PCs or all PCs if less than 20
  pca_matrix <- SingleCellExperiment::reducedDim(object, "PCA")
  if (ncol(pca_matrix) > 20){
    pca_matrix <- pca_matrix[ , 1:20]
  }
  
  # If user has set their own resolutions - generate windows
  if (nres != 40){
    minres <- 20
    if (nres < minres | nres > 100){
      stop(sprintf("nres should be set to a value between %i and 100. Please try again.", minres))
    }
  }
  
  # Generate sliding windows
  windows <- seq((1/nres):1, by = 1/nres)
  
  # Generate a distance matrix from the PCA matrix
  print("Calculating distance matrix...")
  distance_matrix <- stats::dist(pca_matrix)
  print("Generating hclust object...")
  original_tree <- stats::hclust(distance_matrix, method="ward.D2")
  print("Using dynamicTreeCut to generate reference set of clusters...")
  original_clusters <- unname(dynamicTreeCut::cutreeDynamic(original_tree,
                                                            distM=as.matrix(distance_matrix),
                                                            verbose=0))
  
  print("Checking if outliers are present...")
  if (0 %in% original_clusters){
    outlier_barcode_list <- c()
    if (remove.outlier){
      limit <- 0
      while ((0 %in% original_clusters) && (limit <= 5)){
        # If this is the next time the loop has been run, add the old barcodes to
        # the outlier_barcode_list to store in the object
        if (exists("outlier_barcodes")){
          outlier_barcode_list <- c(outlier_barcode_list, outlier_barcodes)
        }
        
        # Get indices of outliers
        outlier_idx <- which(original_clusters == 0)
        print(sprintf("%i outliers detected in dataset. Removing cells and repeating PCA.", length(outlier_idx)))
        outlier_barcodes <- colnames(object)[outlier_idx]
        object <- object[ , !(colnames(object) %in% outlier_barcodes)]
        object <- calculateQC(object)
        object <- runPCA(object)
        pca_matrix <- SingleCellExperiment::reducedDim(object, "PCA")[,1:20]
        distance_matrix <- stats::dist(pca_matrix)
        original_tree <- stats::hclust(distance_matrix, method="ward.D2")
        original_clusters <- unname(dynamicTreeCut::cutreeDynamic(original_tree,
                                                                  distM=as.matrix(distance_matrix), verbose=0))
        if (!(0 %in% original_clusters)){
          outlier_barcode_list <- c(outlier_barcode_list, outlier_barcodes)
          break
        } else{
          limit <- limit + 1
        }
      }
    } else{
      stop("Your dataset may contain cells that are too distinct from the main
           population of cells. We recommend you run this function with
           'remove.outlier = TRUE' or check the cell-cell normalisation of your
           dataset.")
    }
    }
  
  # Generate clustering matrix
  # Features number of clusters produced at varying dendrogram cut heights
  # Perform 40 dynamic tree cuts at varying heights
  print("Generating clusters by running dynamicTreeCut at different heights...")
  
  # Run dynamicTreeCut at different resolutions
  nworkers <- BiocParallel::bpnworkers(BiocParallel::bpparam())
  
  if(nworkers > 1){
    chunked_windows <- split(windows, 1:nworkers)
  } else{
    chunked_windows <- windows
  }
  
  # Parallel Process
  cluster_list <- BiocParallel::bplapply(chunked_windows, retrieveCluster,
                                         hclust_obj = original_tree,
                                         distance_matrix = distance_matrix)
  
  cluster_matrix <- do.call("cbind", cluster_list)
  cluster_matrix <- as.data.frame(cluster_matrix)
  cluster_matrix$REF <- original_clusters
  cluster_matrix <- cluster_matrix[ , c(as.character(windows), "REF")]
  rownames(cluster_matrix) <- original_tree$labels
  
  print("Calculating rand indices...")
  rand_idx_matrix <- calcRandIndexMatrix(cluster_matrix, original_clusters)
  
  
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
  cell_labels <- original_tree$labels
  
  # Add information to the EMSet object
  output_list <- list(
    distanceMatrix = distance_matrix,
    hClust = original_tree,
    putativeClusters = stats::setNames(original_clusters, cell_labels),
    clusteringMatrix = cluster_matrix,
    clusters = stats::setNames(optimal_cluster_list, cell_labels),
    nClusters = as.numeric(optimal_cluster_number),
    optimalTreeHeight = as.numeric(optimal_tree_height),
    keyStats = key_stats
  )
  
  log <- progressLog(object)
  log$clustering <- list(clustering = TRUE,
                         nres = nres,
                         remove_outlier = remove.outlier,
                         conservative = conservative)
  if (exists("outlier_barcode_list")){
    output_list$unclustered <- outlier_barcode_list
    log$clustering$unclustered_cells <- outlier_barcode_list
    object <- object[ , !(colnames(object) %in% outlier_barcode_list)]
  }
  
  # Append all of this information to the EMSet object
  clusterAnalysis(object) <- output_list
  colInfo(object)$cluster <- optimal_cluster_list
  progressLog(object) <- log
  
  return(object)
})
