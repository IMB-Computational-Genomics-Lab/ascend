# Find Optimal Result Function
#' GetConsecutiveSequence
#'
#' Function used to determine consecutive sequences on boundaries.
#' This is called by FindOptimalResult to determine which clustering results are
#' consecutive.
#' 
#' @param x Set of values to find the consecutive sequences in.
#' @param direction Direction to find the consecutive sequences in.
#' @return A list of consecutive numbers.
#' 
GetConsecutiveSequence <- function(x, direction = c("forward", "reverse")){
  consecutive <- c()
  # This is for when we are looking for consecutive numbers in the FORWARD direction
  if (direction == "forward"){
    for (i in 1:length(x)){
      if (i != 1){
        current.idx <- i
        previous.idx <- i - 1
        previous.val <- x[previous.idx]
        current.val <- x[current.idx]
        step.val <- current.val - previous.val
      } else{
        current.idx <- i
        next.idx <- i + 1
        current.val <- x[current.idx]
        next.val <- x[next.idx]
        step.val <- next.val - current.val
      }
      
      if (step.val == 1){
        if (length(consecutive) > 0){
          if (current.val - consecutive[length(consecutive)] == 1){
            consecutive <- c(consecutive, current.val)
          }
        } else{
          consecutive <- c(consecutive, current.val)
        }
      }
    }
  }
  
  # For when we are looking for consecutive numbers in the reverse direction
  if (direction == "reverse"){
    # Do it backwards
    start.point <- length(x)
    end.point <- 1
    
    for(i in start.point:end.point){
      if (i != start.point){
        current.idx <- i
        previous.idx <- i + 1
        previous.val <- x[previous.idx]
        current.val <- x[current.idx]
        step.val <- previous.val - current.val
      } else{
        current.idx <- i
        next.idx <- i - 1
        current.val <- x[current.idx]
        next.val <- x[next.idx]
        step.val <- current.val - next.val
      }
      
      # If this number is consecutive with its neighbour, see if it is
      # consecutive with the current list...
      if (step.val == 1){
        if (length(consecutive) > 0){
          if (consecutive[length(consecutive)] - current.val == 1){
            consecutive <- c(consecutive, current.val)
          }
        } else{
          # If it's the first consecutive value, then just add it to the list
          consecutive <- c(consecutive, current.val)
        }
      }
    }
    consecutive <- rev(consecutive)
  }
  return(consecutive)
}

#' FindOptimalResult
#'
#' Finds the optimal resolution from trends in stability. Called by 
#' \code{\link{RunCORE}}.
#' 
#' @param key.stats A data frame containing clustering results across 40 res.
#' @param conservative In the event of a scenario where there are more than one 
#' equally stable results, choose to use the lower resolution result that 
#' yields less clusters. If set to FALSE, the higher resolution that generates a 
#' greater number of clusters will be used.
#' @param nres Number of resolutions to test.
#' @return An integer.
#' 
FindOptimalResult <- function(key.stats, conservative = TRUE, nres = 40){
  # Get stability values for most clusters and least clusters
  # Which one is greater than the other one?
  # # FindOptimalIndex
  max.min.idx <- c(1,nres)
  stability <- key.stats$Stability
  endpoint.stability.idx <- max.min.idx[which(stability[max.min.idx] ==
                                                max(stability[max.min.idx]))]
  if (length(endpoint.stability.idx) > 1){
    endpoint.stability.idx <- 1
  }
  endpoint.stability.value <- stability[endpoint.stability.idx]
  
  if (endpoint.stability.value >= 0.5){
    optimal.idx <- endpoint.stability.idx
  } else{
    if (stability[1] != stability[nres]){
      left.plateau <- GetConsecutiveSequence(which(stability == stability[1]), direction = "forward")
      right.plateau <- GetConsecutiveSequence(which(stability == stability[nres]), direction = "reverse")
    } else{
      # Remove indices in right plateau that are in left plateau
      left.plateau <- GetConsecutiveSequence(which(stability == stability[1]), direction = "forward")
      right.plateau <- GetConsecutiveSequence(which(stability == stability[nres]), direction = "reverse")
      right.plateau <- right.plateau[-(which(right.plateau %in% left.plateau))]
    }
    
    # Get rand indexes
    rand.idx.values <- key.stats$RandIndex[-c(left.plateau, right.plateau)]
    unique.rand.idx <- unique(rand.idx.values)
    
    consecutive.vals <- list()
    
    # Search for stability between the plateau
    for (idx in unique.rand.idx){
      # Get the indices where the index is present
      indexes <- which(key.stats$RandIndex == idx)
      
      # Remove indices that cross over into the plateau
      indexes <- setdiff(indexes, left.plateau)
      indexes <- setdiff(indexes, right.plateau)
      
      # If the indices plateau, get the consecutive values
      if (length(indexes) > 1){
        consecutive.indexes <- GetConsecutiveSequence(indexes, direction = "forward")
        consecutive.vals[[as.character(idx)]] <- consecutive.indexes
      }
    }
    
    # List of consecutive indexes for each unique rand index
    unique.rand.length <- lapply(consecutive.vals, length)
    max.length <- max(unlist(unique.rand.length))
    max.idx <- names(unique.rand.length)[which(unique.rand.length == max.length)]
    
    if (length(max.idx) > 1){
      # Conservative Check
      if (conservative == FALSE){
        max.idx <- max.idx[1]
      } else{
        max.idx <- max.idx[length(max.idx)]
      }
      optimal.idx <- min(which(key.stats$RandIndex == max.idx))
    } else{
      # Get Minimum
      min.idx <- lapply(max.idx, function(x) return(min(consecutive.vals[[x]])))
      optimal.idx <- min(unlist(min.idx))
    }
  }
  return(optimal.idx)
}

#' BuildKeyStat
#' Generates a bunch of values based on other calculates. Called by 
#' \code{\link{RunCORE}}.
#' @param rand.matrix Rand matrix generated by RunCORE.
#' @param nres Number of resolutions inputted by user.
#' @return A labelled dataframe.
#' 
BuildKeyStat <- function(rand.matrix = NULL, nres = 40){
  step <- signif(1/nres, digits = 3)
  # Set up a key stat dataframe
  rand.matrix <- as.data.frame(rand.matrix)
  key.stats.df <- cbind(as.numeric(rand.matrix$order)*step,
                        rand.matrix$stability_count,
                        rand.matrix$cluster.index.ref,
                        rand.matrix$cluster.index.consec,
                        rand.matrix$cluster_count)
  colnames(key.stats.df) <- c('Height',
                              'Stability',
                              'RandIndex',
                              'ConsecutiveRI',
                              "ClusterCount")
  key.stats.df <- as.data.frame(key.stats.df)
  key.stats.df$Height <-as.character(key.stats.df$Height)
  return(key.stats.df)
}

#' GenerateStabilityValues
#'
#' Generates a series of stability values based on rand matrix result
#' @param rand.idx.matrix Matrix generated by RunCORE
#' @param nres Number of resolutions
#' @return A vector of values that indicate the stability of a set of rand 
#' indices
#'
GenerateStabilityValues <- function(rand.idx.matrix = NULL, nres = 40){
  # WITHIN GENERATE STABILITY VALUES
  stability.values <- rand.idx.matrix$cluster.index.consec
  
  # RATIONALE BEHIND COUNTER STEPS
  # Basically, we are looking for which value remains constant the longest
  # Set up counter to keep track for each column
  general.counter <- rep(0, length(stability.values))
  
  # For first result...
  general.counter[1] <- 1
  
  # Loop over the rest
  for (i in 1:length(stability.values)){
    if (i < nres){
      # Check forward
      if (stability.values[i] == stability.values[i+1]){
        general.counter[i+1] <- general.counter[i]+1
      } else{
        general.counter[i+1] <- 1
      }
    } else{
      if (stability.values[i] == stability.values[i-1]){
        general.counter[i] <- general.counter[i-1]+1
      } else{
        general.counter[i] <- 1
      }
    }
  }
  
  # Reset the counter to find where there is no change
  flat.counter <- general.counter
  for (i in 1:length(general.counter)){
    if (i != nres){
      if (flat.counter[i] == 1 && flat.counter[i+1]==1){
        flat.counter[i] <-0
      }
    } else{
      if (flat.counter[i] == 1 && flat.counter[i-1]==1){
        flat.counter[i] <-0
      }
    }
  }
  
  if (0 %in% flat.counter){
    flat.idx <- which(flat.counter == 0)
  } else{
    flat.idx <- c(1)
  }
  
  # Reset the counter to the count values
  adjusted.counter <- general.counter
  
  # Setup beginning of the counter from the left-most side of the flat point
  if (flat.idx[1] == 1){
    adjusted.counter[1] <- 1
  } else{
    adjusted.counter[1:flat.idx[1] - 1] <- general.counter[flat.idx[1] - 1]
  }
  
  # Set up the end of the counter from the right-most side of the flat point
  # Situation 1 - Where the flat point goes right to the end of the counter
  # Situation 2 - Where the flat point does not end
  if (flat.idx[length(flat.idx)] == length(general.counter)){
    adjusted.counter[flat.idx[length(flat.idx)]] <- 1
  } else {
    adjusted.counter[(flat.idx[length(flat.idx)]+1):length(general.counter)] <- general.counter[length(general.counter)]
    adjusted.counter[flat.idx[length(flat.idx)]] <- 1
  }
  
  # Setup the counter in the middle - if required
  if (length(flat.idx) > 1){
    for (i in 2:length(flat.idx)-1){
      adjusted.counter[(flat.idx[i]+1):(flat.idx[i+1]-1)] <- general.counter[flat.idx[i+1]-1]
      adjusted.counter[flat.idx[i]] <- 1
    }
  }
  
  rand.idx.matrix$stability_count <- adjusted.counter
  return(rand.idx.matrix)
}

#' ChooseNew
#'
#' Function to calculate rand index. Adapted from function by Steve Horvarth and
#' Luohua Jiang (UCLA, 2003).
#' 
#' @param n First clustering result.
#' @param k Second clustering result.
#' @return An integer representing a consecutive or non-consecutive clustering 
#' result
#' 
ChooseNew <- function(n = NULL, k = NULL) {
  n <- c(n); out1 <- rep(0,length(n));
  for (i in c(1:length(n)) ){
    if ( n[i]<k ) {out1[i] <- 0}
    else {out1[i] <- choose(n[i],k) }
  }
  return(out1)
}

#' CalculateRandIndex
#'
#' Calculates rand matrix based on clustering results
#' 
#' @param x Set of clustering results.
#' @param adjust TRUE or FALSE.
#' @return A matrix of rand indices based on \code{x}
#' 
CalculateRandIndex <- function (x, adjust = TRUE)
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
      a <- a + ChooseNew(x[i, j], 2)
      nj <- sum(x[, j])
      c <- c + ChooseNew(nj, 2)
    }
    ni <- sum(x[i, ])
    b <- b + ChooseNew(ni, 2)
    nn <- nn + ni
  }
  if (adjust) {
    d <- ChooseNew(nn, 2)
    adrand <- (a - (b * c)/d)/(0.5 * (b + c) - (b * c)/d)
    return(adrand)
  }
  else {
    b <- b - a
    c <- c - a
    d <- ChooseNew(nn, 2) - a - b - c
    rand <- (a + d)/(a + b + c + d)
    return(rand)
  }
}

#' GenerateRandIndexMatrix
#'
#' Generates columns in KeyStats data frame related to rand index calculations.
#' 
#' @param cluster.matrix A matrix containing clustering results generated by 
#' cutting the hclust object at 40 resolutions.
#' @param original.clusters Clusters generated by dynamicTreeCut.
#' @param cluster.list A list of clustering results to unpack.
#' @return The Key Stats dataframe containing information collated from 
#' \code{cluster.matrix}, \code{original.clusters} and \code{cluster.list}
#' 
GenerateRandIndexMatrix <- function(cluster.matrix, original.clusters, cluster.list){
  # Values to output to
  cluster.index.consec <-list()
  cluster.index.ref <-list()
  
  # Grab column one to seed values for the new table
  input.table <- table(cluster.matrix[,1], original.clusters)
  cluster.index.consec[[1]] <- 1
  cluster.index.ref[[1]] <- CalculateRandIndex(input.table)
  
  # Populate the rest of the matrix
  for (i in 2:(ncol(cluster.matrix) - 1)){
    cluster.index.consec[[i]] <- CalculateRandIndex(table(cluster.matrix[, i], cluster.matrix[ ,i-1]))
    cluster.index.ref[[i]] <- CalculateRandIndex(table(cluster.matrix[, i], original.clusters))  
  }
  
  # Organise rand indexes into a table, assign order/result numbers
  cluster.index.consec <-unlist(cluster.index.consec)
  cluster.index.ref <-unlist(cluster.index.ref)
  rand.idx.matrix <-as.data.frame(cbind(cluster.index.consec, cluster.index.ref))
  rand.idx.matrix$order <-row.names(rand.idx.matrix)
  return(rand.idx.matrix)
}

# Retrieve Cluster Function
#' RetrieveCluster
#'
#' Generates clusters using dynamicTreeCut's unsupervised cut setting.
#' @param height.list List of resolutions to cut at.
#' @param hclust.obj Hclust object generated from the distance matrix.
#' @param distance.matrix Distance matrix generated by RunCORE from PCA matrix.
#' @importFrom dynamicTreeCut cutreeDynamic
#' @return List of clusters.
#'
RetrieveCluster <- function(height.list, hclust.obj = NULL, distance.matrix = NULL){
  clusters <- sapply(height.list, function(height) dynamicTreeCut::cutreeDynamic(hclust.obj,
                                                                                 distM=as.matrix(distance.matrix), 
                                                                                 minSplitHeight=height, verbose=0))
  colnames(clusters) <- height.list
  return(clusters)
}

#' RunCORE
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
#' @param remove_outlier Remove cells that weren't assigned a cluster with 
#' dynamicTreeCut. This is indicative of outlier cells within the sample.
#' Default: FALSE.
#' @return An \code{\linkS4class{EMSet}} with cluster information loaded into the 
#' Clusters slot.
#' @examples
#' \dontrun{
#' clustered_set <- RunCORE(em.set, conservative = TRUE, 
#' windows = 40, remove_outlier = TRUE)
#' }
#' @importFrom stats dist hclust setNames
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom BiocParallel bplapply
#' @export
#'
RunCORE <- function(object, conservative = TRUE, nres = 40, remove_outlier = FALSE){
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
  
  # If user has set their own nresolutions
  if(nres != 40){
    minres <- 20
    if (nres < minres | nres > 100){
      stop(sprintf("nres should be set to a value between %i and 100. Please try again.", minres))
    }
  }
  
  # Generate sliding windows
  windows <- seq((1/nres):1, by = 1/nres)
  
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
  
  # Run dynamicTreeCut at different resolutions
  nworkers <- BiocParallel::bpnworkers(BiocParallel::bpparam())
  
  if(nworkers > 1){
    chunked.windows <- split(windows, 1:nworkers)
  } else{
    chunked.windows <- windows
  }
  
  # Parallel Process
  cluster.list <- BiocParallel::bplapply(chunked.windows, RetrieveCluster,
                                         hclust.obj = original.tree,
                                         distance.matrix = distance.matrix)
  
  # Combine Results, reorder from lowest to highest cut and add barcodes to rows
  cluster.matrix <- do.call("cbind", cluster.list)
  cluster.matrix <- as.data.frame(cluster.matrix)
  cluster.matrix$REF <- original.clusters
  cluster.matrix <- cluster.matrix[ , c(as.character(windows), "REF")]
  rownames(cluster.matrix) <- original.tree$labels
  
  # Generate rand index and stability values
  print("Calculating rand indices...")
  rand.idx.matrix <- GenerateRandIndexMatrix(cluster.matrix, original.clusters)
  
  
  print("Calculating stability values...")
  rand.idx.matrix <- GenerateStabilityValues(rand.idx.matrix, nres)
  
  print("Aggregating data...")
  # Add new cluster counts to rand matrix
  cluster.counts <- apply(cluster.matrix[1:nres], 2, function(x) length(unique(x)))
  rand.idx.matrix$cluster_count <- as.vector(as.numeric(cluster.counts))
  rand.idx.matrix$stability_count <- rand.idx.matrix$stability_count/nres
  
  print("Finding optimal number of clusters...")
  key.stats <- BuildKeyStat(rand.idx.matrix, nres)
  optimal.idx <- FindOptimalResult(key.stats, conservative = conservative, nres)
  optimal.cluster.list <- cluster.matrix[, optimal.idx]
  optimal.tree.height <- colnames(cluster.matrix)[[optimal.idx]]
  optimal.cluster.number <- cluster.counts[[optimal.tree.height]]
  
  print("Optimal number of clusters found! Returning output...")
  cell.labels <- original.tree$labels
  
  # Add information to the EMSet object
  output.list <- list(
    DistanceMatrix = distance.matrix,
    Hclust = original.tree,
    PutativeClusters = stats::setNames(original.clusters, cell.labels),
    ClusteringMatrix = cluster.matrix,
    Clusters = stats::setNames(optimal.cluster.list, cell.labels),
    NumberOfClusters = as.numeric(optimal.cluster.number),
    OptimalTreeHeight = as.numeric(optimal.tree.height),
    KeyStats = key.stats
  )
  
  if (exists("outlier_barcode_list")){
    output.list$UnclusteredCells <- outlier_barcode_list
    object@Log$UnclusteredCells <- outlier_barcode_list
  }
  
  # Append all of this information to the EMSet object
  object@Clusters <- output.list
  object@CellInformation$cluster <- unlist(optimal.cluster.list)
  
  # Update Log
  log <- object@Log
  log$Clustering <- list(Clustering = TRUE, 
                         nres = nres, 
                         remove_outlier = remove_outlier, 
                         conservative = conservative)
  object@Log <- log
  return(object)
}
