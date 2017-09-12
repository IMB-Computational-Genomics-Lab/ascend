# Function used to determine consecutive sequences on boundaries
# Called by FindOptimalResult
GetConsecutiveSequence <- function(x, direction = c("forward", "reverse")){
  consecutive <- c()
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

      # If this number is consecutive with its neighbour, see if it is consecutive with the current list...
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

# Finds the optimal result from trends in stability
# Called by FindOptimalClusters
FindOptimalResult <- function(key.stats){
  # Get stability values for most clusters and least clusters
  # Which one is greater than the other one?
  max.min.idx <- c(1,40)
  stability <- key.stats$Stability
  endpoint.stability.idx <- max.min.idx[which(stability[max.min.idx] == max(stability[max.min.idx]))]
  if (length(endpoint.stability.idx) > 1){
    endpoint.stability.idx <- 1
  }
  endpoint.stability.value <- stability[endpoint.stability.idx]

  if (endpoint.stability.value >= 0.5){
    optimal.idx <- endpoint.stability.idx
  } else{
    if (stability[1] != stability[40]){
      left.plateau <- GetConsecutiveSequence(which(stability == stability[1]), direction = "forward")
      right.plateau <- GetConsecutiveSequence(which(stability == stability[40]), direction = "reverse")
    } else{
      left.plateau <- GetConsecutiveSequence(which(stability == stability[1]), direction = "forward")
      right.plateau <- GetConsecutiveSequence(which(stability == stability[40], direction = "reverse")[-left.plateau])
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
    if (length(max.length) > 1){
      optimal.idx <- min(consecutive.vals[[max.idx]])
    } else{
      # Get Minimun
      min.idx <- lapply(max.idx, function(x) return(min(consecutive.vals[[x]])))
      optimal.idx <- min(unlist(min.idx))
    }
  }
  return(optimal.idx)
}

# Generates a bunch of values based on other calculates.
# Called by FindOptimalClusters
BuildKeyStat <- function(rand.matrix){
  # Set up a key stat dataframe
  rand.matrix <- as.data.frame(rand.matrix)
  key.stats.df <- cbind(as.numeric(rand.matrix$order)*0.025, rand.matrix$stability_count, rand.matrix$cluster.index.ref, rand.matrix$cluster.index.consec, rand.matrix$cluster_count)
  colnames(key.stats.df) <-c('Height', 'Stability', 'RandIndex', 'ConsecutiveRI', "ClusterCount")
  key.stats.df <- as.data.frame(key.stats.df)
  key.stats.df$Height <-as.character(key.stats.df$Height)
  return(key.stats.df)
}

# Generates a series of stability values
GenerateStabilityValues <- function(rand.idx.matrix){
  stability.values <- rand.idx.matrix$cluster.index.consec

  # RATIONALE BEHIND COUNTER STEPS
  # Basically, we are looking for which value remains constant the longest
  # Set up counter to keep track for each column
  general.counter <- rep(0, length(stability.values))

  # For first result...
  general.counter[1] <- 1

  # Loop over the rest
  for (i in 2:length(stability.values)-1){
    if (stability.values[i] == stability.values[i+1]){
      general.counter[i+1] <- general.counter[i]+1
    } else{
      general.counter[i+1] <- 1
    }
  }

  # Reset the counter to find where there is no change
  flat.counter <- general.counter
  for (i in 1:length(general.counter)){
    if (flat.counter[i] == 1 & flat.counter[i+1]==1){
      flat.counter[i] <-0
    }
  }

  if (0 %in% flat.counter){
    flat.index <- which(flat.counter == 0)
  } else{
    flat.index <- c(1)
  }

  # Reset the counter to the count values
  adjusted.counter <- general.counter

  # Setup beginning of the counter from the left-most side of the flat point
  if (flat.index[1] == 1){
    adjusted.counter[1] <- 1
  } else{
    adjusted.counter[1:flat.index[1] - 1] <- general.counter[flat.index[1] - 1]
  }

  # Set up the end of the counter from the right-most side of the flat point
  # Situation 1 - Where the flat point goes right to the end of the counter
  # Situation 2 - Where the flat point does not end
  if (flat.index[length(flat.index)] == length(general.counter)){
    adjusted.counter[flat.index[length(flat.index)]] <- 1
  } else {
    adjusted.counter[(flat.index[length(flat.index)]+1):length(general.counter)] <- general.counter[length(general.counter)]
    adjusted.counter[flat.index[length(flat.index)]] <- 1
  }

  # Setup the counter in the middle - if required
  if (length(flat.index) > 1){
    for (i in 2:length(flat.index)-1){
      adjusted.counter[(flat.index[i]+1):(flat.index[i+1]-1)] <- general.counter[flat.index[i+1]-1]
      adjusted.counter[flat.index[i]] <- 1
    }
  }

  rand.idx.matrix$stability_count <- adjusted.counter
  return(rand.idx.matrix)
}

# Iterator for rand index
ChooseNew <- function(n,k) {
  n <- c(n); out1 <- rep(0,length(n));
  for (i in c(1:length(n)) ){
    if ( n[i]<k ) {out1[i] <- 0}
    else {out1[i] <- choose(n[i],k) }
  }
  return(out1)
}

# Calculates RAND index
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

GenerateRandIndexMatrix <- function(cluster.matrix, original.clusters, cluster.list){
  # Values to output to
  cluster.index.consec <-list()
  cluster.index.ref <-list()

  # Grab column one to seed values for the new table
  input.table <- table(cluster.matrix[[1]], original.clusters)
  cluster.index.consec[[1]] <-1
  cluster.index.ref[[1]] <- CalculateRandIndex(input.table)

  # Populate the rest of the matrix
  for (i in 2:(length(cluster.list))){
    cluster.index.consec[[i]] <- CalculateRandIndex(table(cluster.list[[i]],cluster.list[[i-1]]))
    cluster.index.ref[[i]] <- CalculateRandIndex(table(cluster.list[[i]], original.clusters))
  }

  # Organise rand indexes into a table, assign order/result numbers
  cluster.index.consec <-unlist(cluster.index.consec)
  cluster.index.ref <-unlist(cluster.index.ref)
  rand.idx.matrix <-as.data.frame(cbind(cluster.index.consec, cluster.index.ref))
  rand.idx.matrix$order <-row.names(rand.idx.matrix)
  return(rand.idx.matrix)
}

GenerateClusteringMatrix <- function(original.clusters, cluster.list){
  # Transform clustering results into a data frame, with each cut as a column, each cell as a row and their cluster as the value of the matrix
  cluster.matrix <- matrix(unlist(cluster.list), ncol = 40)
  cluster.matrix <-as.data.frame(cluster.matrix)
  cluster.matrix$ref <- original.clusters
  colnames(cluster.matrix) <-c(seq(1:40),'REF')

  return(cluster.matrix)
}

RetrieveCluster <- function(height, hclust.obj = NULL, distance.matrix = NULL){
  clusters <- unname(dynamicTreeCut::cutreeDynamic(hclust.obj, distM=as.matrix(distance.matrix), minSplitHeight=height, verbose=0))
  return(clusters)
}

#' FindOptimalClusters
#'
#' This function determines the optimal number of clusters for a dataset.
#' This function first generates a distance matrix and a hclust object, and then cuts the tree at different heights.
#' This will return an AEMSet with the following objects:
#' \describe{
#'     \item{DistanceMatrix}{A distance matrix}
#'     \item{Hclust}{A hclust object}
#'     \item{PutativeClusters}{Cluster identities generated by dynamicTreeCut}
#'     \item{ClusteringMatrix}{A matrix containing a cluster identities from cutting at 40 different heights}
#'     \item{Clusters}{Optimum cluster identities for each cell}
#'     \item{NumberOfClusters}{Number of clusters}
#'     \item{OptimalTreeHeight}{Optimal tree height used to generate cluster identities}
#'     \item{RandMatrix}{Rand matrix used to determine optimal number of clusters}
#' }
#' @param object An AEMSet object that has undergone PCA reduction.
#'
FindOptimalClusters <- function(object){
  # User inputs a AEMSet
  if (class(object) == "AEMSet"){
    # Making sure user has run PCA and reduced dimensions
    if ( !object@Log$PCA ){
      stop("Please run RunPCA followed by ReduceDimensions on this object before running this function.")
    }
    if ( is.null(object@PCA$PCA) ){
      stop("Please reduce this dataset with RunPCA and ReduceDimensions.")
    }
  } else {
    stop("Please supply a AEMSet object.")
  }

  # Generate distance matrix and hclust objects as normal
  print("Performing unsupervised clustering...")
  if (ncol(object@PCA$PCA) > 20){
    pca.matrix <- object@PCA$PCA[,1:20]
  } else{
    pca.matrix <- object@PCA$PCA
  }
  distance.matrix <- stats::dist(pca.matrix)
  original.tree <- stats::hclust(distance.matrix, method="ward.D2")
  original.clusters <- unname(dynamicTreeCut::cutreeDynamic(original.tree, distM=as.matrix(distance.matrix), verbose=0))

  # Generate clustering matrix
  # Features number of clusters produced at varying dendrogram cut heights
  # Perform 40 dynamic tree cuts at varying heights
  print("Generating clusters by running dynamicTreeCut at different heights...")
  cluster.list <- lapply(seq(0.025:1, by=0.025), RetrieveCluster, hclust.obj = original.tree, distance.matrix = distance.matrix)
  height.list <- lapply(seq(0.025:1, by=0.025), function(x) x)

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
  optimal.idx <- FindOptimalResult(key.stats)
  optimal.cluster.list <- cluster.list[[optimal.idx]]
  optimal.tree.height <- height.list[[optimal.idx]]
  optimal.cluster.number <- cluster.counts[[optimal.idx]]

  print("Optimal number of clusters found! Returning output...")
  cell.labels <- original.tree$labels

  # Add information to the AEMSet object
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

  # Append all of this information to the AEMSet object
  object@Clusters <- output.list
  return(object)
}
