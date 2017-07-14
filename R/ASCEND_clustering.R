# Function used to determine consecutive sequences on boundaries
# Called by FindOptimalResult
GetConsecutiveSequence <- function(x){
  consecutive <- c()
  for (i in 1:length(x)){
    if (i != length(x)){
      current.idx <- i
      next.idx <- i + 1
      current.val <- x[current.idx]
      next.val <- x[next.idx]
      step.val <- next.val - current.val
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
  return(consecutive)
}

# Finds the optimal result from trends in stability
# Called by FindOptimalClusters
FindOptimalResult <- function(key.stats.df){
  # Set up iterative variables
  optimal.param <- 0
  stability <- key.stats.df$Stability
  unique.stability <-unique(stability)

  # Determine stability centre
  max.stability <- stability[which.max(stability)]
  stability.max.centre <- unique.stability[-which(unique.stability == stability[1])];
  stability.max.centre <- stability.max.centre[-which(stability.max.centre == stability[40])];
  stability.max.centre <- stability.max.centre[which.max(stability.max.centre)]

  # Default stability minus max
  stability.sub.max <- stability - max.stability

  # Split stability into half
  stability.1 <- stability[1:20]
  stability.2 <- stability[21:40]

  # Check if either end does not change
  initial.check <- c(FALSE, FALSE)
  if (length(unique(stability.1)) == 1){
    initial.check[1] <- TRUE
  }

  if (length(unique(stability.2)) == 1){
    initial.check[2] <- TRUE
  }

  # Changes lie in the first 20
  # Or within the last 20
  # Or somewhere in between
  if (identical(initial.check, c(TRUE, FALSE))){
    if (stability.1[20] %in% stability.2){
      stability.2.idx <- which(stability.2 == stability.1[20])
      if (21 %in% stability.2.idx){
        consecutive.index <- GetConsecutiveSequence(stability.2.idx)
        optimal.param <- stability.2[max(consecutive.index)]
      } else{
        optimal.param <- stability.1[2]
      }
    } else{
      optimal.param <-stability.1[20]
    }
  } else if (identical(initial.check, c(FALSE, TRUE))){
    if (stability.2[1] %in% stability.1){
      stability.idx <- c(which(stability.1 == stability.2[1]))
      if (19 %in% stability.idx){
        consecutive.index <- GetConsecutiveSequence(stability.idx)
        optimal.param <- stability.1[min(consecutive.index)]
      }
      optimal.param <- min(stability.idx)
    } else{
      optimal.param <- stability.2[1]
    }
  } else{
    consecutive <- c()
    restart.positions <- list()
    for (i in 1:length(stability)){
      if (i != length(stability)){
        current.idx <- i
        next.idx <- i + 1
        current.val <- stability[current.idx]
        next.val <- stability[next.idx]
        step.val <- next.val - current.val
        if (step.val == 1){
          if (length(consecutive) > 0){
            if (current.val - consecutive[length(consecutive)] == 1){
              consecutive <- c(consecutive, current.val)
            }
          } else{
            consecutive <- c(consecutive, current.val)
          }
        } else{
          restart.positions[[as.character(current.idx)]] <- length(consecutive)
          consecutive <- c()
        }
      }
    }

    # Get the longest continuous flatline
    longest.seq <- max(unlist(restart.positions))
    endpoint <- as.numeric(names(which(restart.positions == longest.seq)))
    optimal.param <- endpoint - longest.seq
  }
  return(optimal.param)
}

# Generates a bunch of values based on other calculates.
# Called by FindOptimalClusters
BuildKeyStat <- function(rand.idx.matrix){
  # Set up a key stat dataframe
  key.stats.df <- as.data.frame(cbind(as.numeric(rand.idx.matrix$order)*0.025, rand.idx.matrix$stability_count, rand.idx.matrix$cluster.index.ref, rand.idx.matrix$cluster.index.consec))
  colnames(key.stats.df) <-c('Height', 'Stability', 'RandIndex', 'ConsecutiveRI')
  key.stats.df$Height <-as.character(key.stats.df$Height)
  key.stats.df$Cluster_Count <- rand.idx.matrix$cluster_count
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

RetrieveCluster <- function(hclust.obj, dist.mtx, height=0){
  clusters <- unname(dynamicTreeCut::cutreeDynamic(hclust.obj, distM=as.matrix(dist.mtx), minSplitHeight=height, verbose=0))
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
      stop("Please run RunPCA followed by GetReducedDimensions on this object before running this function.")
    }
    if ( is.null(object@PCA$ReducedPCA) ){
      stop("Please reduce the PCA dimensions of this object with GetReducedDimensions before running this function.")
    }
  } else {
    stop("Please supply a AEMSet object.")
  }

  # Generate distance matrix and hclust objects as normal
  print("Performing unsupervised clustering...")
  distance.matrix <- stats::dist(object@PCA$ReducedPCA)
  original.tree <- stats::hclust(distance.matrix, method="ward.D2")
  original.clusters <- unname(dynamicTreeCut::cutreeDynamic(original.tree, distM=as.matrix(distance.matrix), verbose=0))

  # Generate clustering matrix
  # Features number of clusters produced at varying dendrogram cut heights
  # Perform 40 dynamic tree cuts at varying heights
  print("Determining best number of clusters...")
  cluster.list <- BiocParallel::bplapply(seq(0.025:1, by=0.025), function(x) RetrieveCluster(original.tree, distance.matrix, height = x))
  height.list <- BiocParallel::bplapply(seq(0.025:1, by=0.025), function(x) x)
  cluster.matrix <- GenerateClusteringMatrix(original.clusters, cluster.list)

  # Generate RAND index
  rand.idx.matrix <- GenerateRandIndexMatrix(cluster.matrix, original.clusters, cluster.list)

  # Generate Stability Values
  rand.idx.matrix <- GenerateStabilityValues(rand.idx.matrix)

  # Add new cluster counts to rand matrix
  cluster.counts <- BiocParallel::bplapply(cluster.list, function(x) length(unique(x)) )
  rand.idx.matrix$cluster_count <- as.vector(as.numeric(cluster.counts))
  rand.idx.matrix$stability_count <- rand.idx.matrix$stability_count/40

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
    RandMatrix = rand.idx.matrix
  )

  # Append all of this information to the AEMSet object
  object@Clusters <- output.list
  return(object)
}
