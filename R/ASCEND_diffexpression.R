# Code for differential expression
# Requires DESeq to run, and parallel

# Called by RunDiffExpression
RunDESeq <- function(data, condition.list = list(), condition.a = NULL, condition.b = NULL){
  library(DESeq)
  count.dataset <- DESeq::newCountDataSet(data, condition.list)
  count.dataset <- DESeq::estimateSizeFactors(count.dataset)
  dispersions <- DESeq::estimateDispersions(count.dataset, method = 'per-condition', fitType = "local")
  de.seq.results <- DESeq::nbinomTest(dispersions, condition.a, condition.b)
  return(de.seq.results)
}

# Called by RunClusterDiffExpression
GenerateConditionList <- function(condition.a = NULL, condition.b = NULL, barcode.list = NULL){
  condition.a.idx <- which(barcode.list == as.character(condition.a))
  condition.b.idx <- which(barcode.list != as.character(condition.a))
  condition.list <- as.vector(barcode.list)
  condition.list[condition.a.idx] <- as.character(condition.a)
  condition.list[condition.b.idx] <- as.character(condition.b)
  return(condition.list)
}

# Called by RunDiffExpression
ProcessDEResults <- function(output.list){
  de.result.df <- dplyr::bind_rows(output.list)

  # Adjust Foldchange
  print("Adjusting fold change values...")
  adjusted.foldchange <- (de.result.df$baseMeanB - 1) / (de.result.df$baseMeanA - 1)
  log2.adjusted.foldchange <- log2(adjusted.foldchange)
  de.result.df$foldChange <- adjusted.foldchange
  de.result.df$log2FoldChange <- log2.adjusted.foldchange
  de.result.df <- de.result.df[order(de.result.df$pval, decreasing = F), ]
  return(de.result.df)
}

# Called by RunDiffExpression
PrepareCountData <- function(x){
  # Set up expression matrix
  if (class(x) == "AEMSet"){
    expression.matrix <- GetExpressionMatrix(x, "matrix")
  } else if(is.data.frame(x)){
    expression.matrix <- as.matrix(x)
  } else if(is.matrix(x)){
    expression.matrix <- x
  } else{
    stop("Please supply an expression matrix in one of the following formats: AEMSet, data.frame, matrix")
  }

  # Filter out genes with zero expression, and add one to make it friendly for DESeq
  print("Rounding expression matrix values...")
  expression.matrix <- expression.matrix[which(rowMeans(expression.matrix) > 0), ]
  expression.matrix <- round(expression.matrix + 1)

  print("Chunking matrix...")
  # We want chunks with 1K rows
  chunk.size <- nrow(expression.matrix) / 1000
  chunked.matrix <- ChunkMatrix(expression.matrix, axis=0, chunks=chunk.size)
  return(chunked.matrix)
}

# Called by RunDiffExpression
VerifyArguments <- function(condition.a = NULL, condition.b = NULL, condition.list = list()){
  # CHECK DATA FIRST
  print("Verifying input...")
  # Check for missing arguments
  if ( missing(condition.a) || missing(condition.list) ){
    stop("Please supply a string for Condition A and Condition B, in addition to a Condition List.")
  }

  # Check there are two conditions in the condition.list, and they match what the user has supplied.
  if (length(unique(condition.list)) != 2){
    stop("Please ensure there are only two conditions in the supplied Condition List.")
  }

  unique.conditions <- unique(condition.list)
  condition.bool <- c((!condition.a %in% unique.conditions), (!condition.b %in% unique.conditions))

  if ( any(condition.bool) ){
    stop("Please ensure both of the specified conditions are present in the Condition List.")
  }

  print("Input is acceptable. Verification complete!")
}

#' RunDiffExpression
#'
#' Runs differential expression on a pre-defined set of conditions. This function can only handle two conditions at a time, where Condition A describes a subset of cells and Condition B describes the remainder.
#'
#' @param x An expression matrix stored in either an \linkS4class{AEMSet} or a data.frame
#' @param condition.a A one-word string describing the main condition
#' @param condition.b A one-word string describing cells that do not fulfil condition A
#' @param condition.list A named list featuring individual cell identifiers and their associated condition.
#'
RunDiffExpression <- function(x, condition.a = NULL, condition.b = "Other", condition.list = NULL){
  # Run Verification
  VerifyArguments(condition.a = condition.a, condition.b = condition.b, condition.list = condition.list)

  # Set up expression matrix
  print("Processing expression matrix...")
  chunked.matrix <- PrepareCountData(x)

  # Convert condition list to factors
  condition.list <- as.factor(unlist(condition.list))

  print(sprintf("Running differential expression on %s vs %s", condition.a, condition.b))

  # Coerce condition variables into characters, so DESeq won't reject it.
  condition.a <- as.character(condition.a)
  condition.b <- as.character(condition.b)


  # Add all of this information into an environment to load into parallel
  print("Running DESeq...")
  result.list <- BiocParallel::bplapply(chunked.matrix, RunDESeq, condition.list = condition.list, condition.a = condition.a, condition.b = condition.b)

  print("Differential expression complete!")
  print("Returning values...")

  print("Combining DE results...")
  de.result.df <- ProcessDEResults(result.list)
  print(sprintf("Condition: %s vs %s complete!", condition.a, condition.b))
  return(de.result.df)
}

#' RunClusterDiffExpression
#'
#' Compare the differential expression of genes in each cluster versus other clusters.
#'
#' @param object A \linkS4class{AEMSet} object that has undergone clustering with the \code{\link{FindOptimalClusters}} function.
#'
RunClusterDiffExpression <- function(object){
  # Object check
  if ( class(object) != "AEMSet" ){
    stop("Please supply a AEMSet object.")
  }
  # Clustering check
  if ( is.null(object@Clusters$Clusters)){
    stop("Please run the FindOptimalClusters function on this object before using this function.")
  }

  # Prepare Clusters
  cluster.list <- as.factor(object@Clusters$Clusters)
  clusters <- sort(unique(cluster.list))
  output <- list()
  if (length(clusters) > 1){
    condition.lists <- list()
    if (length(clusters) > 2){
      condition.lists <- lapply(clusters, function(x) GenerateConditionList(condition.a = x, condition.b = "Others", barcode.list = cluster.list))
      names(condition.lists) <- clusters
    } else{
      condition.lists[[as.character(clusters[1])]] <- cluster.list
    }
    # Loop over condition lists and run DE
    for (x in names(condition.lists)){
      condition.a <- x
      condition.b <- clusters[-which(clusters == x)]
      if (length(condition.b) > 1){
        condition.b <- "Others"
      } else{
        condition.b <- condition.b[1]
      }
      diff.exp <- RunDiffExpression(object, condition.a = x, condition.b = condition.b, condition.list = condition.lists[[x]])
      output[[paste0(as.character(x), "vs", condition.b)]] <- diff.exp
    }
  } else{
    stop("You must have more than one cluster in order to run pairwise comparisons of clusters.")
  }
  return(output)
}
