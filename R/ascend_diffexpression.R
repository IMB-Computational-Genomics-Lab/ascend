# Code for differential expression Requires DESeq to run, and parallel

#' RunDESeq
#'
# Called by RunDiffExpression in parallel. This performs the differential expression part.
RunDESeq <- function(data, condition.list = list(), condition.a = NULL, condition.b = NULL, fitType = NULL, method = NULL) {
    library(DESeq)
    count.dataset <- DESeq::newCountDataSet(data, condition.list)
    count.dataset <- DESeq::estimateSizeFactors(count.dataset)
    dispersions <- DESeq::estimateDispersions(count.dataset, method = method, fitType = fitType)
    de.seq.results <- DESeq::nbinomTest(dispersions, condition.a, condition.b)
    return(de.seq.results)
}

#' ProcessDEREsults
#'
#' Called by RunDiffExpression. Compiles the resultant data into data frames and converts
#' the results to absolute fold change.
#'
ProcessDEResults <- function(output.list) {
    de.result.df <- dplyr::bind_rows(output.list)

    # Adjust Foldchange
    print("Adjusting fold change values...")
    adjusted.foldchange <- (de.result.df$baseMeanB - 1)/(de.result.df$baseMeanA - 1)
    log2.adjusted.foldchange <- log2(adjusted.foldchange)
    de.result.df$foldChange <- adjusted.foldchange
    de.result.df$log2FoldChange <- log2.adjusted.foldchange
    de.result.df <- de.result.df[order(de.result.df$pval, decreasing = F), ]
    return(de.result.df)
}

#' PrepareCountData
#'
#' Called by RunDiffExpression. This chunks up the expression matrix to feed into DESeq.
#'
PrepareCountData <- function(object, cells) {
    # Retrieve information from EMSet
    expression.matrix <- GetExpressionMatrix(object, "matrix")
    expression.matrix <- expression.matrix[, cells]
    
    # Filter out genes with zero expression, and add one to make it friendly for DESeq
    print("Rounding expression matrix values...")
    expression.matrix <- expression.matrix[which(rowMeans(expression.matrix) > 0), ]
    sd.row <- apply(expression.matrix, 1, sd)
    expression.matrix <- expression.matrix[names(sd.row > 0), ]
    expression.matrix <- round(expression.matrix + 1)

    print("Chunking matrix...")
    # We want chunks with 1K rows
    chunk.size <- nrow(expression.matrix)/1000
    chunked.matrix <- ChunkMatrix(expression.matrix, axis = 0, chunks = chunk.size)
    return(chunked.matrix)
}


#' RunDiffExpression
#'
#' Compare the differential expression of genes in each cluster versus other clusters.
#'
#' @param object A \linkS4class{EMSet} object that has undergone clustering
#' with the \code{\link{RunCORE}} function.
#' @param conditions Name of the column in the CellInformations lot where you have
#' defined the conditions you would like to test. eg cluster to compare clusters
#' identified by RunCORE.
#' @param condition.a Condition of the group you want to use as the baseline
#' @param condition.b Condition of the group you want to compare to the baseline. 
#' @param fitType Method used to fit a dispersion-mean relation by \pkg{DESeq}. 
#' Options: parametric, local (Default)
#' @param method Method used by \pkg{DESeq} to compute emperical dispersion.
#' Options: pooled, pooled-CR, per-condition (Default), blind 
#' @export
#'
RunDiffExpression <- function(object, conditions = NULL, condition.a = NULL, condition.b = NULL, fitType = c("parametric", "local"), method = c("pooled", "pooled-CR", "per-condition", "blind")) {
  # Object check
  if (class(object) != "EMSet") {
    stop("Please supply a EMSet object.")
  }
  
  # If user wants to compare clusters but hasn't run it
  if ((conditions == "cluster") & (is.null(object@CellInformation[, "cluster"]))) {
    stop("Please run the RunCORE function on this object before using this function.")
  }
  
  # Check for missing variables 
  if (missing(conditions) | missing(condition.a) | missing(condition.b)) {
    stop("Please supply your conditions and try again.")
  }
  
  # Check your conditions
  if (!(condition.a %in% object@CellInformation[, conditions])){
    stop("Please make sure Condition A is in your conditions column.")
  }
  if (!(condition.b %in% object@CellInformation[, conditions]) & (condition.b != "Others")){
    stop("Please make sure Condition B is in your conditions column.")
  }
  
  # Use default values for DESeq if user hasn't defined them.
  if (missing(fitType)){
    fitType <- "local"
  }
  if (missing(method)){
    method <- "per-condition"
  }
  
  # Prepare conditions for input into DESeq
  cell.info <- GetCellInfo(object)
  
  # Identify relevent conditions
  barcodes.a <- as.character(cell.info[,1][which(cell.info[, conditions] == condition.a)])
  
  # Subset cells if they match Condition B - Either a specific condition or
  # others
  if (condition.b == "Others"){
    barcodes.b <- as.character(cell.info[,1][which(cell.info[, conditions] != condition.a)])
  } else{
    barcodes.b <- as.character(cell.info[,1][which(cell.info[, conditions] == condition.b)])      
  }
  
  condition.a <- as.character(condition.a)
  condition.b <- as.character(condition.b)
  cells <- c(barcodes.a, barcodes.b)
  condition.list <- cell.info[which(cell.info[,1] %in% cells), conditions]
  
  if (condition.b == "Others"){
    condition.list[which(condition.list != condition.a)] <- "Others"
  }
  
  condition.list <- as.factor(condition.list)
  # Prepare data for differential expression
  print("Processing expression matrix...")
  chunked.matrix <- PrepareCountData(object, cells)
  
  print("Running DESeq...")    
  result.list <- BiocParallel::bplapply(chunked.matrix, RunDESeq, 
                                        condition.list = condition.list, 
                                        condition.a = condition.a, 
                                        condition.b = condition.b, 
                                        fitType = fitType,
                                        method = method)
  
  print("Differential expression complete!")
  
  print("Combining DE results...")
  output <- ProcessDEResults(result.list)
  print(sprintf("Condition: %s vs %s complete!", condition.a, condition.b))
  return(output)
}