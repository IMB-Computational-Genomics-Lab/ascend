################################################################################
#
# ascend_DESeq.R
# description: Wrappers for DESeq
#
################################################################################

#' func_DESeq
#' 
#' Function called by parallel version of runDESeq.
#' 
#' @param x List of genes
#' @param expression_matrix Matrix of counts with a pseudocount added
#' @param condition_list Vector of conditions linked to cells
#' @param condition.a Condition A to test
#' @param condition.b Condition B to test
#' @param fitType DESeq argument - fitType to use
#' @param method DESeq argument - method to use
#' @export
func_DESeq <- function(x, 
                       expression_matrix = NULL, 
                       condition_list = NULL, 
                       condition.a = NULL, 
                       condition.b = NULL,
                       fitType = NULL,
                       method = NULL){
  loadNamespace("DESeq")
  loadNamespace("locfit")
  subset_matrix <- expression_matrix[x, ]
  count_dataset <- DESeq::newCountDataSet(subset_matrix, conditions = condition_list)
  count_dataset <- DESeq::estimateSizeFactors(count_dataset)
  dispersions <- DESeq::estimateDispersions(count_dataset, method = method, fitType = fitType)
  de_seq_results <- DESeq::nbinomTest(dispersions, condA = condition.a, condB = condition.b)
  return(de_seq_results)
}

#' runDESeq
#' 
#' This wrapper runs differential expression analysis via the \pkg{DESeq} 
#' package.
#' 
#' @param object An \linkS4class{EMSet} object that has undergone 
#' filtering and normalisation.
#' @param group Name of the column in the colInfo dataframe where you have
#' defined the conditions you would like to test. eg cluster to compare clusters
#' identified by \code{\link[ascend:runCORE]{runCORE}}.
#' @param condition.a Condition of the group you want to use as the baseline.
#' @param condition.b Conditions of the group you want to compare to the baseline.
#' @param ngenes Perform differential expression analysis using top number of genes.
#' If omitted, this function will run analysis on ALL genes.
#' @param fitType Method used to fit a dispersion-mean relation by \pkg{DESeq}. 
#' Options: parametric, local (Default).
#' @param method Method used by \pkg{DESeq} to compute emperical dispersion.
#' Options: pooled, pooled-CR, per-condition (Default), blind.
#' @param parallel Run DESeq through parallelised wrapper (Default: TRUE)
#' @return A dataframe containing \pkg{DESeq} results
#' @examples
#' \dontrun{
#' cluster1_vs_others <- runDESeq(EMSet, group = "cluster", condition.a = "1",
#' condition.b = c("2", "3"), ngenes = 1500, fitType = "local", method = "per-condition")
#' }
#' @importFrom BiocParallel bplapply bpnworkers bpparam
#' @importFrom dplyr intersect bind_rows
#' @importFrom SingleCellExperiment rowData normcounts
#' @export
#'
runDESeq <- function(object, group = NULL, condition.a = NULL,
                     condition.b = NULL, ngenes = NULL,
                     fitType = c("parametric", "local"), 
                     method = c("pooled", "pooled-CR", "per-condition", "blind"),
                     parallel = TRUE){
  if (missing(fitType)){
    fitType <- "local"
  }
  
  if (missing(method)){
    method <- "per-condition"
  }
  
  if (missing(ngenes)){
    ngenes <- nrow(object)  
  }
  
  # Check if EMSet
  if (!(is(object, "EMSet"))){
    stop("Please supply an EMSet.")
  }
  
  # Retrieve information from EMSet
  col_info <- colInfo(object)
  row_data <- SingleCellExperiment::rowData(object)
  
  # Check if group exists in col_info
  if (group %in% colnames(col_info)){
    if (!all(c(condition.a, condition.b) %in% col_info[, group])){
      stop(sprintf("Not all conditions are present in %s. Please check the data in 
                   colInfo and try again.", group))
    }
    } else{
      print(sprintf("%s not found in colInfo. Please choose another group.", group))
  }
  
  # Get condition information
  condition_df <- col_info[, c("cell_barcode", group)]
  
  # Separate into conditions a and conditions b
  a_df <- condition_df[which(condition_df[, group] %in% condition.a), ]
  b_df <- condition_df[which(condition_df[, group] %in% condition.b), ]
  
  
  # Order in input order...
  if (length(condition.a) > 1){
    a_df <- a_df[order(match(a_df[, group], condition.a)), ]  
  }
  
  if (length(condition.b) > 1){
    b_df <- b_df[order(match(b_df[, group], condition.b)), ]  
  }
  
  # Set order for vectors
  cell_barcodes <- c(a_df$cell_barcode, b_df$cell_barcode)
  condition_list <- c(a_df[, group], b_df[, group])
  
  # Extract expression matrix
  expression_matrix <- SingleCellExperiment::normcounts(object)
  expression_matrix <- expression_matrix[, cell_barcodes]
  
  # Calculate most variable genes
  print("Identifying genes to retain...")
  top_genes <- row_data[order(row_data$qc_topgeneranking),1][1:ngenes]
  nonzero_genes <- rownames(expression_matrix)[which(Matrix::rowMeans(expression_matrix) > 0)]
  variable_genes <- rownames(expression_matrix)[which(apply(expression_matrix, 1, stats::sd) > 0)]
  gene_list <- dplyr::intersect(dplyr::intersect(top_genes, nonzero_genes), variable_genes)
  
  # Add pseudo-count
  expression_matrix <- round(expression_matrix + 1)
  
  # Chunk expression matrix by chunks that are equal in size to the nunmber of
  # workers
  nworkers <- BiocParallel::bpnworkers(BiocParallel::bpparam())
  chunked_gene_list <- split(gene_list, 1:nworkers)
  
  # Reformat conditions into a readable string
  reformatCondition <- function(y, condition_list = NULL){
    replace_idx <- which(condition_list %in% y)
    
    if (length(y) > 2){
      new_condition <- paste(y[1:(length(y) - 1)], collapse = ",")
      condition <- paste0(new_condition, " & ", y[length(y)])
    } else{
      condition <- paste(y, collapse = " & ")
    }
    
    condition_list[replace_idx] <- condition
    return(list(condition = condition, condition_list = condition_list))
  }
  
  if (length(condition.a) > 1){
    reformatted <- reformatCondition(condition.a, condition_list = condition_list)  
    condition.b <- reformatted$condition
    condition_list <- reformatted$condition_list
  } else{
    replace_idx <- which(condition_list %in% condition.a)
    condition.a <- as.character(condition.a)
    condition_list[replace_idx] <- condition.a
  }
  
  if (length(condition.b > 1)){
    reformatted <- reformatCondition(condition.b, condition_list = condition_list)
    condition.b <- reformatted$condition
    condition_list <- reformatted$condition_list
  } else{
    replace_idx <- which(condition_list %in% condition.b)
    condition.b <- as.character(condition.b)
    condition_list[replace_idx] <- condition.b
  }
  
  # Factor condition list
  condition_list <- factor(condition_list, levels = c(condition.a, condition.b))
  
  print("Running DESeq...")
  if (parallel){
    de_list <- bplapply(chunked_gene_list, func_DESeq, 
                        expression_matrix = expression_matrix,
                        condition_list = condition_list,
                        condition.a = condition.a,
                        condition.b = condition.b,
                        fitType = fitType,
                        method = method)
    
    print("DESeq complete! Adjusting results...")
    de_results <- dplyr::bind_rows(de_list)    
  } else{
    loadNamespace("DESeq")
    loadNamespace("locfit")
    subset_matrix <- expression_matrix[gene_list, ]
    count_dataset <- DESeq::newCountDataSet(subset_matrix, conditions = condition_list)
    count_dataset <- DESeq::estimateSizeFactors(count_dataset)
    dispersions <- DESeq::estimateDispersions(count_dataset, method = method, fitType = fitType)
    de_results <- DESeq::nbinomTest(dispersions, condA = condition.a, condB = condition.b)
  }
  
  adjusted_foldchange <- (de_results$baseMeanB-1)/(de_results$baseMeanA - 1)
  log2_adjustedfoldchange <- log2(adjusted_foldchange)
  de_results$foldChange <- adjusted_foldchange
  de_results$log2FoldChange <- -log2_adjustedfoldchange
  de_results <- de_results[order(de_results$pval, decreasing = FALSE), ]
  return(de_results)  
} 
