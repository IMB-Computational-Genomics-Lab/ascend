################################################################################
#
# ascend_DESeq.R
# description: Wrappers for DESeq
#
################################################################################

#' @export
func_DESeq <- function(x, expression_matrix = NULL, condition_df = NULL, group = NULL,
                       condition.a = NULL, condition.b = NULL,
                       fitType = NULL, method = NULL){
  subset_matrix <- expression_matrix[x, ]
  
  # Merge conditions into one if there are more than one
  if (length(condition.a) > 1){
    replace_a_idx <- which(condition_df[, group] %in% condition.a)
    condition.a <- paste0(condition.a, collapse = ",")
    condition_df[replace_a_idx, group] <- condition.a
  }
  
  if(length(condition.b) > 1){
    replace_b_idx <- which(condition_df[, group] %in% condition.b)
    if (length(condition.b) > 2){
      new_condition.b <- paste(condition.b[1:(length(condition.b) - 1)], collapse = ",")
      condition.b <- paste0(new_condition.b, " & ", condition.b[length(condition.b)])
    } else{
      condition.b <- paste(condition.b, collapse = " & ")
    }
    condition_df[replace_b_idx, group] <- condition.b
  }
  
  condition_list <- factor(condition_df[, group], levels = c(condition.a, condition.b))
  # Have to load DESeq as there appears to be an issue with the locfit export...
  library(DESeq)
  count_dataset <- DESeq::newCountDataSet(subset_matrix, condition_list)
  count_dataset <- DESeq::estimateSizeFactors(count_dataset)
  dispersions <- DESeq::estimateDispersions(count_dataset, method = method, fitType = fitType)
  de_seq_results <- DESeq::nbinomTest(dispersions, condition.a, condition.b)
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
#' identified by \code{\link{runCORE}}.
#' @param condition.a Condition of the group you want to use as the baseline.
#' @param condition.b Conditions of the group you want to compare to the baseline.
#' @param ngenes Perform differential expression analysis using top number of genes.
#' If omitted, this function will run analysis on ALL genes.
#' @param fitType Method used to fit a dispersion-mean relation by \pkg{DESeq}. 
#' Options: parametric, local (Default).
#' @param method Method used by \pkg{DESeq} to compute emperical dispersion.
#' Options: pooled, pooled-CR, per-condition (Default), blind.
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
                     method = c("pooled", "pooled-CR", "per-condition", "blind")){
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
  
  # Retrieve data for each condition
  condition_df <- col_info[, c("cell_barcode", group)]
  
  # Group data into conditions
  a_df <- condition_df[which(condition_df[, group] %in% condition.a), ]
  b_df <- condition_df[which(condition_df[, group] %in% condition.b), ]
  
  # Conver to factor
  a_df[, group] <- factor(a_df[, group], levels = unique(as.vector(a_df[, group])))
  b_df[, group] <- factor(b_df[, group], levels = unique(as.vector(b_df[, group])))
  
  # Get expression matrix
  expression_matrix <- SingleCellExperiment::normcounts(object)

  print("Identifying genes to retain...")
  top_genes <- row_data[order(row_data$qc_topgeneranking),1][1:ngenes]
  nonzero_genes <- rownames(expression_matrix)[which(rowMeans(expression_matrix) > 0)]
  variable_genes <- rownames(expression_matrix)[which(apply(expression_matrix, 1, stats::sd) > 0)]
  gene_list <- dplyr::intersect(dplyr::intersect(top_genes, nonzero_genes), variable_genes)
  expression_matrix <- round(expression_matrix + 1)
  
  # Chunk expression matrix by chunks that are equal in size to the nunmber of
  # workers
  nworkers <- BiocParallel::bpnworkers(BiocParallel::bpparam())
  chunked_gene_list <- split(gene_list, 1:nworkers)
  
  print("Running DESeq...")
  de_list <- BiocParallel::bplapply(chunked_gene_list, func_DESeq,
                                    expression_matrix = expression_matrix,
                                    condition_df = condition_df,
                                    group = group,
                                    condition.a = condition.a,
                                    condition.b = condition.b,
                                    fitType = fitType,
                                    method = method)
  
  print("DESeq complete! Adjusting results...")
  de_results <- dplyr::bind_rows(de_list)
  
  adjusted_foldchange <- (de_results$baseMeanB-1)/(de_results$baseMeanA - 1)
  log2_adjustedfoldchange <- log2(adjusted_foldchange)
  de_results$foldChange <- adjusted_foldchange
  de_results$log2FoldChange <- -log2_adjustedfoldchange
  de_results <- de_results[order(de_results$foldChange, decreasing = TRUE), ]
  return(de_results)  
} 
