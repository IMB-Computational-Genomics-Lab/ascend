################################################################################
#
# ascend_DESeq2.R
# description: Wrappers for DESeq2
#
################################################################################

#' runDESeq2
#' 
#' Wrapper for DESeq2. 
#' 
#' This is a wrapper for DESeq2. The function handles the conversion of the 
#' \linkS4class{EMSet} to a DESeqDataset object for use with the DESeq2 
#' functions.
#' 
#' @param object An \linkS4class{EMSet}
#' @param group Column in colInfo that contains a set of conditions you want to test
#' @param condition.a Condition(s) in the group column that you want to test for.
#' @param condition.b Condition(s) in the group column that you want to test against.
#' @param fitType Type of fit to use with DESeq2 - parametric, local (Default) and mean.
#' @param ngenes Number of genes you want to test (Default: all genes)
#' 
#' @importFrom dplyr intersect
#' @importFrom SingleCellExperiment rowData counts
#' @export
#'
runDESeq2 <- function(object, 
                      group = NULL,
                      condition.a = NULL,
                      condition.b = NULL,
                      fitType = c("parametric", "local", "mean"),
                      ngenes = NULL){
  # Retrieve information from EMSet
  col_info <- colInfo(object)
  row_data <- SingleCellExperiment::rowData(object)
  
  if (is.null(ngenes)){
    ngenes <- nrow(object)
  }
  
  if (missing(fitType)){
    fitType <- "local"
  }
  
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
  
  # Merge conditions into one if there are more than one
  if (length(condition.a) > 1){
    condition.a <- paste0(condition.a, collapse = "_")
    a_df[ , group] <- condition.a
  }
  
  if(length(condition.b) > 1){
    condition.b <- paste(condition.b, collapse = "_")
    b_df[, group] <- condition.b
  }
  
  # Convert grouped conditions
  # Convert to factor
  a_df[, group] <- factor(a_df[, group], levels = unique(as.vector(a_df[, group])))
  b_df[, group] <- factor(b_df[, group], levels = unique(as.vector(b_df[, group])))
  
  # Use the simple way to make the datasets as we want to chunk it...
  condition_df <- S4Vectors::DataFrame(rbind(a_df, b_df))
  barcodes <- condition_df$cell_barcode
  conditions <- condition_df[, group]
  
  # Get expresison matrix
  expression_matrix <- SingleCellExperiment::counts(object)
  expression_matrix <- expression_matrix[, barcodes]
  top_genes <- row_data[match(sort(row_data$qc_topgeneranking), row_data$qc_topgeneranking), "gene_id"][1:ngenes]
  
  print("Identifying genes to retain...")
  nonzero_genes <- rownames(expression_matrix)[which(rowMeans(expression_matrix) > 0)]
  variable_genes <- rownames(expression_matrix)[which(apply(expression_matrix, 1, stats::sd) > 0)]
  gene_list <- dplyr::intersect(dplyr::intersect(top_genes, nonzero_genes), variable_genes)
  expression_matrix <- round(expression_matrix + 1)
  
  # Feed chunks into DESeq
  deseq_obj <- DESeq2::DESeqDataSetFromMatrix(countData = expression_matrix, 
                                              colData = condition_df,
                                              design = as.formula(paste0("~ ", group)))
  
  deseq_result <- DESeq2::DESeq(deseq_obj, 
                                fitType = fitType, 
                                minReplicatesForReplace = Inf, 
                                parallel = TRUE)
  
  # Unpack DESeq2 results
  results <- DESeq2::results(deseq_result, 
                             contrast = c(group, condition.a, condition.b), 
                             parallel = TRUE)
  return(results)
}
  