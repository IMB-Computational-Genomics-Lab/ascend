################################################################################
#
# ascend_dimreduction.R
# description: Functions related to the analysis of differential expression
#
################################################################################
minMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

bimodLikData_FixNoVar <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- minMax(
    data = length(x = x2) / length(x = x),
    min = 1e-5,
    max = (1 - 1e-5)
  )
  
  likA <- length(x = x1) * log(x = 1 - xal)
  #likelihood of positive cells 
  likB <- length(x = x2) *
    log(x = xal)
  return(likA + likB)
}

bimodLikData <- function(x, xmin = 0) {
  #x1 and x2 are 2 vectors representing 2 modes
  #x1 for 0 values -> on/off distribution model 
  #x2 for positive values -> normal, continuous distribution 
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  #estimate proportion of positive cells 
  #use 1e-5 as min and 1-1e-5 as max (i.e. if there is only 1 nonzero among 100K cells)
  xal <- minMax(
    data = length(x = x2) / length(x = x),
    min = 1e-5,
    max = (1 - 1e-5)
  )
  #likelihood for observing x1, 1-xal is expected ratio of 0 values  
  likA <- length(x = x1) * log(x = 1 - xal)
  #calculate variabce for x2, to be used in dnorm to calculate prob distr
  if (length(x = x2) < 2) {
    mysd <- 1
  } else {
    mysd <- stats::sd(x = x2)
  }
  #likelihood for observing x2 
  likB <- length(x = x2) *
    log(x = xal) +
    sum(stats::dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
  return(likA + likB)
}


differentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimodLikData(x = x)
  lrtY <- bimodLikData(x = y)
  lrtZ <- bimodLikData(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  
  # Check to account for results that do not conform to expected model
  if (is.infinite(lrt_diff) || (lrt_diff < 0) || is.nan(lrt_diff) || is.na(lrt_diff)){
    lrtX <- bimodLikData_FixNoVar(x = x)
    lrtY <- bimodLikData_FixNoVar(x = y)
    lrtZ <- bimodLikData_FixNoVar(x = c(x, y))
    lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  }
  
  return(stats::pchisq(q = lrt_diff, df = 3, lower.tail = FALSE))
}

runGeneLRT <- function(gene, test_matrix = NULL, control_matrix = NULL){
  # Retrieve counts
  test_counts <- as.numeric(unlist(test_matrix[gene, ]))
  control_counts <- as.numeric(unlist(control_matrix[gene, ]))
  
  # Perform Combined LRT
  de_result <- differentialLRT(test_counts, control_counts)
  return(de_result)
}

worker_LRT <- function(x, test_matrix = NULL, control_matrix = NULL){
  lrt_results <- sapply(x, runGeneLRT, test_matrix = test_matrix, control_matrix = control_matrix)
  return(lrt_results)
}

#' runDiffExpression
#' 
#' This function uses a combined Likelihood-Ratio Test (LRT) for 
#' discrete/continuous models to examine differentially expressed genes, on a 
#' gene-gene level. This method was adapted from the method published by [McDavid
#' et al. 2013](https://doi.org/10.1093/bioinformatics/bts714).
#' 
#' @param object An \linkS4class{EMSet} that has undergone filtering and 
#' normalisation
#' @param group A column in colInfo that contains a vector of conditions
#' @param condition.a Condition(s) of first group of cells
#' @param condition.b Condition(s) of second group of cells
#' @param subsampling TRUE or FALSE (Default). Whether or not to subsample 
#' from larger group of cells if cell populations are uneven
#' @param ngenes Test this number of the most variable genes in the dataset
#' 
#' @examples 
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Compare cluster 1 vs cluster 2
#' cluster1_vs_cluster2 <- runDiffExpression(em_set, group = "cluster",
#' condition.a = 1, condition.b = 2, subsampling = FALSE, ngenes = 1500)
#' 
#' @return A data frame with DE analysis results
#' @importFrom methods is
#' @importFrom SingleCellExperiment normcounts
#' @importFrom BiocParallel bpvec
#' @importFrom stats p.adjust
#' @export 
#' 
runDiffExpression <- function(object, 
                              group = NULL,
                              condition.a = NULL, 
                              condition.b = NULL,
                              subsampling = FALSE,
                              ngenes = NULL){
  # Check if EMSet
  if (!(is(object, "EMSet"))){
    stop("Please supply an EMSet.")
  } else{
    if(!("normcounts" %in% SummarizedExperiment::assayNames(object))){
      stop("Please normalise your dataset.")
    }
  }
  
  if (is.null(ngenes)){
    ngenes <- nrow(object)
  }
  
  # Retrieve information from EMSet
  col_info <- colInfo(object)
  
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
  a_cells <- condition_df[which(condition_df[, group] %in% condition.a), "cell_barcode"]
  b_cells <- condition_df[which(condition_df[, group] %in% condition.b), "cell_barcode"]
  
  # Get expression matrix
  expression_matrix <- SingleCellExperiment::normcounts(object)
  
  # If user wishes to subsample data, randomly select barcodes from larger list
  if (subsampling){
    if (length(a_cells) > length(b_cells)){
      a_cells <- sample(a_cells, length(b_cells), replace = FALSE)
    }
    else if (length(b_cells) > length(a_cells)){
      b_cells <- sample(b_cells, length(a_cells), replace = FALSE)
    }
  }
  
  # If user has specified a number of genes, select top n most variable genes
  if (ngenes != nrow(expression_matrix)){
    print("Identifying genes to retain...")
    nonzero_genes <- rownames(expression_matrix)[which(Matrix::rowMeans(expression_matrix) > 0)]
    variable_genes <- rownames(expression_matrix)[which(apply(expression_matrix, 1, stats::sd) > 0)]
    top_genes <- SummarizedExperiment::rowData(object)[order(SummarizedExperiment::rowData(object)$qc_topgeneranking),1][1:ngenes]
    gene_list <- dplyr::intersect(dplyr::intersect(top_genes, nonzero_genes), variable_genes)
  } else{
    gene_list <- rownames(expression_matrix)
  }
  
  # Separate the matrices
  a_matrix <- expression_matrix[gene_list , a_cells]
  b_matrix <- expression_matrix[gene_list, b_cells]
  
  print("Running LRT...")
  lrt_results <- BiocParallel::bpvec(rownames(a_matrix), worker_LRT, test_matrix = a_matrix, control_matrix = b_matrix)
  
  # Format output
  # Create a row of results with the following information:
  # 1. Raw LRT Pvalue (unadjusted, generated by LRT)
  # 2. Mean Test Expression (mean expression of genes in guide-assigned cells)
  # 3. Mean Control Expression (mean expression of genes in control-assigned cells)
  
  print("LRT complete! Returning results...")
  de_result <- data.frame(a_mean = Matrix::rowMeans(a_matrix),
                          b_mean = Matrix::rowMeans(b_matrix),
                          pval = lrt_results,
                          padj = stats::p.adjust(lrt_results, method = "BH"),
                          foldChange = Matrix::rowMeans(a_matrix)/Matrix::rowMeans(b_matrix),
                          log2FoldChange = log2(Matrix::rowMeans(a_matrix)/Matrix::rowMeans(b_matrix)))
  de_result$id <- rownames(de_result)
  de_result <- de_result[ , c("id", "a_mean", "b_mean", 
                              "pval", "padj",
                              "foldChange", "log2FoldChange")]
  de_result <- de_result[order(abs(de_result$foldChange), decreasing = TRUE), ]
  rownames(de_result) <- NULL
  return(de_result)
}
