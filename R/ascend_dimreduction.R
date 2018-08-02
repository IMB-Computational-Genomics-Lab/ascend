################################################################################
#
# ascend_dimreduction.R
# description: Functions related to the dimensional reduction of data.
#
################################################################################

#' @export
setGeneric("runTSNE", def = function(object, ..., PCA, dimensions, seed, 
                                     perplexity, theta) {
  standardGeneric("runTSNE")  
})

#' runTSNE
#' 
#' Wrapper for the \code{\link[Rtsne]{Rtsne}} function. Users may call this 
#' directly or call it through the \code{\link{plotTSNE}} function and also pass
#' additional arguments related to the \code{\link[Rtsne]{Rtsne}} function.
#' 
#' @param object An expression matrix or a PCA-reduced matrix.
#' @param ... Additional arguments to pass on to \code{\link[Rtsne]{Rtsne}}
#' function.
#' @param PCA Set this PCA flag to TRUE if the object is a PCA-reduced matrix. 
#' Default: FALSE.
#' @param dimensions Number of dimensions you would like to reduce to. 
#' Default: 2.
#' @param seed (Optional) Set to a specific value for reproducible TSNE plots. 
#' Default: 0.
#' @param perplexity (Optional) Numeric; perplexity parameter. Default: 30.
#' @param theta (Optional) Numeric; Speed/accuracy trade-off. 
#' (increase for less accuracy). Default: 0.5.

#' @return A dataframe containing expression data for each cell reduced to 
#' selected number of dimensions.
#' 
#' @include ascend_objects.R
#' @importFrom SingleCellExperiment reducedDimNames reducedDims normcounts
#' @importFrom methods is
#' @importFrom Rtsne Rtsne
#' @export
#'
setMethod("runTSNE", signature("EMSet"), function(object, ...,
                                                  PCA = FALSE, 
                                                  dimensions = 2, 
                                                  seed = 0, 
                                                  perplexity = 30, 
                                                  theta = 0.5
                                                  ){
  # Fill in missing values
  if (missing(PCA)){
    PCA <- FALSE
  }
  if (missing(dimensions)){
    dimensions <- 2
  }
  if (missing(seed)){
    seed <- 0
  }
  if (missing(perplexity)){
    perplexity <- 30
  }
  if (missing(theta)){
    theta <- 0.5
  }
  
  # If PCA is true, check PCA matrix has been generated
  if (PCA){
    if (!("PCA" %in% SingleCellExperiment::reducedDimNames(object))){
      stop("Please generate a PCA matrix with runPCA before using this option.")
    } else{
      raw_matrix <- SingleCellExperiment::reducedDim(object, "PCA")
    }
  } else{
    raw_matrix <- Matrix::t(SingleCellExperiment::normcounts(object))
  }

  print("Running Rtsne...")
  set.seed(seed)
  tsne <- Rtsne::Rtsne(raw_matrix, dims = dimensions, 
                       pca = PCA, 
                       perplexity = perplexity, 
                       theta = theta, 
                       ignore_duplicates = TRUE, ...)
  tsne_matrix <- as.matrix(tsne$Y)
  
  print("Rtsne complete! Returning matrix...")
  # Add cell names back to the results
  rownames(tsne_matrix) <- rownames(raw_matrix)
  SingleCellExperiment::reducedDim(object, "TSNE") <- tsne_matrix
  return(object)
})

#' @export
calcVariance <- function(x, axis = c("row", "column")){
  if (axis == "row"){
    variance <- sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
  }
  if (axis == "column"){
    variance <- sqrt(colSums((x - colMeans(x))^2)/(dim(x)[1]-1))
  }
  return(variance)
}

#' @export
setGeneric("runPCA", def = function(object, ...., ngenes, scaling) {
             standardGeneric("runPCA")
           })

#' runPCA
#'
#' Reduce the dimensions of an expression matrix stored in an 
#' \code{\linkS4class{EMSet}} based on the most variable genes. Datasets must
#' be reduced prior to clustering analysis as it is used to construct the 
#' distance matrix.
#'  
#' @param object An \code{\linkS4class{EMSet}} that has undergone filtering and 
#' normalisation.
#' @param ngenes The top number of genes you would like to perform the reduction by. 
#' Default: 1500.
#' @param scaling Boolean - set to FALSE if you do not want to scale your values. 
#' Default: TRUE.
#' @return An \code{\linkS4class{EMSet}} with a PCA-reduced matrix stored in the PCA
#' slot.
#' @include ascend_objects.R
#' @importFrom stats prcomp
#' @importFrom SingleCellExperiment normcounts reducedDim
#' @export
setMethod("runPCA", signature("EMSet"), function(object, ngenes = 1500, scaling = TRUE){
  # Check for ngenes
  if (missing(ngenes)){
    ngenes <- 1500
  }
  if (missing(scaling)){
    scaling <- TRUE
  }
  
  if (ngenes > nrow(object)){
    print(sprintf("There are less genes than specified. Setting number of genes to %i.", nrow(object)))
    ngenes <- nrow(object)
  }
  
  # Check if normalised
  if (is.null(progressLog(object)$NormalisationMethod)){
    stop("Please normalise your dataset before using this function.")
  }
  
  # Extract normalised counts
  expression_matrix <- as.matrix(SingleCellExperiment::normcounts(object))
  gene_variance <- calcVariance(expression_matrix, axis = "row")
  sorted_gene_variance <- gene_variance[order(unlist(gene_variance), decreasing = TRUE)]
  top_genes <- sorted_gene_variance[1:ngenes]
  
  # Subset matrix
  subset_matrix <- expression_matrix[names(top_genes), ]
  scaled_matrix <- scale(t(subset_matrix), scale = scaling)
  #scaled_transposed_matrix <- scale(t(subset_matrix), scale = scaling)
  pca_input_matrix <- scaled_matrix[, calcVariance(scaled_matrix, axis = "column") > 0]
  
  print("Computing PCA values...")
  pca_result <- stats::prcomp(pca_input_matrix)
  pca_percent_var <- pca_result$sdev^2/sum(pca_result$sdev^2)
  
  print("PCA complete! Loading PCA into EMSet...")
  pca_matrix <- as.matrix(pca_result$x)
  SingleCellExperiment::reducedDim(object, "PCA") <- pca_matrix
  
  # Update log
  log <- progressLog(object)
  log$PCAVariance <- pca_percent_var
  progressLog(object) <- log
  
  return(object)
})