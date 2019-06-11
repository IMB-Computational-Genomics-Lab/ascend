################################################################################
#
# ascend_dimreduction.R
# description: Functions related to the dimensional reduction of data.
#
################################################################################

#' runUMAP
#' 
#' Wrapper for the \code{\link[umap]{umap}} function. 
#' 
#' @param object An EMSet
#' @param method Method to use with the UMAP function - naive (default) or 
#' umap-learn (Python package required)
#' @param config Configuration to use with UMAP function (Optional) 
#' @param ... Additional arguments to pass to the UMAP function
#' @return UMAP matrix stored in "UMAP" slot in reducedDims
#' 
#' @examples 
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Run UMAP
#' em_set <- runUMAP(em_set, method = "naive")
#' 
#' @importFrom umap umap
#' @importFrom SummarizedExperiment assayNames
#' @export
setGeneric("runUMAP", def = function(object, ..., method, config){
  standardGeneric("runUMAP")
})


#' @rdname runUMAP
setMethod("runUMAP", signature = "EMSet", function(object,
                                                   ...,
                                                   method = c("naive", "umap-learn"),
                                                   config = NULL){
  loadNamespace("umap")
    # Check if normcounts are in arrayNames
  if (!("normcounts" %in% SummarizedExperiment::assayNames(object))){
    stop("Please normalise your data before proceeding.")
  }
  
  if (missing(method)){
    method <- "naive"
  }
  
  if (missing(config)){
    config <- umap::umap.defaults
  }
  
  # Run UMAP
  umap_obj <- umap::umap(as.matrix(Matrix::t(normcounts(object))), method = method, config = config, ...)
  
  # Get generated values from UMAP object
  umap_matrix <- umap_obj$layout
  
  # Store in object
  reducedDim(object, "UMAP") <- umap_matrix
  log <- progressLog(object)
  log$UMAP <- umap_obj$config
  progressLog(object) <- log
  return(object)
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
#' @param dims Number of dimensions you would like to reduce to. 
#' Default: 2.
#' @param seed (Optional) Set to a specific value for reproducible TSNE plots. 
#' Default: 0.
#' @param perplexity (Optional) Numeric; perplexity parameter. Default: 30.
#' @param theta (Optional) Numeric; Speed/accuracy trade-off. 
#' (increase for less accuracy). Default: 0.5.
#' @return t-SNE coordinates in the "TSNE" slot of reducedDims
#' 
#' @examples
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Run tSNE
#' em_set <- runTSNE(em_set, PCA = TRUE, dims = 2, seed = 1)
#' 
#' @include ascend_objects.R
#' @importFrom Rtsne Rtsne
#' @importFrom SingleCellExperiment reducedDimNames reducedDims normcounts
#' @importFrom methods is
#' @export
setGeneric("runTSNE", def = function(object, ..., PCA, dims, seed, 
                                     perplexity, theta) {
  standardGeneric("runTSNE")  
})

#' @rdname runTSNE
setMethod("runTSNE", signature("EMSet"), function(object, ...,
                                                  PCA = FALSE, 
                                                  dims = 2, 
                                                  seed = 0, 
                                                  perplexity = 30, 
                                                  theta = 0.5
                                                  ){
  # Fill in missing values
  if (missing(PCA)){
    PCA <- FALSE
  }
  if (missing(dims)){
    dims <- 2
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

  raw_matrix <- Matrix::as.matrix(raw_matrix)
  
  print("Running Rtsne...")
  set.seed(seed)
  loadNamespace("Rtsne")
  tsne <- Rtsne::Rtsne(raw_matrix, dims = dims, 
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

calcVariance <- function(x, axis = c("row", "column")){
  direction <- list(row = 1, column = 2)
  variance <- apply(x, direction[[axis]], stats::var)
  return(variance)
}

#' runPCA
#'
#' Reduce the dimensions of an expression matrix stored in an 
#' \linkS4class{EMSet} based on the most variable genes. Datasets must
#' be reduced prior to clustering analysis as it is used to construct the 
#' distance matrix.
#'  
#' @param object An \linkS4class{EMSet} that has undergone filtering and 
#' normalisation.
#' @param ngenes The top number of genes you would like to perform the reduction 
#' by. Default: 1500.
#' @param scaling Boolean - set to FALSE if you do not want to scale your values. 
#' Default: TRUE.
#' @param ... Additional arguments to pass to irlba
#' @return An \linkS4class{EMSet} with a PCA-reduced matrix stored under PCA in
#' the reducedDims slot
#' 
#' @include ascend_objects.R
#' 
#' @examples 
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Run PCA
#' em_set <- runPCA(em_set, scaling = TRUE, ngenes = 1500)
#' 
#' 
#' @importFrom irlba prcomp_irlba
#' @importFrom SingleCellExperiment normcounts reducedDim
#' @export
setGeneric("runPCA", def = function(object, ..., ngenes, scaling) {
             standardGeneric("runPCA")
           })

#' @rdname runPCA
setMethod("runPCA", signature("EMSet"), function(object,
                                                 ...,
                                                 ngenes = 1500,
                                                 scaling = TRUE){
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
  
  print("Identifying variable genes...")
  matrix <- normcounts(object)
  gene_variance <- calcVariance(matrix, axis = "row")
  sorted_gene_variance <- gene_variance[order(unlist(gene_variance), decreasing = TRUE)]
  top_genes <- sorted_gene_variance[1:ngenes]
  
  # Subset matrix
  matrix <- matrix[names(top_genes), ]
  print("Computing 50 PCs with irlba...")
  loadNamespace("irlba")
  matrix_irlba <- irlba::prcomp_irlba(Matrix::t(matrix), 
                                      n = 50, 
                                      retx = TRUE, 
                                      scale. = scaling,
                                      ...)
  pca_percent_var <- (matrix_irlba$sdev^2/sum(matrix_irlba$sdev^2))
  print("PCA complete! Loading PCA into EMSet...")
  matrix <- as.matrix(matrix_irlba$x)
  
  SingleCellExperiment::reducedDim(object, "PCA") <- matrix
  log <- progressLog(object)
  log$PCAVariance <- pca_percent_var
  progressLog(object) <- log
  
  remove(matrix)
  return(object)
})
