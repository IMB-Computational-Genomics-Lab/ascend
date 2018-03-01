#' RunTSNE
#' 
#' Wrapper for the \code{\link[Rtsne]{Rtsne}} function. Users may call this 
#' directly or call it through the \code{\link{PlotTSNE}} function and also pass
#' additional arguments related to the \code{\link[Rtsne]{Rtsne}} function.
#' 
#' @param object An expression matrix or a PCA-reduced matrix.
#' @param PCA Set this PCA flag to TRUE if the object is a PCA-reduced matrix. 
#' Default: FALSE.
#' @param dimensions Number of dimensions you would like to reduce to. 
#' Default: 2.
#' @param seed (Optional) Set to a specific value for reproducible TSNE plots. 
#' Default: 0.
#' @param perplexity (Optional) Numeric; perplexity parameter. Default: 30.
#' @param theta (Optional) Numeric; Speed/accuracy trade-off. 
#' (increase for less accuracy). Default: 0.5.
#' @param ... Additional arguments to pass on to \code{\link[Rtsne]{Rtsne}}
#' function.
#' @return A dataframe containing expression data for each cell reduced to 
#' selected number of dimensions.
#' @examples
#' \dontrun{
#' tsne_matrix <- RunTSNE(em.set, PCA = TRUE, dimensions = 2, seed = 1, 
#' perplexity = 30, theta = 0.5)
#' }
#' @importFrom methods is
#' @importFrom Rtsne Rtsne
#' @export
#'
RunTSNE <- function(object, PCA = FALSE, dimensions = 2, seed = 0, perplexity = 30, theta = 0.5, ...) {
    if (class(object) == "EMSet") {
        if (PCA) {
            if (is.null(object@PCA$PCA)) {
                x <- object@PCA$PCA
                PCA <- TRUE
            } else {
                x <- object@ExpressionMatrix
                PCA <- FALSE
            }
        }
    } else {
        if (!(is.matrix(expression.matrix) || is.data.frame(expression.matrix) || is(expression.matrix, "sparseMatrix"))) {
            stop("Please supply an EMSet, matrix or dataframe.")
        } else {
            x <- object
            PCA = FALSE
        }
    }
    
    # Will load a matrix to work with, regardless
    if (PCA) {
        transposed.matrix <- as.matrix(x)
    } else {
        transposed.matrix <- as.matrix(Matrix::t(x))
    }
    
    print("Running Rtsne...")
    set.seed(seed)
    tsne <- Rtsne::Rtsne(transposed.matrix, dims = dimensions, pca = PCA, perplexity = perplexity, theta = theta, ignore_duplicates = TRUE, ...)
    tsne.mtx <- data.frame(tsne$Y)
    print("Rtsne complete! Returning matrix...")
    # Add cell names back to the results
    if (PCA) {
        rownames(tsne.mtx) <- rownames(x)
    } else {
        rownames(tsne.mtx) <- colnames(x)
    }
    return(tsne.mtx)
}

#' GetReducedDimensions
#'
#' @param object An \code{\linkS4class{EMSet}} object that has undergone PCA reduction.
#' @param n The number of PC dimensions you would like to select. Refer to
#' vignette on how to select this value. Default: 10.
#' @return An \code{\linkS4class{EMSet}} with a PCA matrix of dimensions ncells by 
#' ndimensions.
#' @examples
#' \dontrun{
#' # Reduce a dataset with RunPCA
#' pca_set <- RunPCA(em.set)
#' 
#' # View scree plot
#' print(PlotPCAVariance(pca_set))
#' 
#' # Reduce number of Principle Components stored in EMSet
#' reduced_pca <- ReduceDimensions(pca_set, n = 10)
#' 
#' }
#' @export
#'
ReduceDimensions <- function(object, n = 10) {
    if (length(object@PCA) == 0) {
        stop("Please run RunPCA on this object before using this function.")
    }
    
    pca.matrix <- object@PCA$PCA
    reduced.pca <- as.data.frame(pca.matrix[, 1:n])
    object@PCA$PCA <- reduced.pca
    return(object)
}

# Functions for PCA
#' CalcColVariance
#' 
#' Internal function called by RunPCA function
#' 
#' @param x A matrix to calculate the column-wise variance of.
#' @return A vector of numeric values.
#' 
CalcColVariance <- function(x) {
    col.variance <- sqrt(colSums((x - colMeans(x))^2)/(dim(x)[1] - 1))
    return(col.variance)
}

#' CalcRowVariance
#' 
#' Internal function called by RunPCA function
#' 
#' @param x A matrix to calculate the row-wise variance of.
#' @return A vector of numeric values.
#' 
CalcRowVariance <- function(x) {
    row.variance <- sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
    return(row.variance)
}

#' RunPCA
#'
#' Reduce the dimensions of an expression matrix stored in an 
#' \code{\linkS4class{EMSet}}. This should be used prior to clustering.
#' 
#' @param object An \code{\linkS4class{EMSet}} that has undergone filtering and 
#' normalisation.
#' @param ngenes The top number of genes you would like to perform the reduction by. 
#' Default: 1500.
#' @param scaling Boolean - set to FALSE if you do not want to scale your values. 
#' Default: TRUE.
#' @return An \code{\linkS4class{EMSet}} with a PCA-reduced matrix stored in the PCA
#' slot.
#' @examples
#' \dontrun{
#' pca_set <- RunPCA(em.set)
#' }
#' @importFrom stats prcomp
#' @export
#'
RunPCA <- function(object, ngenes = 1500, scaling = TRUE) {
    print("Retrieving data...")
    expression.matrix <- GetExpressionMatrix(object, format = "matrix")
    
    # Check if ngenes is too large, just use all the genes to stop it from
    # spazzing out.
    if (ngenes > nrow(expression.matrix)){
      ngenes <- nrow(expression.matrix)
    }
    
    # Selecting top ngenes by variance
    print("Calculating variance...")
    gene.variance <- CalcRowVariance(expression.matrix)
    names(gene.variance) <- rownames(expression.matrix)
    sorted.gene.variance <- gene.variance[order(unlist(gene.variance), decreasing = TRUE)]
    top.genes <- sorted.gene.variance[1:ngenes]
    
    # Subsetting matrix
    subset.matrix <- expression.matrix[names(top.genes), ]
    
    # Skip additional scaling as this was already done with RLE
    scaled.transposed.matrix <- scale(t(subset.matrix), scale = scaling)
    pca.input.matrix <- scaled.transposed.matrix[, CalcColVariance(scaled.transposed.matrix) > 0]
    
    # Compute PCA
    print("Computing PCA values...")
    pca.result <- stats::prcomp(pca.input.matrix)
    pca.percent.var <- pca.result$sdev^2/sum(pca.result$sdev^2)
    
    # Output back to EMSet
    print("PCA complete! Returning object...")
    pca.result.matrix <- as.data.frame(pca.result$x)
    object@PCA <- list(PCA = pca.result.matrix, PCAPercentVariance = pca.percent.var)
    object@Log <- c(object@Log, list(PCA = TRUE))
    return(object)
}
