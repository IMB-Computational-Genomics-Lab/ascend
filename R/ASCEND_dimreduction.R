#' GetReducedDimensions
#'
#' @param object An \linkS4class{AEMSet} object that has undergone PCA reduction.
#' @param n The number of PC dimensions you would like to select. Refer to vignette on how to select this value.
#'
ReduceDimensions <- function(object, n = 10){
  pca.matrix <- object@PCA$PCA
  reduced.pca <- as.data.frame(pca.matrix[, 1:n])
  object@PCA <- c(object@PCA, list(ReducedPCA = reduced.pca))
  return(object)
}

# Functions for PCA
CalcColVariance <- function(x){
  col.variance <- sqrt(colSums((x - colMeans(x))^2)/(dim(x)[1] - 1))
  return(col.variance)
}

CalcRowVariance <- function(x){
  row.variance <- sqrt(rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1))
  return(row.variance)
}

#' RunTSNE
#'
#' @param object An \linkS4class{AEMSet} object, or an expression matrix than can either be a full count matrix or a PCA-reduced matrix.
#' @param PCA Set this PCA flag to true if you are using a PCA-reduced matrix.
#' @param dimensions Number of dimensions you would like to reduce to
#'
RunTSNE <- function(object, PCA=FALSE, dimensions = 2){
  if (class(object) == "AEMSet"){
    if (PCA){
      if (!is.null(object@PCA$ReducedPCA)){
        x <- object@PCA$ReducedPCA
        PCA <- TRUE
      } else if(is.null(object@PCA$PCA)){
        x <- object@PCA$PCA
        PCA <- TRUE
      } else{
        x <- object@ExpressionMatrix
        PCA <- FALSE
      }
    }
  } else if(is.matrix(object) || is.data.frame(object)){
    x <- object
  } else{
    stop("Please supply an AEMSet, matrix or dataframe.")
  }
  # Will load a matrix to work with, regardless
  if (PCA){
    transposed.matrix <- as.matrix(x)
  } else{
    transposed.matrix <- as.matrix(Matrix::t(x))
  }

  tsne.3d <- Rtsne::Rtsne(transposed.matrix,dims=dimensions, pca = PCA, ignore_duplicates = TRUE)
  tsne.mtx <- data.frame(tsne.3d$Y)

  # Add cell names back to the results
  if (PCA){
    rownames(tsne.mtx) <- rownames(x)
  } else{
    rownames(tsne.mtx) <- colnames(x)
  }

  return(tsne.mtx)
}

#' RunPCA
#'
#' @param object A \linkS4class{AEMSet} object that has undergone filtering and normalisation.
#' @param ngenes The top number of genes you would like to perform the reduction by.
#' @param scaling Boolean - set to FALSE if you do not want to scale your values
#'
RunPCA <- function(object, ngenes = 1500, scaling = TRUE){
  print("Retrieving data...")
  expression.matrix <- LoadMatrix(object)

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
  pca.percent.var <- pca.result$sdev ^ 2 / sum(pca.result$sdev ^ 2)

  # Output back to AEMSet
  print("PCA complete! Returning object...")
  pca.result.matrix <- as.data.frame(pca.result$x)
  object@PCA <- list(PCA = pca.result.matrix, PCAPercentVariance = pca.percent.var)
  object@Log <- c(object@Log, list(PCA=TRUE))
  return(object)
}

