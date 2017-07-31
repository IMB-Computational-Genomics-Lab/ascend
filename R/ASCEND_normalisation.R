#' scranNormalise
#'
#' Normalise an \linkS4class{AEMSet} with \pkg{scran}'s deconvolution method by Lun et al. 2016.
#'
#' @param object An \linkS4class{AEMSet} that has not undergone normalisation.
#'
scranNormalise <- function(object){
  if(!is.null(object@Log$NormalisationMethod)){
    stop("This data is already normalised.")
  }

  # Convert a scryeR object to a SCRAN object
  print("Converting AEMSet to SCESet...")
  sce.obj <- ConvertToScater(object)

  # Remove controls from SCESet object
  mt.idx <- match(object@Controls$Mt, rownames(sce.obj))
  rb.idx <- match(object@Controls$Rb, rownames(sce.obj))

  # Remove ERCC sequences if in dataset
  ercc.idx <- grep("^ercc", ignore.case=TRUE, rownames(object@ExpressionMatrix))
  if (length(ercc.idx) == 0){
    ercc.idx <- c()
  }
  
  # Remove these items
  remove.idx <- c(mt.idx, rb.idx, ercc.idx)
  sce.obj_rmMtRb <- sce.obj[-remove.idx,]

  # Run computeSumFactors based on number of samples
  # If we have more than 10,000 samples, we will run quickCluster
  # Otherwise we supply a set list of sizes
  if (ncol(sce.obj) > 10000){
    print(sprintf("%i cells detected. Running quickCluster to feed into computeSumFactors...", ncol(sce.obj)))
    quick.cluster <- scran::quickCluster(sce.obj_rmMtRb)
    factored.sce.obj <- scran::computeSumFactors(sce.obj_rmMtRb, clusters = quick.cluster, positive = T)
  } else {
    print(sprintf("%i cells detected. Running computeSumFactors with preset sizes of 40, 60, 80, 100...", ncol(sce.obj)))
    preset.sizes <- c(40, 60, 80, 100)
    factored.sce.obj <- scran::computeSumFactors(sce.obj_rmMtRb, sizes = preset.sizes, positive=T)
  }

  print("scran's computeSumFactors complete. Removing zero sum factors from dataset...")
  zero.size.factors <- which(factored.sce.obj@phenoData@data$size_factor == 0)
  if (any(zero.size.factors)){
    min.size.factor <- min(factored.sce.obj@phenoData@data$size_factor[-zero.size.factors])
    factored.sce.obj@phenoData@data$size_factor[zero.size.factors] <- min.size.factor
  }

  print("Running scater's normalize method...")
  dcvl.sce.obj <- scater::normalize(factored.sce.obj)

  # Convert log-transformed results back into counts
  print("Normalisation complete. Converting SCESet back to AEMSet...")
  dcvl.matrix <- as.matrix(scater::norm_exprs(dcvl.sce.obj))
  unlog.dcvl.matrix <- UnLog2Matrix(dcvl.matrix)
  normalised.obj <- ReplaceExpressionMatrix(unlog.dcvl.matrix, object)
  normalised.obj@Log <- c(normalised.obj@Log, list(NormalisationMethod = "scranNormalise"))

  remove(sce.obj)
  return(normalised.obj)
}

NormWithinBatch <- function(batch.id, expression.matrix = NULL, batch.list = NULL){
  # Function called by NormaliseBatches
  barcodes <- names(batch.list[batch.list == batch.id])
  sub.mtx <- expression.matrix[,barcodes]

  # Collapse all cells into one cell
  collapsed.mtx <- rowSums(sub.mtx)
  norm.factor <- sum(collapsed.mtx)

  ## Scale sub-matrix to normalisation factor
  cpm.mtx <- Matrix::t(Matrix::t(sub.mtx)/norm.factor)

  # Output to list for further calculations
  output.list <- list(NormalisationFactor = norm.factor, CpmMatrix = cpm.mtx)
  return(output.list)
}

#' NormaliseBatches
#'
#' Normalise counts to remove batch effects. This normalisation method is for
#' experiments where data from batches of different samples are combined without
#' undergoing library equalisation.
#'
#' This step must be done prior to any filtering and normalisation between cells.
#'
#' @param object An \linkS4class{AEMSet} object.
#'
NormaliseBatches <- function(object){
  if(!is.null(object@Log$NormalisationMethod)){
    stop("This data is already normalised.")
  }

  # Retrieve variables from AEMSet object
  exprs.mtx <- as.matrix(object@ExpressionMatrix)
  batch.list <- object@BatchInformation
  unique.batch.identifiers <- unique(batch.list)

  # Loop to get batch-specific data
  # PARALLEL
  print("Retrieving batch-specific data...")
  batch.data <- BiocParallel::bplapply(unique.batch.identifiers, NormWithinBatch, expression.matrix = exprs.mtx, batch.list = batch.list)

  # Unpacking results
  print("Scaling data...")
  norm.factors <- unlist(lapply(batch.data, function(x) x$NormalisationFactor))
  median.size <- median(norm.factors)
  sub.matrix.list <- lapply(batch.data, function(x) x$CpmMatrix)
  cpm.matrix <- data.frame(sub.matrix.list)
  scaled.matrix <- cpm.matrix * median.size

  # Load back into AEMSet and write metrics
  print("Returning object...")
  colnames(scaled.matrix) <- names(batch.list)
  object@ExpressionMatrix <- LoadSparseMatrix(scaled.matrix)
  object@Log <- c(object@Log, list(NormaliseBatches=TRUE))
  updated.object <- GenerateMetrics(object)
  return(updated.object)
}

#' NormaliseLibSize
#'
#' Normalise library sizes by scaling.
#' @param object A \linkS4class{AEMSet} object. Please remove spike-ins from the expression matrix before normalising.
#'
NormaliseLibSize <- function(object){
  expression.matrix <- as.matrix(object@ExpressionMatrix)
  norm.factor <- colSums(expression.matrix)
  median.size <- median(norm.factor)

  unscaled.matrix <- Matrix::t(Matrix::t(expression.matrix)/norm.factor)
  normalised.exprs.mtx <- unscaled.matrix * median.size
  object@Log <- c(object@Log, list(NormaliseLibSize=TRUE))
  object@ExpressionMatrix <- LoadSparseMatrix(normalised.exprs.mtx)
  new.object <- GenerateMetrics(object)
  new.object@Log <- c(object@Log, list(NormalisationMethod="NormaliseLibSize"))
  return(new.object)
}


## SUBFUNCTIONS FOR NormaliseByRLE
CalcNormFactor <- function(x, geo.means){
  x.geo.means <- cbind(x, geo.means);
  x.geo.means <-x.geo.means[(x.geo.means[,1] >0),]
  non.zero.median <- median(apply(x.geo.means,1, function(y) { y <-as.vector(y); y[1]/y[2] }))
  return(non.zero.median)
}

CalcGeoMeans <- function(x){
  x <- x [ x > 0]
  x <- exp(mean(log(x)))
  return(x)
}

#' NormaliseByRLE
#'
#' Normalisation of expression between cells, by scaling to relative log expression.
#' This method assumes all genes express a pseudo value higher than 0, and also assumes most genes are not differentially expressed. Only counts that are greater than zero are considered in this normalisation method.
#' @param object A \linkS4class{AEMSet} object that has undergone filtering. Please ensure spike-ins have been removed before using this function.
#'
NormaliseByRLE <- function(object){
  if(!is.null(object@Log$NormalisationMethod)){
    stop("This data is already normalised.")
  }
  expression.matrix <- as.matrix(object@ExpressionMatrix)

  print("Calculating geometric means...")
  geo.means <- apply(expression.matrix, 1, CalcGeoMeans)

  print("Calculating normalisation factors...")
  norm.factor <- apply(expression.matrix, 2, function(x) CalcNormFactor(x, geo.means))

  print("Normalising data...")
  normalised.matrix <- Matrix::t(Matrix::t(expression.matrix)/norm.factor)
  object@ExpressionMatrix <- LoadSparseMatrix(normalised.matrix)
  object@Log <- c(object@Log, list(NormalisationMethod="NormaliseByRLE"))

  return(object)
}
