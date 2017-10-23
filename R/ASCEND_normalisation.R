#' scranNormalise
#'
#' Normalise an \linkS4class{EMSet} with \pkg{scran}'s deconvolution method by Lun et al. 2016.
#'
#' @details Pooling method of cells is influenced by the size of the cell population.
#' For datasets containing less than 20000 cells, this function will run \code{\link[scran]{computeSumFactors}}
#' with preset sizes of 40, 60, 80 and 100.
#' For datasets containing over 20000 cells, you have the option to use \code{\link[scran]{quickCluster}} to
#' form pools of cells to feed into \code{\link[scran]{computeSumFactors}}. This method is slower and may fail
#' with very large datasets containing over 40000 cells. If \code{\link[scran]{quickCluster}} is not selected,
#' cells will be randomly assigned to a group.
#'
#' @param object An \linkS4class{EMSet} that has not undergone normalisation.
#' @param quickCluster TRUE: Use scran's quickCluster method FALSE: Use randomly-assigned groups. Default: FALSE
#' @export
#'
scranNormalise <- function(object, quickCluster = FALSE) {
    if (!is.null(object@Log$NormalisationMethod)) {
        stop("This data is already normalised.")
    }
    
    # Convert a scryeR object to a SCRAN object
    print("Converting EMSet to SCESet...")
    sce.obj <- ConvertToScater(object)
    
    # Remove controls from SCESet object
    present.mt <- which(rownames(sce.obj) %in% object@Controls$Mt)
    present.rb <- which(rownames(sce.obj) %in% object@Controls$Rb)

    remove.idx <- c(present.mt, present.rb)
    
    # Remove ERCC sequences if in dataset
    ercc.idx <- grep("^ercc-", ignore.case = TRUE, rownames(object@ExpressionMatrix))
    if (length(ercc.idx) > 0) {
        remove.idx <- c(remove.idx, ercc.idx)
    }
    
    # Remove these items
    sce.obj_rmMtRb <- sce.obj[-remove.idx, ]
    
    # Run computeSumFactors based on number of samples If we have more than 10,000 samples, we will run quickCluster Otherwise we supply a set list of sizes
    if (ncol(sce.obj_rmMtRb) > 10000) {
        if (quickCluster) {
            print(sprintf("%i cells detected. Running quickCluster to feed into computeSumFactors...", ncol(sce.obj)))
            quick.cluster <- scran::quickCluster(sce.obj_rmMtRb, method = "hclust")
        } else {
            print(sprintf("%i cells detected. Randomly grouping cells to feed into computeSumFactors...", ncol(sce.obj)))
            cell.identifiers <- colnames(sce.obj_rmMtRb)
            
            # Assign pseudocluster
            chunked.idx <- split(sample(1:length(cell.identifiers)), 1:10)
            quick.cluster <- cell.identifiers
            
            for (cluster.id in names(chunked.idx)) {
                cell.idx <- chunked.idx[[cluster.id]]
                quick.cluster[cell.idx] <- cluster.id
            }
        }
        
        # Feed rough clusters into computeSumFactors
        factored.sce.obj <- scran::computeSumFactors(sce.obj_rmMtRb, clusters = quick.cluster, positive = T)
    } else {
        print(sprintf("%i cells detected. Running computeSumFactors with preset sizes of 40, 60, 80, 100...", ncol(sce.obj)))
        preset.sizes <- c(40, 60, 80, 100)
        factored.sce.obj <- scran::computeSumFactors(sce.obj_rmMtRb, sizes = preset.sizes, positive = T)
    }
    
    print("scran's computeSumFactors complete. Removing zero sum factors from dataset...")
    zero.size.factors <- which(factored.sce.obj@phenoData@data$size_factor == 0)
    if (any(zero.size.factors)) {
        min.size.factor <- min(factored.sce.obj@phenoData@data$size_factor[-zero.size.factors])
        factored.sce.obj@phenoData@data$size_factor[zero.size.factors] <- min.size.factor
    }
    
    print("Running scater's normalize method...")
    dcvl.sce.obj <- scater::normalize(factored.sce.obj)
    
    # Convert log-transformed results back into counts
    print("Normalisation complete. Converting SCESet back to EMSet...")
    dcvl.matrix <- as.matrix(scater::norm_exprs(dcvl.sce.obj))
    unlog.dcvl.matrix <- UnLog2Matrix(dcvl.matrix)
    
    # Replace the object
    normalised.obj <- ReplaceExpressionMatrix(object, unlog.dcvl.matrix)
    normalised.obj <- SyncSlots(normalised.obj)
    normalised.obj@Log$NormalisationMethod <- "scranNormalise"
    normalised.obj@Log$Controls <- FALSE
    normalised.obj@Log$ExcludeControls <- list(Mt = TRUE, Rb = TRUE)
    
    remove(sce.obj)
    return(normalised.obj)
}

#' NormWithinBatch
#'
#' Called by NormaliseBatches. Performs normalisation within a batch.
#'
NormWithinBatch <- function(batch.id, expression.matrix = NULL, cell.info = NULL) {
    # Function called by NormaliseBatches
    barcodes <- cell.info[, 1][which(cell.info[, 2] == batch.id)]
    sub.mtx <- expression.matrix[, barcodes]
    
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
#' @param object An \linkS4class{EMSet} object.
#' @export
#'
NormaliseBatches <- function(object) {
    if (!is.null(object@Log$NormaliseBatches)) {
        stop("This data is already normalised.")
    }
    
    # Retrieve variables from EMSet object
    exprs.mtx <- GetExpressionMatrix(object, "matrix")
    cell.info <- GetCellInfo(object)
    unique.batch.identifiers <- unique(cell.info[, 2])
    
    # Loop to get batch-specific data PARALLEL
    print("Retrieving batch-specific data...")
    batch.data <- BiocParallel::bplapply(unique.batch.identifiers, NormWithinBatch, expression.matrix = exprs.mtx, cell.info = cell.info)
    
    # Unpacking results
    print("Scaling data...")
    norm.factors <- unlist(lapply(batch.data, function(x) x$NormalisationFactor))
    median.size <- median(norm.factors)
    sub.matrix.list <- lapply(batch.data, function(x) x$CpmMatrix)
    cpm.matrix <- data.frame(sub.matrix.list)
    scaled.matrix <- cpm.matrix * median.size
    
    # Load back into EMSet and write metrics
    print("Returning object...")
    colnames(scaled.matrix) <- cell.info[, 1]
    object@ExpressionMatrix <- ConvertMatrix(scaled.matrix, format = "sparse.matrix")
    object@Log$NormaliseBatches <- TRUE
    updated.object <- GenerateMetrics(object)
    return(updated.object)
}

#' NormaliseLibSize
#'
#' Normalise library sizes by scaling.
#' @param object A \linkS4class{EMSet} object. Please remove spike-ins from the expression matrix before normalising.
#' @export
#'
NormaliseLibSize <- function(object) {
    expression.matrix <- as.matrix(object@ExpressionMatrix)
    norm.factor <- colSums(expression.matrix)
    median.size <- median(norm.factor)
    
    unscaled.matrix <- Matrix::t(Matrix::t(expression.matrix)/norm.factor)
    normalised.exprs.mtx <- unscaled.matrix * median.size
    object@Log <- c(object@Log, list(NormaliseLibSize = TRUE))
    object@ExpressionMatrix <- ConvertMatrix(normalised.exprs.mtx, format = "sparse.matrix")
    new.object <- GenerateMetrics(object)
    new.object@Log <- c(object@Log, list(NormalisationMethod = "NormaliseLibSize"))
    return(new.object)
}


## SUBFUNCTIONS FOR NormaliseByRLE
#' CalculateNormFactor
#'
#' Calculate the normalisation factor between all cells.
#'
CalcNormFactor <- function(x, geo.means) {
    x.geo.means <- cbind(x, geo.means)
    x.geo.means <- x.geo.means[(x.geo.means[, 1] > 0), ]
    non.zero.median <- median(apply(x.geo.means, 1, function(y) {
        y <- as.vector(y)
        y[1]/y[2]
    }))
    return(non.zero.median)
}

#' CalcGeoMeans
#'
#' Calculate the geometric mean around a point.
#'
CalcGeoMeans <- function(x) {
    x <- x[x > 0]
    x <- exp(mean(log(x)))
    return(x)
}

#' NormaliseByRLE
#'
#' Normalisation of expression between cells, by scaling to relative log expression.
#' This method assumes all genes express a pseudo value higher than 0, and also assumes most genes
#'are not differentially expressed. Only counts that are greater than zero are considered in this normalisation method.
#' @param object A \linkS4class{EMSet} object that has undergone filtering. Please ensure spike-ins have been removed before using this function.
#' @export
#'
NormaliseByRLE <- function(object) {
    if (!is.null(object@Log$NormalisationMethod)) {
        stop("This data is already normalised.")
    }
    expression.matrix <- GetExpressionMatrix(object, format = "matrix")
    
    print("Calculating geometric means...")
    geo.means <- apply(expression.matrix, 1, CalcGeoMeans)
    
    print("Calculating normalisation factors...")
    norm.factor <- apply(expression.matrix, 2, function(x) CalcNormFactor(x, geo.means))
    
    print("Normalising data...")
    normalised.matrix <- Matrix::t(Matrix::t(expression.matrix)/norm.factor)
    object@ExpressionMatrix <- ConvertMatrix(normalised.matrix, "sparse.matrix")
    object@Log <- c(object@Log, list(NormalisationMethod = "NormaliseByRLE"))
    return(object)
}
