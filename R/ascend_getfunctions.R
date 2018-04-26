#' GetExpressionMatrix
#' 
#' Returns a data frame containing the expression matrix from an 
#' \linkS4class{EMSet}.
#'
#' @param object An \linkS4class{EMSet}
#' @param format Format of the returned matrix ("data.frame", "matrix")
#' @return Returns the expression matrix in the chosen format (data.frame or 
#' matrix)
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", 
#' "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve expression matrix as a data frame
#' data_frame <- GetExpressionMatrix(EMSet, "data.frame")
#' 
#' # Retrieve expression matrix as a matrix
#' matrix <- GetExpressionMatrix(EMSet, "matrix")
#' 
#' @include ascend_objects.R
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "GetExpressionMatrix", def = function(object, format) {
    standardGeneric("GetExpressionMatrix")
})

#' GetExpressionMatrix
#'
#' Returns a data frame containing the expression matrix from an 
#' \linkS4class{EMSet} object.
#'
#' @param object An \linkS4class{EMSet}
#' @param format Format of the returned matrix ("data.frame", "matrix")
#' @return Returns the expression matrix in the chosen format (data.frame or 
#' matrix)
#'
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", 
#' "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve expression matrix as a data frame
#' data_frame <- GetExpressionMatrix(EMSet, "data.frame")
#' 
#' # Retrieve expression matrix as a matrix
#' matrix <- GetExpressionMatrix(EMSet, "matrix")
#' 
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetExpressionMatrix", signature("EMSet"), 
          function(object, format = c("data.frame", "matrix", "sparseMatrix")) {
    if (missing(format)) {
        format <- "data.frame"
    }
    sparse.expression.matrix <- object@ExpressionMatrix
    output <- ConvertMatrix(sparse.expression.matrix, format = format)
    return(output)
})

#' GetControls
#' 
#' Generic for \code{\link{GetControls}}.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return This function returns a list of controls defined by the user.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleEMSet.rds"))
#' 
#' # Retrieve controls from the object
#' controls <- GetControls(EMSet)
#'
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "GetControls", def = function(object) {
    standardGeneric("GetControls")
})

#' GetControls
#'
#' Retrieve list of controls from \code{\linkS4class{EMSet}}.
#'
#' @param object An \code{\linkS4class{EMSet}} object.
#' @return This function returns a list of controls defined by the user.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleEMSet.rds"))
#' 
#' # Retrieve controls from the object
#' controls <- GetControls(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetControls", signature("EMSet"), function(object) {
    control.list <- object@Controls
    return(control.list)
})

#' GetCellInfo
#' 
#' Generic for \code{\link{GetCellInfo}}.
#' @importFrom methods setGeneric
#' @param object An \code{\linkS4class{EMSet}} object
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleEMSet.rds"))
#' 
#' # Retrieve data frame containing cell information from the object
#' cell.info <- GetCellInfo(EMSet)
#'
#' @export
setGeneric(name = "GetCellInfo", def = function(object) {
    standardGeneric("GetCellInfo")
})

#' GetCellInfo
#'
#' Retrieve cell information from an \code{\linkS4class{EMSet}}.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return A data frame with cell identifiers and associated information.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve data frame containing cell information from the object
#' cell.info <- GetCellInfo(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetCellInfo", signature("EMSet"), function(object) {
    return(object@CellInformation)
})

#' GetGeneInfo
#' 
#' Generic for \code{\link{GetGeneInfo}}.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return This function returns a data frame containing information about 
#' the genes in \code{object}.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve data frame containing gene information from the object
#' gene.info <- GetGeneInfo(EMSet)
#'
#' @importFrom methods setGeneric
setGeneric(name = "GetGeneInfo", def = function(object) {
    standardGeneric("GetGeneInfo")
})

#' GetGeneInfo
#'
#' Retrieve gene information from an \code{\linkS4class{EMSet}}.
#' @param object A \code{\linkS4class{EMSet}} object.
#' @return This function returns a data frame containing information about 
#' the genes in \code{object}.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve data frame containing gene information from the object
#' gene.info <- GetGeneInfo(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetGeneInfo", signature("EMSet"), function(object) {
    gene.annotation <- object@GeneInformation
    return(gene.annotation)
})

#' GetBatchMatrix
#' 
#' Generic for \code{\link{GetBatchMatrix}}.
#' @param object An \code{\linkS4class{EMSet}} object
#' @param batch.id Batch identifier that you would like to retrieve
#' @return Gene/transcript counts for cells in \code{batch.id}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve data frame containing gene information from the object
#' matrix1 <- GetBatchMatrix(EMSet, 1)
#'
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "GetBatchMatrix", def = function(object, batch.id) {
    standardGeneric("GetBatchMatrix")
})

#' GetBatchMatrix
#'
#' Retrieve a portion of the matrix by batch label
#' @param object An \code{\linkS4class{EMSet}} object
#' @param batch.id Batch identifier that you would like to retrieve
#' @return Gene/transcript counts for cells in \code{batch.id}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Retrieve data frame containing gene information from the object
#' matrix1 <- GetBatchMatrix(EMSet, 1)
#'
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetBatchMatrix", signature("EMSet"), function(object, batch.id) {
    expression.matrix <- GetExpressionMatrix(object, "matrix")
    cell.info <- GetCellInfo(object)
    batch.list <- cell.info[, 2]
    barcodes <- cell.info[,1][which(batch.list == batch.id)]
    batch.matrix <- expression.matrix[, barcodes]
    return(batch.matrix)
})

#' GetRandMatrix
#' 
#' Generic for \code{\link{GetRandMatrix}}.
#' @param object An \code{\linkS4class{EMSet}} object that has undergone 
#' clustering
#' 
#' @return This function returns a data frame containing cluster stability 
#' information.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve controls from the object
#' randMatrix <- GetRandMatrix(EMSet)
#'
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "GetRandMatrix", def = function(object) {
    standardGeneric("GetRandMatrix")
})

#' GetRandMatrix
#'
#' Retrieve the rand matrix used to determine clusters from an 
#' \code{\linkS4class{EMSet}} object.
#' 
#' @param object An \code{\linkS4class{EMSet}} object that has undergone 
#' clustering
#' @return This function returns a data frame containing cluster stability 
#' information.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve rand matrix from the object
#' randMatrix <- GetRandMatrix(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetRandMatrix", signature("EMSet"), function(object) {
    rand.matrix <- object@Clusters$KeyStats
    return(rand.matrix)
})

#' GetPCA
#' 
#' Generic for \code{\link{GetPCA}}.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return PCA matrix stored in \code{object}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve PCA matrix from the object
#' PCA <- GetPCA(EMSet)
#'
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "GetPCA", def = function(object) {
    standardGeneric("GetPCA")
})

#' GetPCA
#'
#' Retrieve the PCA matrix from an \code{\linkS4class{EMSet}} object that has 
#' undergone PCA reduction.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return PCA matrix stored in \code{object}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve PCA matrix from the object
#' PCA <- GetPCA(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetPCA", signature = "EMSet", function(object) {
    if (is.null(object@PCA$PCA)) {
        stop("Please use the RunPCA function on this object before using this function.")
    } else {
        pca.matrix <- object@PCA$PCA
        return(pca.matrix)
    }
})

#' GetDistanceMatrix
#' 
#' Generic for \code{\link{GetDistanceMatrix}}.
#' @importFrom methods setGeneric
#' @param object An \code{\linkS4class{EMSet}} object
#' @return Distance matrix stored in \code{object}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve distance matrix from the object
#' distanceMatrix <- GetDistanceMatrix(EMSet)
#'
#' @export
setGeneric(name = "GetDistanceMatrix", def = function(object) {
    standardGeneric("GetDistanceMatrix")
})

#' GetDistanceMatrix
#'
#' Retrieve the distance matrix from an \code{\linkS4class{EMSet}} object that 
#' has undergone clustering.
#' 
#' @param object An \code{\linkS4class{EMSet}} object
#' @return Distance matrix stored in \code{object}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve distance matrix from the object
#' distanceMatrix <- GetDistanceMatrix(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetDistanceMatrix", signature = "EMSet", function(object) {
    if (is.null(object@Clusters$DistanceMatrix)) {
        stop("Please use the RunCORE function on this object before using this function.")
    } else {
        distance.matrix <- object@Clusters$DistanceMatrix
        return(distance.matrix)
    }
})

#' GetHclust
#' 
#' Generic for \code{\link{GetHclust}}.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return This function returns the hclust object stored in \code{object}
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve hclust from the object
#' hClust <- GetHclust(EMSet)
#'
#' @importFrom methods setGeneric
#' @export
setGeneric(name = "GetHclust", def = function(object) {
    standardGeneric("GetHclust")
})

#' GetHclust
#'
#' Retrieve the Hclust object from an \code{\linkS4class{EMSet}} object that has 
#' undergone clustering.
#' @param object An \code{\linkS4class{EMSet}} object
#' @return This function returns the hclust object stored in \code{object}
#'
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", 
#' "ExampleClusteredEMSet.rds"))
#' 
#' # Retrieve hclust from the object
#' hClust <- GetHclust(EMSet)
#'
#' @include ascend_objects.R
#' @importFrom methods setGeneric setMethod
#' @export
setMethod("GetHclust", signature = "EMSet", function(object) {
    if (is.null(object@Clusters$Hclust)) {
        stop("Please use the RunCORE function on this object before using this function.")
    } else {
        hclust.obj <- object@Clusters$Hclust
        return(hclust.obj)
    }
})


