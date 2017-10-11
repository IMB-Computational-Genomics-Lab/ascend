# Functions for retrieving functions from AEMSets
#' GetExpressionMatrix
#'
#' Returns a data frame containing the expression matrix from a \linkS4class{AEMSet} object.
#'
#' @param object A \linkS4class{AEMSet} to retrieve the expression matrix from
#' @param format Format of the returned matrix - 'data.frame' or 'matrix'
#' @return Returns the expression matrix in the chosen format (data.frame or matrix).
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetExpressionMatrix", def = function(object, format) {
    standardGeneric("GetExpressionMatrix")
})

setMethod("GetExpressionMatrix", signature("AEMSet"), function(object, format = c("data.frame", "matrix", "sparse.matrix")) {
    if (missing(format)) {
        format <- "data.frame"
    }
    sparse.expression.matrix <- object@ExpressionMatrix
    output <- ConvertMatrix(sparse.expression.matrix, format = format)
    return(output)
})

#' GetControls
#'
#' Retrieve list of controls from \linkS4class{AEMSet}.
#'
#' @param object An \linkS4class{AEMSet} object.
#' @return This function returns a list of controls defined by the user.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetControls", def = function(object) {
    standardGeneric("GetControls")
})

setMethod("GetControls", signature("AEMSet"), function(object) {
    control.list <- object@Controls
    return(control.list)
})

#' GetCellInfo
#'
#' Retrieve cell information from an \linkS4class{AEMSet}.
#' @return A data frame with cell identifiers and associated information.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetCellInfo", def = function(object) {
    standardGeneric("GetCellInfo")
})

setMethod("GetCellInfo", signature("AEMSet"), function(object) {
    return(object@CellInformation)
})

#' GetGeneInfo
#'
#' Retrieve gene information from an \linkS4class{AEMSet} object.
#' @param object A \linkS4class{AEMSet} object.
#' @return This function returns a data frame.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetGeneInfo", def = function(object) {
    standardGeneric("GetGeneInfo")
})

setMethod("GetGeneInfo", signature("AEMSet"), function(object) {
    gene.annotation <- object@GeneInformation
    return(gene.annotation)
})

#' GetBatchMatrix
#'
#' Retrieve a portion of the matrix by batch label
#' @param object \linkS4class{AEMSet} object
#' @param batch.id Batch identifier that you would like to retrieve
#' @export
#'
setGeneric(name = "GetBatchMatrix", def = function(object, batch.id) {
    standardGeneric("GetBatchMatrix")
})

setMethod("GetBatchMatrix", signature("AEMSet"), function(object, batch.id) {
    expression.matrix <- GetExpressionMatrix(object, "matrix")
    cell.info <- GetCellInfo(object)
    batch.list <- cell.info[, 2]
    barcodes <- names(batch.list[batch.list == batch.id])
    batch.matrix <- expression.matrix[, barcodes]
    return(batch.matrix)
})

#' GetRandMatrix
#'
#' Retrieve the rand matrix used to determine clusters from an \linkS4class{AEMSet} object.
#' @param object An \linkS4class{AEMSet} object that has undergone clustering.
#' @return This function returns a data frame.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetRandMatrix", def = function(object) {
    standardGeneric("GetRandMatrix")
})

setMethod("GetRandMatrix", signature("AEMSet"), function(object) {
    rand.matrix <- object@Clusters$KeyStats
    return(rand.matrix)
})

#' GetPCA
#'
#' Retrieve the PCA matrix from an \linkS4class{AEMSet} object that has undergone PCA reduction.
#' @param object An \linkS4class{AEMSet} object.
#' @return This function returns a data frame.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetPCA", def = function(object) {
    standardGeneric("GetPCA")
})

setMethod("GetPCA", signature = "AEMSet", function(object) {
    if (is.null(object@PCA$PCA)) {
        stop("Please use the RunPCA function on this object before using this function.")
    } else {
        pca.matrix <- object@PCA$PCA
        return(pca.matrix)
    }
})

#' GetDistanceMatrix
#'
#' Retrieve the distance matrix from an \linkS4class{AEMSet} object that has undergone clustering.
#' @param object An \linkS4class{AEMSet} object.
#' @return This function returns a distance matrix.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetDistanceMatrix", def = function(object) {
    standardGeneric("GetDistanceMatrix")
})

setMethod("GetDistanceMatrix", signature = "AEMSet", function(object) {
    if (is.null(object@Clusters$DistanceMatrix)) {
        stop("Please use the RunCORE function on this object before using this function.")
    } else {
        distance.matrix <- object@Clusters$DistanceMatrix
        return(distance.matrix)
    }
})

#' GetHclust
#'
#' Retrieve the Hclust object from an \linkS4class{AEMSet} object that has undergone clustering.
#' @param object An \linkS4class{AEMSet} object.
#' @return This function returns a hclust object.
#' @include ASCEND_objects.R
#' @export
setGeneric(name = "GetHclust", def = function(object) {
    standardGeneric("GetHclust")
})

setMethod("GetHclust", signature = "AEMSet", function(object) {
    if (is.null(object@Clusters$Hclust)) {
        stop("Please use the RunCORE function on this object before using this function.")
    } else {
        hclust.obj <- object@Clusters$Hclust
        return(hclust.obj)
    }
})


