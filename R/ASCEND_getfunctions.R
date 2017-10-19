# Functions for retrieving functions from EMSets
#' GetExpressionMatrix
#'
#' Returns a data frame containing the expression matrix from a \linkS4class{EMSet} object.
#'
#' @param object A \linkS4class{EMSet} to retrieve the expression matrix from
#' @param format Format of the returned matrix - 'data.frame' or 'matrix'
#' @return Returns the expression matrix in the chosen format (data.frame or matrix).
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetExpressionMatrix", def = function(object, format) {
    standardGeneric("GetExpressionMatrix")
})

setMethod("GetExpressionMatrix", signature("EMSet"), function(object, format = c("data.frame", "matrix", "sparse.matrix")) {
    if (missing(format)) {
        format <- "data.frame"
    }
    sparse.expression.matrix <- object@ExpressionMatrix
    output <- ConvertMatrix(sparse.expression.matrix, format = format)
    return(output)
})

#' GetControls
#'
#' Retrieve list of controls from \linkS4class{EMSet}.
#'
#' @param object An \linkS4class{EMSet} object.
#' @return This function returns a list of controls defined by the user.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetControls", def = function(object) {
    standardGeneric("GetControls")
})

setMethod("GetControls", signature("EMSet"), function(object) {
    control.list <- object@Controls
    return(control.list)
})

#' GetCellInfo
#'
#' Retrieve cell information from an \linkS4class{EMSet}.
#' @return A data frame with cell identifiers and associated information.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetCellInfo", def = function(object) {
    standardGeneric("GetCellInfo")
})

setMethod("GetCellInfo", signature("EMSet"), function(object) {
    return(object@CellInformation)
})

#' GetGeneInfo
#'
#' Retrieve gene information from an \linkS4class{EMSet} object.
#' @param object A \linkS4class{EMSet} object.
#' @return This function returns a data frame.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetGeneInfo", def = function(object) {
    standardGeneric("GetGeneInfo")
})

setMethod("GetGeneInfo", signature("EMSet"), function(object) {
    gene.annotation <- object@GeneInformation
    return(gene.annotation)
})

#' GetBatchMatrix
#'
#' Retrieve a portion of the matrix by batch label
#' @param object \linkS4class{EMSet} object
#' @param batch.id Batch identifier that you would like to retrieve
#' @export
#'
setGeneric(name = "GetBatchMatrix", def = function(object, batch.id) {
    standardGeneric("GetBatchMatrix")
})

setMethod("GetBatchMatrix", signature("EMSet"), function(object, batch.id) {
    expression.matrix <- GetExpressionMatrix(object, "matrix")
    cell.info <- GetCellInfo(object)
    batch.list <- cell.info[, 2]
    barcodes <- names(batch.list[batch.list == batch.id])
    batch.matrix <- expression.matrix[, barcodes]
    return(batch.matrix)
})

#' GetRandMatrix
#'
#' Retrieve the rand matrix used to determine clusters from an \linkS4class{EMSet} object.
#' @param object An \linkS4class{EMSet} object that has undergone clustering.
#' @return This function returns a data frame.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetRandMatrix", def = function(object) {
    standardGeneric("GetRandMatrix")
})

setMethod("GetRandMatrix", signature("EMSet"), function(object) {
    rand.matrix <- object@Clusters$KeyStats
    return(rand.matrix)
})

#' GetPCA
#'
#' Retrieve the PCA matrix from an \linkS4class{EMSet} object that has undergone PCA reduction.
#' @param object An \linkS4class{EMSet} object.
#' @return This function returns a data frame.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetPCA", def = function(object) {
    standardGeneric("GetPCA")
})

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
#' Retrieve the distance matrix from an \linkS4class{EMSet} object that has undergone clustering.
#' @param object An \linkS4class{EMSet} object.
#' @return This function returns a distance matrix.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetDistanceMatrix", def = function(object) {
    standardGeneric("GetDistanceMatrix")
})

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
#' Retrieve the Hclust object from an \linkS4class{EMSet} object that has undergone clustering.
#' @param object An \linkS4class{EMSet} object.
#' @return This function returns a hclust object.
#' @include ascend_objects.R
#' @export
setGeneric(name = "GetHclust", def = function(object) {
    standardGeneric("GetHclust")
})

setMethod("GetHclust", signature = "EMSet", function(object) {
    if (is.null(object@Clusters$Hclust)) {
        stop("Please use the RunCORE function on this object before using this function.")
    } else {
        hclust.obj <- object@Clusters$Hclust
        return(hclust.obj)
    }
})


