################################################################################
#
# ascend_setters.R
# description: Methods related to the addition of data into an EMSet.
#
################################################################################

#' @include ascend_objects.R
#' @export
setGeneric("clusterAnalysis<-", function(x, ..., value) standardGeneric("clusterAnalysis<-")) 

#' @include ascend_objects.R
#' @export
setReplaceMethod("clusterAnalysis", "EMSet", function(x, value) {
  # Sync colInfo, matrix and subsequent objects
  x@clusterAnalysis <- value
  x
})

#' @include ascend_objects.R
#' @export
setGeneric("progressLog<-", function(x, ..., value) standardGeneric("progressLog<-")) 

#' @include ascend_objects.R
#' @export
setReplaceMethod("progressLog", "EMSet", function(x, value) {
  # Sync colInfo, matrix and subsequent objects
  x@log <- value
  x
})

#' @include ascend_objects.R
#' @export
setGeneric("colInfo<-", function(x, ..., value) standardGeneric("colInfo<-")) 

#' @include ascend_objects.R
#' @export
setReplaceMethod("colInfo", "EMSet", function(x, value) {
  # Sync colInfo, matrix and subsequent objects
  x <- BiocGenerics:::replaceSlots(x, colInfo = value, check = FALSE)
  x <- x[ , value[,1]]
  x <- calculateQC(x)
  x
})

#' @include ascend_objects.R
#' @export
setGeneric("rowInfo<-", function(x, ..., value) standardGeneric("rowInfo<-"))

#' @include ascend_objects.R
#' @export
setReplaceMethod("rowInfo", "EMSet",  function(x, value) {
  x <- BiocGenerics:::replaceSlots(x, rowInfo = value, check = FALSE)
  x <- x[value[,1], ]
  x <- calculateQC(x)
  x
})

#' @include ascend_objects.R
#' @export
setReplaceMethod("[", c("EMSet", "ANY", "ANY", "EMSet"),
                 function(x, i, j, ..., value) {
                   # Extract infos from original dataset
                   col_info <- colInfo(x, withDimnames=FALSE)
                   row_info <- rowInfo(x, withDimnames=FALSE)
                   
                   if (!missing(i)) {
                     if (is.character(i)) {
                       fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
                       i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                         i, rownames(x), fmt
                       )
                     }
                     i <- as.vector(i)
                     row_info[i] <- rowInfo(value, withDimnames=FALSE)
                   }
                   
                   if (!missing(j)) {
                     if (is.character(j)) {
                       fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
                       j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                         j, colnames(x), fmt
                       )
                     }
                     j <- as.vector(j)
                     col_info[j] <- colInfo(value, withDimnames=FALSE)
                   }
                   
                   out <- callNextMethod()
                   out <- BiocGenerics:::replaceSlots(out, rowInfo=row_info, 
                                                      colInfo=col_info,
                                                      check=FALSE)
                   out <- calculateQC(out)
                   out
                 })