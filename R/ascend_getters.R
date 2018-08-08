################################################################################
#
# ascend_getters.R
# description: Methods related to the retrieval of data into an EMSet.
#
################################################################################

#' @include ascend_objects.R
#' @export
setGeneric("clusterAnalysis", function(x, ...) standardGeneric("clusterAnalysis"))

#' clusterAnalysis
#' 
#' Convenience method to retrieve cluster analysis data from the slot.
#' 
#' @include ascend_objects.R
#' @export
setMethod("clusterAnalysis", "EMSet", function(x){
  output <- x@clusterAnalysis
  return(output)
})

#' @include ascend_objects.R
#' @export
setGeneric("progressLog", function(x, ...) standardGeneric("progressLog"))

#' progressLog
#' 
#' Covenience method to retrieve the progress log.
#' 
#' @include ascend_objects.R
#' @export
setMethod("progressLog", "EMSet", function(x){
  output <- x@log
  return(output)
})

#' @include ascend_objects.R
#' @export
setGeneric("colInfo", function(x, ...) standardGeneric("colInfo"))

#' colInfo
#' 
#' Convenience method to retrieve cell information.
#' 
#' @include ascend_objects.R
#' @importFrom BiocGenerics rownames
#' @export
setMethod("colInfo", "EMSet", function(x, withDimnames = TRUE){
  output <- x@colInfo
  if (withDimnames){
    BiocGenerics::rownames(output) <- colnames(x)
  }
  return(output)
})
 
#' @include ascend_objects.R
#' @export
setGeneric("rowInfo", function(x, ...) standardGeneric("rowInfo"))


#' rowInfo
#' 
#' Convenience method to retrieve gene information.
#' 
#' @include ascend_objects.R
#' @importFrom BiocGenerics rownames
#' @export
setMethod("rowInfo", "EMSet", function(x, withDimnames = TRUE){
  output <- x@rowInfo
  if (withDimnames){
    BiocGenerics::rownames(output) <- rownames(x)
  }
  return(output)
})

#' @include ascend_objects.R
#' @export
setMethod("[", "EMSet", function(x, i, j, drop=TRUE) {
  # Extract infos from original dataset
  col_info <- colInfo(x, withDimnames=FALSE)
  row_info <- rowInfo(x, withDimnames=FALSE)
  
  # If row is defined
  if (!missing(i)) {
    if (is.character(i)) {
      fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
      i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        i, rownames(x), fmt
      )
    }
    i <- as.vector(i)
    row_info <- row_info[i, , drop = FALSE]
  }
  
  # If column is defined
  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        j, colnames(x), fmt
      )
    }
    j <- as.vector(j)
    col_info <- col_info[j, ,drop=FALSE]
  }
  
  out <- callNextMethod()
  out <- BiocGenerics:::replaceSlots(out, colInfo = col_info, rowInfo = row_info, check=FALSE)
  out <- calculateQC(out)
  out
})
