################################################################################
#
# ascend_getters.R
# description: Methods related to the retrieval of data into an EMSet.
#
################################################################################

#' clusterAnalysis
#' 
#' Convenience method to access cluster analysis data from an 
#' \linkS4class{EMSet}.
#' 
#' @param x \linkS4class{EMSet} 
#' @param ... ...
#' 
#' @return List of cluster-related analyses
#' @include ascend_objects.R
#' @export
setGeneric("clusterAnalysis", function(x, ...) standardGeneric("clusterAnalysis"))

#' @rdname clusterAnalysis
setMethod("clusterAnalysis", "EMSet", function(x){
  output <- x@clusterAnalysis
  return(output)
})

#' controls
#' 
#' Access set controls.
#' 
#' @param x \linkS4class{EMSet}
#' @param ... ...
#' 
#' @include ascend_objects.R
#' @export
setGeneric("controls", function(x, ...) standardGeneric("controls"))

#' @rdname controls
setMethod("controls", "EMSet", function(x){
  controls <- x@log$set_controls
  return(controls)
})

#' progressLog
#' 
#' Access log stored in the \linkS4class{EMSet}.
#' 
#' @param x \linkS4class{EMSet}
#' @param ... ...
#' @return A list of information.
#' 
#' @include ascend_objects.R
#' @export
setGeneric("progressLog", function(x, ...) standardGeneric("progressLog"))

#' @rdname progressLog
setMethod("progressLog", signature(x = "EMSet"), function(x){
  output <- x@log
  return(output)
})

#' colInfo
#' 
#' Access cell information stored in colInfo.
#' 
#' @param x \linkS4class{EMSet}
#' @param ... ...
#' @param withDimnames Retrieve rownames of DataFrame Default: TRUE
#' 
#' @return colInfo DataFrame
#' 
#' @include ascend_objects.R
#' @importFrom BiocGenerics rownames
#' @export
setGeneric("colInfo", function(x, ...) standardGeneric("colInfo"))

#' @rdname colInfo
setMethod("colInfo", "EMSet", function(x, withDimnames = TRUE){
  output <- x@colInfo
  
  # If user has set withDimnames
  if (withDimnames){
    BiocGenerics::rownames(output) <- colnames(x)
  }
  return(output)
})
 
#' rowInfo
#' 
#' Access gene information stored in rowInfo.
#'
#' @param x \linkS4class{EMSet}
#' @param ... ...
#' @param withDimnames TRUE or FALSE
#' 
#' @return rowInfo DataFrame
#' 
#' @include ascend_objects.R
#' @importFrom BiocGenerics rownames
#' @export
setGeneric("rowInfo", function(x, ...) standardGeneric("rowInfo"))

#' @rdname rowInfo
setMethod("rowInfo", "EMSet", function(x, withDimnames = TRUE){
  output <- x@rowInfo
  if (withDimnames){
    BiocGenerics::rownames(output) <- rownames(x)
  }
  return(output)
})
 
#' Subsetter for EMSet
#' 
#' @param x An \linkS4class{EMSet}
#' @param i Row index
#' @param j Col index
#' @param drop Drop stuff
#' 
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
    row_info <- row_info[i, , drop = drop]
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
    col_info <- col_info[j, ,drop=drop]
  }
  
  out <- callNextMethod()
  out <- BiocGenerics:::replaceSlots(out, colInfo = col_info, rowInfo = row_info, check=FALSE)
  out <- calculateQC(out)
  out
})
