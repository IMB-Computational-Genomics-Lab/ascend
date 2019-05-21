################################################################################
#
# ascend_setters.R
# description: Methods related to the addition of data into an EMSet.
#
################################################################################

#' @rdname regcounts
#' @param value A matrix with confounding factors regressed out
#' @include ascend_objects.R
#' @include ascend_confoundingfactors.R
#' @export
setGeneric("regcounts<-", function(x, ..., value) standardGeneric("regcounts<-"))

#' @rdname regcounts
setReplaceMethod("regcounts", signature(x = "EMSet"), function(x, value){
  assays <- SummarizedExperiment::assays(x)
  assays$regcounts <- value
  
  # Add to assays
  SummarizedExperiment::assays(x) <- assays
  
  # Return
  return(x)
})


#' @rdname controls
#' @param value Named list of control genes
#' @include ascend_objects.R
#' @export
setGeneric("controls<-", function(x, ...,  value) standardGeneric("controls<-"))

#' @rdname controls
setReplaceMethod("controls", signature(x = "EMSet"), function(x, value){
  # Get row information
  row_info <- rowInfo(x)
  
  # Set control group to NULL by default. This is for non-control genes
  if (is.null(row_info$control_group)){
    row_info$control_group <- NA    
  }
  
  # For each control group...
  # Control group check in case user did not group the controls
  if (length(names(value)) >= 1){
    for (control_name in names(value)){
      gene_set <- value[[control_name]]
      row_info$control_group[which(row_info[ ,1] %in% gene_set)] <- control_name
    }
  } else{
    row_info$control_group[which(row_info[,1] %in% unlist(value))] <- "Control"
  }
  
  # Replace control information
  rownames(row_info) <- row_info[, 1]
  rowInfo(x) <- S4Vectors::DataFrame(row_info)
  
  # Get log information
  log <- progressLog(x)
  
  # Update log information with controls
  log$set_controls <- value
  log$controls <- TRUE
  progressLog(x) <- log
  
  x <- calculateQC(x)
  return(x)
})

#' @rdname clusterAnalysis
#' @param value List to store in the clusterAnalysis slot
#' @include ascend_objects.R
#' @export
setGeneric("clusterAnalysis<-", function(x, ..., value) standardGeneric("clusterAnalysis<-")) 

#' @rdname clusterAnalysis
setReplaceMethod("clusterAnalysis", "EMSet", function(x, value) {
  # Sync colInfo, matrix and subsequent objects
  x@clusterAnalysis <- value
  x
})

#' @rdname progressLog
#' @param value List to store in the progressLog slot
#' @include ascend_objects.R
#' @export
setGeneric("progressLog<-", function(x, ..., value) standardGeneric("progressLog<-")) 

#' @rdname progressLog
setReplaceMethod("progressLog", "EMSet", function(x, value) {
  # Sync colInfo, matrix and subsequent objects
  x@log <- value
  x
})

#' @rdname colInfo
#' @param value DataFrame to store in colInfo slot.
#' @include ascend_objects.R
#' @importClassesFrom S4Vectors DataFrame
#' @export
setGeneric("colInfo<-", function(x, ..., value) standardGeneric("colInfo<-")) 

#' @rdname colInfo
setReplaceMethod("colInfo", signature(x = "EMSet"), function(x, value) {
  # If it's a data frame
  if (is.data.frame(value)){
    value <- S4Vectors::DataFrame(value)
  }
  
  # Make rownames equal colInfo[ , 1]
  rownames(value) <- value[, 1]
  
  # Replace slot directory
  x@colInfo <- value
  x
})


#' @rdname rowInfo
#' @param value DataFrame to store in rowInfo slot.
#' @include ascend_objects.R
#' @importClassesFrom S4Vectors DataFrame
#' @export
setGeneric("rowInfo<-", function(x, ..., value) standardGeneric("rowInfo<-"))

#' @rdname rowInfo
setReplaceMethod("rowInfo", signature(x = "EMSet"),  function(x, value) {
  # Make rownames of value equal value[,1]
  if (is.data.frame(value)){
    value <- S4Vectors::DataFrame(value)
  }
  rownames(value) <- value[, 1]
  
  # Replace item in slot
  x@rowInfo <- value
  x
})

#' Replace entries of EMSet
#' 
#' @param x EMSet
#' @param i Row index
#' @param j Column index
#' @param value Replacement dataframe
#' @param ... ...
#' 
setReplaceMethod("[", c("EMSet", "ANY", "ANY", "EMSet"),
                 function(x, i, j, ..., value) {
                   # Extract infos from original dataset
                   col_info <- colInfo(x, withDimnames=FALSE)
                   row_info <- rowInfo(x, withDimnames=FALSE)
                   
                   if (!missing(i)) {
                     if (is.character(i)) {
                       fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
                       i <- .SummarizedExperiment.charbound(
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
