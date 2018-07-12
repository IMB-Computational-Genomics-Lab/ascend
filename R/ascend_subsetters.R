################################################################################
#
# ascend_subsetters.R
# description: Methods related to the subsetting of EMSets.
#
################################################################################

#' @export
setGeneric(name = "excludeControl", def = function(object, ...., control) {
  standardGeneric("excludeControl")
})

#' excludeControl
#' 
#' Excludes a control group from an \linkS4class{EMSet}.
#' 
#' @param object An \linkS4class{EMSet}.
#' @param control Name of the control group you would like to exclude. It must
#' be in the `control` column of rowInfo.
#' 
#' @return An \linkS4class{EMSet} without the specified control group.
#' @export
setMethod("excludeControl", signature("EMSet"), function(object, 
                                                          control = c()){
  # Retrieve rowInfo
  row_info <- rowInfo(object)
  
  # Retrieve log
  log <- progressLog(object)
  control_list <- log$set_controls
  control_bool <- log$controls
  control_names <- names(control_list)
  
  if (!any(sapply(control, function(x) x %in% names(control_list)))){
    stop("Please make sure the selected control group has been defined in your dataset.")
  }
  
  # Check if controls are defined.
  if (!("control_group" %in% colnames(row_info))){
    stop("Please make sure controls are defined for your dataset.")
  }
  
  # Check if selected controls are present.
  if (!(all(control %in% row_info[, "control_group"]))){
    stop("Please make sure all selected controls are present in the control_group column of the rowInfo dataframe.")
  }
  
  # Identify rows to remove
  remove_rows <- which(row_info$control_group %in% control)
  
  control_list <- lapply(control_names, function(x) if(!x %in% control){control_list[[x]]})
  names(control_list) <- control_names
  control_list[sapply(control_list, is.null, USE.NAMES = TRUE)] <- NULL
  
  if (length(control_list) == 0){
    log$controls <- FALSE
  }
  log$set_controls <- control_list
  progressLog(object) <- log
  # Remove from EMSet
  return_set <- object[-remove_rows, ]
  
  return(return_set)
})

#' @export
setGeneric(name = "subsetCondition", def = function(object, ..., by, conditions) {
  standardGeneric("subsetCondition")
})

#' subsetCondition
#'
#' Subset cells of specific conditions from an \code{\linkS4class{EMSet}}.
#'
#' @param object An \code{\linkS4class{EMSet}}
#' @param by A list containing parameter to subset object by. These should be 
#' present as a column in colInfo and names in the named list `conditions`.
#' @param conditions  Conditions you would like to select. They should be 
#' organised into a named list, with the values specified in `by` as sublist
#' names.
#' 
#' @return An \code{\linkS4class{EMSet}} containing cells that have these 
#' conditions.
#' 
#' @examples
#' # Load a pre-existing EMSet
#' 
#' @include ascend_objects.R
#' @export
setMethod("subsetCondition", signature("EMSet"), function(object, 
                                                          by = c(), 
                                                          conditions = list()) {
  # Retrieve colInfo
  col_info <- colInfo(object, withDimnames = TRUE)
  
  # Check if by value is in colInfo
  if (! by %in% names(conditions)){
    stop("Please ensure your conditions are organised into a named list.")
  }
  if (!all(by %in% colnames(col_info))){
    stop("Please check if your defined condition(s) is/are present in colInfo.")
  }
  
  # Check if all conditions are present
  if (!all(sapply(by, function(x) all(conditions[[x]] %in% col_info[, x])))){
    stop("Please check if your defined condition(s) is/are present in colInfo.")
  }
  
  # Retrieve list of cells from col_info
  cell_list <- lapply(by, function(x) which(col_info[, x] %in% conditions[[x]]))
  cell_list <- unique(unlist(cell_list))
  
  # Subset EMSet
  subset_set <- object[ , cell_list]
  return(subset_set)
})