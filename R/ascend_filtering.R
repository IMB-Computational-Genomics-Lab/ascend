#' FilterLowAbundanceGenes
#'
#' Removes genes if they are expressed in less than a set percentage of cells . 
#' This step is usually done after the other filtering steps and prior to 
#' normalisation. As this step may remove rare transcripts, this filtering step
#' is optional.
#' 
#' @param object An \code{\linkS4class{EMSet}} that has been filtered by 
#' \code{\link{FilterByOutliers}} and \code{\link{FilterByControl}}.
#' @param pct.value Percentage threshold as a whole number. Default: 1.
#' @return An \code{\linkS4class{EMSet}} with low abundance genes removed from 
#' the dataset.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Filter out low abundance genes
#' filtered_EMSet <- FilterLowAbundanceGenes(EMSet, pct.value = 1)
#' 
#' @importFrom Matrix rowSums
#' @export
#'
FilterLowAbundanceGenes <- function(object, pct.value = 1){
  filtered.object <- object
  expression.matrix <- filtered.object@ExpressionMatrix
  
  # Generate list of genes and cells to keep
  cells.per.gene <- Matrix::rowSums(expression.matrix)
  remove.genes <- cells.per.gene < (ncol(expression.matrix) * (pct.value/100))
  
  if (any(remove.genes)){
    removed.genes <- names(which(remove.genes))
    remove.idx <- which(rownames(expression.matrix) %in% removed.genes)
    filtered.matrix <- expression.matrix[-remove.idx,]
  } else{
    removed.genes <- list()
    filtered.matrix <- expression.matrix
  }
  
  # Updating the matrix
  filtered.object <- ReplaceExpressionMatrix(filtered.object, filtered.matrix)
  
  # Updating the log
  current.log <- filtered.object@Log
  
  # Update the data frame
  if (is.null(current.log$FilteringLog)){
    filtered.df <- data.frame(FilteredLowAbundanceGenes = length(removed.genes))
  } else{
    new.df <- data.frame(FilteredLowAbundanceGenes = length(removed.genes))
    filtered.df <- current.log$FilteringLog
    filtered.df <- cbind(filtered.df, new.df)
  }
  
  current.log$FilteringLog <- filtered.df
  current.log$FilteredLowAbundanceGenes <- removed.genes
  filtered.object@Log <- current.log 
  
  # Sync Object
  filtered.object <- SyncSlots(filtered.object)
  
  return(filtered.object)
}

#' FilterByControl
#'
#' Filter cells in an expression matrix based on the expression levels of a
#' specific control.
#' This function should be used AFTER the cells have undergone general filtering
#' with the \code{\link{FilterByOutliers}} function.
#'
#' @param control.name Name of the control group, as used in the named list 
#' supplied to the EMSet object.
#' @param pct.threshold Percentage threshold to filter cells by, as a whole 
#' number. Default: 20.
#' @param object An \code{\linkS4class{EMSet}}.
#' @return An \code{\linkS4class{EMSet}} with filtered controls.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Filter out cells where mitochondrial genes account for at least 20% of expression
#' filtered_EMSet <- FilterByControl(control.name = "Mt", pct.threshold = 20,
#' EMSet)
#' 
#' @export
#'

FilterByControl <- function(control.name = NULL, pct.threshold = 20, object){
  # Check in case user hasn't defined any controls.
  if(!object@Log$Controls){
    stop("Please define controls before attempting to filter this dataset.")
  }
  
  if(missing(control.name)){
    stop("Please specify a control name before using this function.")
  }
  
  # Get values related to this
  filtered.object <- object
  percentage.counts <- unlist(filtered.object@Metrics$PercentageTotalCounts[[control.name]])
  expression.matrix <- filtered.object@ExpressionMatrix
  
  # Perform test
  # Remove barcodes from the matrix
  discard.barcodes.test <- percentage.counts > pct.threshold
  
  if (any(discard.barcodes.test)){
    discard.barcodes <- names(which(discard.barcodes.test))
    discard.idx <- as.vector(which(discard.barcodes.test))
    custom.filtered.matrix <- expression.matrix[, -discard.idx]
  } else{
    discard.barcodes <- list()
    discard.idx <- c()
    custom.filtered.matrix <- expression.matrix
  }
  
  filtered.object <- ReplaceExpressionMatrix(filtered.object, custom.filtered.matrix)
  filtered.object <- SyncSlots(filtered.object)
  
  # Update the log
  log.entry <- list()
  log.entry[[control.name]] <- discard.barcodes
  
  if (is.null(filtered.object@Log$FilterByControl)){
    log <- list()
  } else{
    log <- filtered.object@Log$FilterByControl
  }
  
  if (length(log[[control.name]] > 0)){
    log[[control.name]] <- c(log[[control.name]], as.vector(unlist(discard.barcodes)))
  } else{
    log[[control.name]] <- as.vector(unlist(discard.barcodes))
  }
  
  filtered.object@Log$FilterByControl <- log
  
  # Update the dable
  log.entry.name <- paste0("CellsFilteredBy", control.name)
  column <- list()
  column[[log.entry.name]] <- length(log[[control.name]])
  
  # Get Existing Table
  if (is.null(filtered.object@Log$FilteringLog)){
    filtering.log <- data.frame(column)
  } else{
    filtering.log <- filtered.object@Log$FilteringLog
    filtering.log[[log.entry.name]] <- length(log[[control.name]])
  }
  filtered.object@Log$FilteringLog <- filtering.log
  return(filtered.object)
}

#' FilterControl
#'
#' Called by \code{\link{FilterByOutliers}}. This function identifies cells to
#' remove based on expression levels of control genes.
#' 
#' @param control.group Name of control group.
#' @param total.counts List of total counts of each cell.
#' @param expression.matrix An expression matrix.
#' @return Percentage expression of controls per cell.
#' @importFrom Matrix colSums
FilterControl <- function(control.group, total.counts, expression.matrix ){
  # Get transcript counts for the controls
  control.bool <- rownames(expression.matrix) %in% control.group
  control.transcript.counts <- expression.matrix[control.bool, ]
  control.transcript.total.counts <- colSums(control.transcript.counts)
  control.pt.matrix <- (control.transcript.total.counts/total.counts)*100
  return(control.pt.matrix)
}

#' FindOutliers
#'
#' Adapted from \code{\link[scater]{isOutlier}}. Determines outliers based on 
#' Mean Absolute Deviation (MAD) value.
#' @param values List of values.
#' @param nmads Mean Absolute Deviation value threshold. Default: 3.
#' @param type Direction to find outliers in - both (Default), lower, upper.
#' @param na.rm Remove NA values if present (Default: False).
#' @return A boolean of values that fall between the ranges set by \code{nmads} 
#' and \code{type}. 
#' @importFrom stats median mad
FindOutliers <- function(values, nmads = 3, type = c("both", "lower", "upper"), na.rm = FALSE) {
  med.val <- median(values, na.rm = na.rm)
  mad.val <- mad(values, center = med.val, na.rm = na.rm)
  upper.limit <- med.val + nmads * mad.val
  lower.limit <- med.val - nmads * mad.val
  if (type == "lower"){
    upper.limit <- Inf
  } else if (type == "higher") {
    lower.limit <- -Inf
  }
  return(values < lower.limit | upper.limit < values)
}

#' FilterByOutliers
#' Automatically filter cells based on expression levels
#'
#' These values are then used to filter out cells based on the following criteria:
#' \itemize{
#' \item{Low overall gene expression.}
#' \item{Low number of expressed genes.}
#' \item{Expression of control genes beyond set threshold.}
#' }
#'
#' This function then loads the filtered expression matrix into the EMSet object.
#'
#' @param object An \code{\linkS4class{EMSet}}.
#' @param cell.threshold  Mean Absolute Deviation (MAD) value to filter cells by 
#' library size. Default: 3.
#' @param control.threshold  Mean Absolute Deviation (MAD) value to filter cells 
#' by proportion of control genes. Default: 3.
#' @return An \code{\linkS4class{EMSet}} with outlier cells filtered out.
#' @examples
#' # Load EMSet
#' EMSet <- readRDS(system.file(package = "ascend", "extdata", "ExampleEMSet.rds"))
#' 
#' # Filter outliers for cells and controls by 3 MAD
#' filtered_EMSet <- FilterByOutliers(EMSet, cell.threshold = 3, 
#' control.threshold = 3)
#' 
#' @export
#'
FilterByOutliers <- function(object, cell.threshold = 3, control.threshold = 3) {
  filtered.object <- object
  
  # Input check
  if (!is.numeric(cell.threshold)){
    stop("Please set your Cell Threshold (NMAD value) to a valid integer.")
  }
  if(!is.numeric(control.threshold)){
    stop("Please set your Control Threshold (NMAD value) to a valid integer.")
  }
  # Stop if the user hasn't set any controls
  if(!object@Log$Controls){
    stop("Please define controls before filtering this dataset.")
  }
  
  # Retrieve required objects from EMSet
  expression.matrix <- filtered.object@ExpressionMatrix
  control.list <- filtered.object@Controls
  
  # Retrieve values from the object
  total.counts <- filtered.object@Metrics$TotalCounts
  log10.total.counts <- log10(total.counts)
  total.features.counts.per.cell <- filtered.object@Metrics$TotalFeatureCountsPerCell
  log10.total.features.counts.per.cell <- log10(total.features.counts.per.cell)
  percentage.lists.counts <- filtered.object@Metrics$PercentageTotalCounts
  
  ## Start identifying cells by barcodes
  cells.libsize <- FindOutliers(log10.total.counts, nmads=cell.threshold, type="lower") ## Remove cells with low expression
  cells.feature <- FindOutliers(log10.total.features.counts.per.cell, nmads=cell.threshold, type="lower") ## Remove cells with low number of genes
  
  ## Extract Indexes
  if (any(cells.libsize)){
    drop.barcodes.libsize <- which(cells.libsize)
  } else {
    drop.barcodes.libsize <- list()
  }
  
  if (any(cells.feature)){
    drop.barcodes.feature <- which(cells.feature)
  } else {
    drop.barcodes.feature <- list()
  }
  
  ## Identify cells to remove based on proportion of expression
  print("Identifying outliers...")
  controls.counts <- BiocParallel::bplapply(percentage.lists.counts, FindOutliers, nmads=control.threshold, type="higher") ## Use nmad to identify outliers
  print("Removing cells by library size...")
  drop.barcodes.controls <- BiocParallel::bplapply(controls.counts, which) ## Identify cell barcodes to remove
  print("Updating object information...")
  drop.barcodes.controls.names <- BiocParallel::bplapply(drop.barcodes.controls, names)
  
  ### Barcode master list of cells to remove
  remove.cell.barcodes <- c(names(drop.barcodes.libsize), names(drop.barcodes.feature))
  
  for (control in drop.barcodes.controls.names){
    remove.cell.barcodes <- c(remove.cell.barcodes, control) 
  }
  
  remove.cells.bool <- !(colnames(expression.matrix) %in% unique(remove.cell.barcodes))
  filtered.expression.matrix <- expression.matrix[, which(remove.cells.bool)]
  filtered.object <- ReplaceExpressionMatrix(filtered.object, filtered.expression.matrix)
  
  ### Loading filtering log
  if (is.null(filtered.object@Log$FilterByOutliers)){
    filtering.log <- list(CellsFilteredByLibSize = names(drop.barcodes.libsize),
                          CellsFilteredByLowExpression = names(drop.barcodes.feature),
                          CellsFilteredByControls = drop.barcodes.controls.names)    
  } else{
    filtering.log <- filtered.object@Log$FilterByOutliers
    filtering.log$CellsFilteredByLibSize <- c(filtering.log$CellsFilteredByLibSize, names(drop.barcodes.libsize))
    filtering.log$CellsFilteredByLowExpression <- c(filtering.log$CellsFilteredByLowExpression, names(drop.barcodes.feature))
    filtering.log$CellsFilteredByControls <- c(filtering.log$CellsFilteredByControls, names(drop.barcodes.controls.names))
  }
  
  
  # To go into the dataframe
  if (is.null(filtered.object@Log$FilteringLog)){
    filtering.df <- data.frame(
      CellsFilteredByLibSize = length(filtering.log$CellsFilteredByLibSize),
      CellsFilteredByExpression = length(filtering.log$CellsFilteredByLowExpression),
      CellsFilteredByControls = length(unlist(filtering.log$CellsFilteredByControls)))
  } else{
    filtering.df <- filtered.object@Log$FilteringLog
    filtering.df$CellsFilteredByLibSize <- length(filtering.log$CellsFilteredByLibSize)
    filtering.df$CellsFilteredByExpression <- length(filtering.log$CellsFilteredByLowExpression)
    filtering.df$CellsFilteredByControls <- length(unlist(filtering.log$CellsFilteredByControls))
  }
  
  # Add to object
  filtered.object@Log$FilterByOutliers <- filtering.log
  filtered.object@Log$FilteringLog <- filtering.df
  
  # Finished!
  return(filtered.object)
}
