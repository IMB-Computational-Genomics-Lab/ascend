#' FilterByExpressedGenesPerCell
#'
#' Filtered cells by the proportion of expressed genes in the cell. This step is
#' usually done after the other filtering steps and prior to normalisation.
#' @param object AEMSet (ASCEND) object that has been filtered by
#' \code{\link{FilterByOutliers}} and \code{\link{FilterByCustomControl}}.
#' @param pct.value Percentage threshold as a whole number. Default: 1
#' @export
#'
FilterByExpressedGenesPerCell <- function(object, pct.value = 1){
  filtered.object <- object
  expression.matrix <- filtered.object@ExpressionMatrix

  # Generate list of genes and cells to keep
  genes.per.cell <- Matrix::colSums(expression.matrix != 0)
  remove.genes <- genes.per.cell < ncol(expression.matrix) * (pct.value/100)

  if (any(remove.genes)){
    remove.gene.list <- names(which(remove.genes))
    remove.idx <- as.vector(unlist(which(remove.genes)))
    filtered.matrix <- expression.matrix[-c(remove.idx),]
  } else{
    remove.barcodes <- list()
    remove.idx <- list()
    filtered.matrix <- expression.matrix
  }

  # Updating the matrix
  filtered.object <- ReplaceExpressionMatrix(filtered.object, filtered.matrix)
  filtered.object <- SyncSlots(filtered.object)

  # Updating the log
  filtered.object@Log$FilterByExpressedGenesPerCell <- remove.barcodes

  # Update the data frame
  if (is.null(filtered.object@Log$FilteringLog)){
    filtered.df <- data.frame()
  } else{
    filtered.df <- filtered.object@Log$FilteringLog
  }

  output.list <- list(FilterByExpressedGenesPerCell = length(remove.barcodes))
  filtered.df <- cbind(filtered.df, output.list)
  filtered.object@Log$FilteringLog <- filtered.df

  return(filtered.object)
}

#' FilterByCustomControl
#'
#' Filter cells in an expression matrix based on the expression levels of a
#' specific control.
#' This function should be used AFTER the cells have undergone general filtering
#' with the \code{\link{FilterByOutliers}} function.
#'
#' @param control.name Name of the control group, as used in the named list supplied to the AEMSet object
#' @param pct.threshold Percentage threshold to filter cells by, as a whole number. Default: 20
#' @param object A \linkS4class{AEMSet} object.
#' @export
#'

FilterByCustomControl <- function(control.name = NULL, pct.threshold = 20, object){
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
    custom.filtered.matrix <- expression.matrix[, -c(discard.idx)]
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

  if (is.null(filtered.object@Log$FilterByCustomControl)){
    log <- list()
  } else{
    log <- filtered.object@Log$FilterByCustomControl
  }

  log[[control.name]] <- discard.barcodes
  filtered.object@Log$FilterByCustomControl <- log

  # Update the dable
  log.entry.name <- paste0("CellsFilteredBy", control.name)
  column <- list()
  column[[log.entry.name]] <- length(discard.barcodes)

  # Get Existing Table
  if (is.null(filtered.object@Log$FilteringLog)){
    filtering.log <- data.frame()
  } else{
    filtering.log <- filtered.object@Log$FilteringLog
  }

  filtering.log <- cbind(filtering.log, column)
  filtered.object@Log$FilteringLog <- filtering.log

  return(filtered.object)
}

#' FilterByControl
#'
#' Called by \code{\link{FilterByOutliers}}. This function identifies cells to
#' remove based on expression levels of control genes.
#' @export
#'
FilterByControl <- function(control.group, total.counts, expression.matrix ){
  # Get transcript counts for the controls
  control.bool <- rownames(expression.matrix) %in% control.group
  control.transcript.counts <- expression.matrix[control.bool, ]
  control.transcript.total.counts <- colSums(control.transcript.counts)
  control.pt.matrix <- (control.transcript.total.counts/total.counts)*100
  return(control.pt.matrix)
}

#' FindOutliers
#'
#' Adapted from \code{\link[scater]{isOutlier}}. Determines outliers based on MAD value.
#' @param values List of values
#' @param nmads Mean Absolute Deviation value threshold. Default: 3
#' @param type Direction to find outliers in - both, lower, upper. Default: both
#' @export
#'
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
#' \item{Low overall gene expression}
#' \item{Low number of expressed genes}
#' \item{Expression of control genes beyond set threshold}
#' }
#'
#' This function then loads the filtered expression matrix into the AEMSet object.
#'
#' @param object An \linkS4class{AEMSet} object
#' @param cell.threshold  Mean Absolute Deviation (MAD) value to filter cells by library size. Default: 3
#' @param control.threshold  Mean Absolute Deviation (MAD) value to filter cells by proportion of control genes. Default: 3
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

  # Retrieve required objects from AEMSet
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
  remove.cell.barcodes <- c(drop.barcodes.libsize, drop.barcodes.feature, drop.barcodes.controls)
  remove.cells.bool <- colnames(expression.matrix) %in% unique(names(remove.cell.barcodes))
  filtered.expression.matrix <- expression.matrix[, !remove.cells.bool]
  filtered.object@ExpressionMatrix <- filtered.expression.matrix

  ### Loading filtering log
  filtering.log <- list(CellsFilteredByLibSize = names(drop.barcodes.libsize),
                        CellsFilteredByLowExpression = names(drop.barcodes.feature),
                        CellsFilteredByControls = drop.barcodes.controls.names)

  # To go into the dataframe
  if (is.null(filtered.object@Log$FilteringLog)){
    filtering.df <- data.frame(
                              CellsFilteredByLibSize = length(filtering.log$CellsFilteredByLibSize),
                              CellsFilteredByExpression = length(filtering.log$CellsFilteredByLowExpression),
                              CellsFilteredByControls = length(unlist(filtering.log$CellsFilteredByControls)))
  } else{
    filtering.df <- filtered.object@Log$FilteringLog
    libsize <- filtering.df$CellsFilteredByLibSize + length(filtering.log$BarcodesFilteredByLibsize)
    expression <- filtering.df$CellsFilteredByExpression + length(filtering.log$BarcodesFilteredByLowExpression)
    controls <- filtering.df$CellsFilteredByControls + length(unlist(filtering.log$BarcodesFilteredByControls))
    filtering.df$CellsFilteredByLibSize <- libsize
    filtering.df$CellsFilteredByExpression <- expression
    filtering.df$CellsFilteredByControls <- controls
  }

  # Add to object
  filtered.object@Log$FilterByOutliers <- filtering.log
  filtered.object@Log$FilteringLog <- filtering.df

  # Finished!
  return(filtered.object)
}
