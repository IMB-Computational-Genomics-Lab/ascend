#' ExcludeControl
#'
#' Removes the specified control from the expression matrix.
#' @param object A \linkS4class{AEMSet} object. It is recommended that you run this step after this object has undergone filtering.
#' @param control.name Name of the control set you want to remove from the dataset.
#'
ExcludeControl <- function(object, control.name){
  # Identify indices of control genes in the current expression matrix
  # Keep this matrix sparse for faster processing power
  expression.matrix <- object@ExpressionMatrix

  # Convert the control list into a boolean so we can remove rows from the sparse matrix
  control.list <- object@Controls[[ control.name ]]
  control.in.mtx <- rownames(object@ExpressionMatrix) %in% control.list

  # Remove control genes from the matrix by identified matrices
  endogenous.exprs.mtx <- expression.matrix[ !control.in.mtx, ]

  # Reload the expression matrix with the updated matrix
  object@ExpressionMatrix <- endogenous.exprs.mtx

  # Update the log
  updated.log <- list()
  updated.log[[control.name]] <- TRUE

  if (is.null(object@Log$Filtering[["RemoveControls"]])){
    remove.log <- list(RemoveControls = updated.log)
    object@Log$Filtering$RemoveControls <- remove.log
  } else {
    object@Log$Filtering$RemoveControls <- c(object@Log$Filtering$RemoveControls, updated.log)
  }

  # Regenerate metrics and return the updated object
  return.object <- GenerateMetrics(object)
  return(return.object)
}

#' FilterByExpressedGenesPerCell
#'
#' Filtered cells by the proportion of expressed genes in the cell. This step is usually done after the other filtering steps and prior to normalisation.
#' @param object AEMSet (scryeR) object that has been filtered by \code{\link{FilterByOutliers}} and \code{\link{FilterByCustomControl}}.
#' @param pct.value Percentage threshold as a decimal.
#'
FilterByExpressedGenesPerCell <- function(object, pct.value){
  expression.matrix <- object@ExpressionMatrix

  # Generate list of genes and cells to keep
  genes.per.cell <- Matrix::colSums(expression.matrix != 0)
  keep.genes <- genes.per.cell >= ncol(expression.matrix) * pct.value
  remove.genes <- genes.per.cell < ncol(expression.matrix) * pct.value

  keep.bool.arr <- colnames(expression.matrix) %in% names(keep.genes)
  filtered.matrix <- expression.matrix[keep.bool.arr,]

  # Updating the matrix
  object@ExpressionMatrix <- filtered.matrix

  # Updating the log
  removed.gene.list <- which(remove.genes)
  output.list <- list(FilterByExpressedGenesPerCell = length(removed.gene.list))
  updated.log <- list(FilterByExpressedGenesPerCell = removed.gene.list, FilteringList = output.list)
  object@Log$Filtering <- updated.log

  # Updating the metrics
  object <- UpdateBatchInfo(object)
  object <- GenerateMetrics(object)
  remove(expression.matrix)
  return(object)
}

#' FilterByCustomControl
#'
#' Filter cells in an expression matrix based on the expression levels of a specific control.
#' This function should be used AFTER the cells have undergone general filtering with the \code{\link{FilterByOutliers}} function.
#'
#' @param control.name Name of the control group, as used in the named list supplied to the AEMSet object
#' @param percentage.threshold Percentage threshold to filter cells by, as a decimal.
#' @param object A \linkS4class{AEMSet} object.
#'
FilterByCustomControl <- function(control.name, percentage.threshold, object){
  # Retrieve values required to run this function
  percentage.counts <- unlist(object@Metrics$PercentageTotalCounts[control.name])
  replace.names <- gsub(paste0(control.name, "."), "", names(percentage.counts))
  names(percentage.counts) <- replace.names

  # Perform test
  keep.barcodes.test <- percentage.counts <= percentage.threshold
  discard.barcodes.test <- percentage.counts > percentage.threshold

  if (any(keep.barcodes.test)){
    keep.barcodes <- names(which(keep.barcodes.test))
  } else {
    keep.barcodes <- list()
  }
  remove.barcodes <- names(which(keep.barcodes.test == FALSE))

  # Remove barcodes from the matrix
  expression.matrix <- object@ExpressionMatrix
  keep.barcodes.bool <- colnames(expression.matrix) %in% keep.barcodes
  custom.filtered.matrix <- expression.matrix[, keep.barcodes.bool]

  # Update the log
  log.entry.name <- paste0("BarcodesFilteredBy", control.name)
  updated.log <- list()
  updated.log[[log.entry.name]] <- remove.barcodes
  temp.log <- updated.log
  temp.log[[log.entry.name]] <- length(remove.barcodes)
  output.df <- data.frame(temp.log)
  updated.log <- c(updated.log, FilteringLog = output.df)

  # Update the object and return
  object@Log$Filtering <- updated.log
  object@ExpressionMatrix <- custom.filtered.matrix
  object <- UpdateBatchInfo(object)
  object <- GenerateMetrics(object)
  return(object)
}

#' FilterByControl
#'
#' Called by \code{\link{FilterByOutliers}}. This function identifies cells to remove based on expression levels of control genes.
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
#' @param nmads Mean Absolute Deviation value threshold
#' @param type Direction to find outliers in - both, lower, upper
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
#' @param object A \linkS4class{AEMSet} scryeR object
#' @param CellThreshold  Mean Absolute Deviation (MAD) value to filter cells by library size (Default Value: 3)
#' @param ControlThreshold  Mean Absolute Deviation (MAD) value to filter cells by proportion of control genes (Default Value: 3)
#'
FilterByOutliers <- function(object, CellThreshold = 3, ControlThreshold = 3) {
  filtered.object <- object

  # Input check
  if (!is.numeric(CellThreshold)){
    stop("Please set your Cell Threshold (NMAD value) to a valid integer.")
  }
  if(!is.numeric(ControlThreshold)){
    stop("Please set your Control Threshold (NMAD value) to a valid integer.")
  }

  # Retrieve required objects from AEMSet
  expression.matrix <- object@ExpressionMatrix
  control.list <- object@Controls

  # Retrieve values from the object
  total.counts <- object@Metrics$TotalCounts
  log10.total.counts <- log10(total.counts)
  total.features.counts.per.cell <- object@Metrics$TotalFeatureCountsPerCell
  log10.total.features.counts.per.cell <- log10(total.features.counts.per.cell)
  percentage.lists.counts <- object@Metrics$PercentageTotalCounts

  ## Start identifying cells by barcodes
  cells.libsize <- FindOutliers(log10.total.counts, nmads=CellThreshold, type="lower") ## Remove cells with low expression
  cells.feature <- FindOutliers(log10.total.features.counts.per.cell, nmads=CellThreshold, type="lower") ## Remove cells with low number of genes

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
  controls.counts <- BiocParallel::bplapply(percentage.lists.counts, function(x) FindOutliers(x, nmads=ControlThreshold, type="higher")) ## Use nmad to identify outliers
  print("Removing cells by library size...")
  drop.barcodes.controls <- BiocParallel::bplapply(controls.counts, function(x) which(x)) ## Identify cell barcodes to remove
  print("Updating object information...")
  drop.barcodes.controls.names <- BiocParallel::bplapply(drop.barcodes.controls, function(x) names(x))

    ### Barcode master list of cells to remove
  remove.cell.barcodes <- c(drop.barcodes.libsize, drop.barcodes.feature, drop.barcodes.controls)
  remove.cells.bool <- colnames(expression.matrix) %in% unique(names(remove.cell.barcodes))
  filtered.expression.matrix <- expression.matrix[,!remove.cells.bool]
  filtered.object@ExpressionMatrix <- filtered.expression.matrix

  ### Loading filtering log
  filtering.log <- list(BarcodesFilteredByLibSize = names(drop.barcodes.libsize), BarcodesFilteredByLowExpression = names(drop.barcodes.feature), BarcodesFilteredByControls = drop.barcodes.controls.names)

  # Prepare table to output
  writing.log.control <- lapply(filtering.log$BarcodesFilteredByControls, length)
  names(writing.log.control) <- lapply(names(writing.log.control), function(x) paste0("NumberOfBarcodesFilteredBy", x))

  writing.log <- list()
  # Loop to add metrics
  for ( x in names(filtering.log) ){
    if( is.null(filtering.log[[x]]) ){
      writing.log[[x]] <- 0
    } else {
      writing.log[[x]] <- length(filtering.log[[x]])
    }
  }
  writing.log.output <- c(writing.log, writing.log.control)
  writing.log.df <- as.data.frame(writing.log.output)
  filtered.object@Log <- list(Filtering=filtering.log, FilteringLog = writing.log.df)

  ### Rerun metrics
  filtered.object <- UpdateBatchInfo(filtered.object)
  filtered.object <- GenerateMetrics(filtered.object)
  return(filtered.object)
}
