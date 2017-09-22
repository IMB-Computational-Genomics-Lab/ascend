# Functions for subsetting AEMSets
#' SubsetCondition
#'
#' Subset specific conditions from a \linkS4class{AEMSet} object.
#'
#' @param object A \linkS4class{AEMSet} object
#' @param conditions Name of the columns which contain boolean information on whether to keep or discard a cell.
#'
#' @return An \linkS4class{AEMSet} containing only this batch.
#' @include ASCEND_objects.R
#' @export
setGeneric(
  name = "SubsetCondition",
  def = function(object, conditions) {
    standardGeneric("SubsetCondition")
  }
)

setMethod("SubsetCondition", signature("AEMSet"), function(object, conditions = c()) {
  # Check if selected batches are present in the batches column
  if (!any(conditions %in% colnames(object@CellInformation))){
    stop("None of your selected batches are present in the dataset.")
  }

  # Create a new object to output, ensures original object does not get overwritten.
  subset.obj <- object

  # Retrieve data from relevant slots
  expression.matrix <- GetExpressionMatrix(subset.obj, format = "data.frame")
  cell.info <- subset.obj@CellInformation

  # Select out the specified columns
  subset.cell.info <- cell.info[ , which(colnames(cell.info) %in% conditions)]

  # Select out rows that are true in at least one column
  keep.list <- c()

  for (condition in conditions){
    keep.list <- c(keep.list, which(subset.cell.info[, condition]))
  }

  subset.cell.info <- subset.cell.info[which(subset.cell.info[,1] %in% keep.list),]

  # Subset the expression matrix.
  subset.matrix <- expression.matrix[,keep.list]
  subset.obj <- ReplaceExpressionMatrix(subset.obj, subset.matrix)
  subset.obj <- ReplaceCellInfo(subset.cell.info)
  subset.obj <- SyncSlots(subset.obj)

  # Clean up the object
  subset.obj@PCA <- list()
  subset.obj@Clusters <- list()
  subset.obj@Log <- c(subset.obj@Log, list(SubsetByBatches = TRUE, SubsettedBatches = batches))
  return(subset.obj)
})


#' SubsetBatch
#'
#' Subset a specific batch from a \linkS4class{AEMSet} object. This data is already normalised, but if you wish to recluster the data, you will need to use the RunPCA function again.
#'
#' @param object A \linkS4class{AEMSet} object
#' @param batches Name or number of the batch(es) you would like to subset.
#'
#' @return An \linkS4class{AEMSet} containing only this batch.
#' @include ASCEND_objects.R
#' @export
setGeneric(
  name = "SubsetBatch",
  def = function(object, batches) {
    standardGeneric("SubsetBatch")
  }
)

setMethod("SubsetBatch", signature("AEMSet"), function(object, batches = c()) {
  # Check if selected batches are present in the batches column
  if (!any(batches %in% object@CellInformation[,2])){
    stop("None of your selected batches are present in the dataset.")
  }

  # Create a new object to output, ensures original object does not get overwritten.
  subset.obj <- object

  # Retrieve data from relevant slots
  expression.matrix <- GetExpressionMatrix(subset.obj, format = "data.frame")
  cell.info <- subset.obj@CellInformation

  # Get barcodes for selected clusters
  subset.cell.info <- cell.info[which(cell.info[,2] %in% batches),]

  # Subset the expression matrix.
  subset.matrix <- expression.matrix[,unlist(subset.cell.info[,1])]
  subset.obj <- ReplaceExpressionMatrix(subset.obj, subset.matrix)
  subset.obj <- SyncSlots(subset.obj)

  # Clean up the object
  subset.obj@PCA <- list()
  subset.obj@Clusters <- list()
  subset.obj@Log <- c(subset.obj@Log, list(SubsetByBatches = TRUE, SubsettedBatches = batches))
  return(subset.obj)
})

#' SubsetCluster
#'
#' Subset a \linkS4class{AEMSet} object by cluster. Please make sure you have clustered with \code{\link{FindOptimalClusters}} before using this function.
#' This data is already normalised, but if you wish to recluster the data, you will need to use the RunPCA function again.
#'
#' @param object A \linkS4class{AEMSet} object
#' @param clusters Clusters to subset from the dataset.
#' @return Returns an \linkS4class{AEMSet} containing only this cluster.
#' @include ASCEND_objects.R
#' @export
#'
setGeneric(
  name = "SubsetCluster",
  def = function(object, clusters) {
    standardGeneric("SubsetCluster")
  }
)

setMethod("SubsetCluster", signature("AEMSet"), function(object, clusters = c()){
  # Check if data has been clustered
  if (is.null(object@CellInformation[,"cluster"])){
    stop("Please cluster your data before using this function.")
  }

  # Check if selected clusters are present in the clustered column
  if (!any(clusters %in% object@CellInformation[,"cluster"])){
    stop("None of your selected clusters are present in the dataset.")
  }

  # Create a new object to output, ensures original object does not get overwritten.
  subset.obj <- object

  # Retrieve data from relevant slots
  expression.matrix <- GetExpressionMatrix(subset.obj, format = "data.frame")
  cell.info <- subset.obj@CellInformation

  # Get barcodes for selected clusters
  subset.cell.info <- cell.info[which(cell.info$cluster %in% clusters),]

  # Subset the expression matrix.
  subset.matrix <- expression.matrix[,unlist(subset.cell.info[,1])]
  subset.obj <- ReplaceExpressionMatrix(subset.obj, subset.matrix)
  subset.obj <- SyncSlots(subset.obj)

  # Clean up the object
  subset.obj@PCA <- list()
  subset.obj@Clusters <- list()
  subset.obj@Log <- c(subset.obj@Log, list(SubsetByCluster = TRUE, SubsettedClusters = clusters))
  return(subset.obj)
})

#' SubsetCells
#'
#' Subset cells in the supplied list from an \linkS4class{AEMSet}.
#'
#' @param object An \linkS4class{AEMSet}
#' @param cell_barcodes A list of cell identifiers to subset from the \linkS4class{AEMSet}
#' @return An \linkS4class{AEMSet}
#' @include ASCEND_objects.R
#' @export
#'
setGeneric(
  name = "SubsetCells",
  def = function(object, cell_barcodes) {
    standardGeneric("SubsetCells")
  }
)

setMethod("SubsetCells", signature("AEMSet"), function(object, cell_barcodes = c()){
  # Check cells are in the expression matrix.
  expression.matrix <- GetExpressionMatrix(object, "data.frame")
  cell.info <- GetCellInfo(object)
  gene.info <- GetGeneInfo(object)
  control.list <- GetControls(object)
  present.cells <- cell_barcodes[cell_barcodes %in% colnames(expression.matrix)]

  # Missing cell catch
  if (length(present.cells) == 0){
    stop("All listed cells were not present in the AEMSet.")
  } else{
    # Subset out information
    subset.matrix <- expression.matrix[,present.cells]
    subset.obj <- ReplaceExpressionMatrix(object, subset.matrix)
    subset.obj <- SyncSlots(subset.obj)

    # Clear Clustering and PCA
    subset.obj@Clusters <- list()
    subset.obj@PCA <- list()

    return(subset.object)
  }
})
