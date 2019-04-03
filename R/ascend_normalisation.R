################################################################################
#
# ascend_normalisation.R
# description: Functions related to the normalisation of data.
#
################################################################################

#' normaliseBatches
#' 
#' Normalise counts to remove batch effects. This normalisation method is for
#' experiments where data from batches of different samples are combined without
#' undergoing library equalisation.
#'
#' This step should be done prior to any filtering and normalisation between cells.
#'
#' @param object An \code{\linkS4class{EMSet}} with cells from more than one batch.
#' @return An \code{\linkS4class{EMSet}} with batch-normalised expression values.
#' 
#' @examples
#' # Generate example matrix
#' count_matrix <- matrix(sample(0:1, 900, replace=TRUE),30,30)
#' colnames(count_matrix) <- paste0("Cell-", 1:ncol(count_matrix))
#' rownames(count_matrix) <- paste0("Gene-", 1:nrow(count_matrix))
#' 
#' # Generate example colInfo
#' col_info <- data.frame(cell_barcode = colnames(count_matrix))
#' col_info$batch <- sample(1:4, nrow(col_info), replace = TRUE)
#' 
#' # Create test EMSet
#' em_set <- EMSet(list(counts = count_matrix), colInfo = col_info)
#' 
#' # Normalise by RLE
#' norm_set <- normaliseBatches(em_set)
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom stats median
#' @importFrom S4Vectors merge
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment counts sizeFactors
#' 
#' @export
#'
normaliseBatches <- function(object){
  if (!is.null(progressLog(object)$normaliseBatches)) {
    stop("This data is already normalised.")
  }
  # 1. Check if there are more than one batches
  if ("batch" %in% colnames(colInfo(object))){
    col_info <- colInfo(object)
    batches <- unique(col_info$batch)
    if (length(batches) == 1){
      stop("Only one sample detected. Batch normalisation is not required.")
    }
  } else{
    stop("Please define batches in a column called 'batch' in colInfo")
  }
  
  # Extract data from object
  counts <- SingleCellExperiment::counts(object)
  col_data <- SummarizedExperiment::colData(object)
  cell_info <- S4Vectors::merge(col_info, col_data, by = "cell_barcode")
  cell_info <- cell_info[match(colnames(counts), cell_info$cell_barcode), ]
  
  # Pick out important information
  info_df <- cell_info[, c("cell_barcode", "batch", "qc_libsize")]
  batchSizeFactors <- sapply(batches, function(x) stats::median(info_df[which(info_df$batch == x), "qc_libsize"]))
  betweenBatchSizeFactor <- stats::median(batchSizeFactors)
  
  print("Scaling counts...")
  
  chunked_matrix <- lapply(batches, function(x){
    cell_barcodes <- info_df[info_df$batch == x, "cell_barcode"]
    return(counts[, cell_barcodes])
  })
  
  names(chunked_matrix) <- batches
  
  scaleBatch <- function(x, size_factor = NULL){
    batch_counts <- stats::median(Matrix::colSums(x))
    batch_norm_factor <- batch_counts / size_factor
    scaled_matrix <- Matrix::t(Matrix::t(x) / batch_norm_factor)
    return(list(norm_factor = batch_norm_factor, scaled_matrix = scaled_matrix))
  }
  
  scaled_data <- BiocParallel::bplapply(chunked_matrix, scaleBatch, size_factor = betweenBatchSizeFactor)
  
  # Collect data
  scaled_count_list <- lapply(scaled_data, function(x) x$scaled_matrix)
  scaled_counts <- do.call(cbind, scaled_count_list)
  size_factors <- sapply(scaled_data, function(x) x$norm_factor)
  size_factor_vector <- lapply(colnames(scaled_counts), function(x) size_factors[[cell_info$batch[cell_info$cell_barcode == x]]])
  names(size_factor_vector) <- colnames(scaled_counts)
  
  # Load data back into matrix
  SingleCellExperiment::counts(object) <- scaled_counts[rownames(object), colnames(object)]
  SingleCellExperiment::sizeFactors(object, "batch") <- size_factor_vector 
  
  print("Re-calculating QC metrics...")
  progress_log <- progressLog(object)
  progress_log$normaliseBatches <- TRUE
  progressLog(object) <- progress_log
  object <- calculateQC(object)
  
  return(object)
}

#' normaliseByRLE
#'
#' Normalisation of expression between cells, by scaling to relative log 
#' expression (RLE). Only counts that are greater than zero are considered in 
#' this normalisation method.
#' 
#' @param object An \code{\linkS4class{EMSet}} set that has undergone filtering. 
#' Please ensure spike-ins have been removed before using this function.
#' @return An \code{\linkS4class{EMSet}} with normalised expression values.
#'
#' @examples
#' # Generate example matrix
#' count_matrix <- matrix(sample(0:1, 900, replace=TRUE),30,30)
#' colnames(count_matrix) <- paste0("Cell-", 1:ncol(count_matrix))
#' rownames(count_matrix) <- paste0("Gene-", 1:nrow(count_matrix))
#' 
#' # Generate example colInfo
#' col_info <- data.frame(cell_barcode = colnames(count_matrix))
#' col_info$batch <- sample(1:4, nrow(col_info), replace = TRUE)
#' 
#' # Create test EMSet
#' em_set <- EMSet(list(counts = count_matrix), colInfo = col_info)
#' 
#' # Normalise by RLE
#' norm_set <- normaliseByRLE(em_set)
#' 
#' @importFrom BiocParallel bpvec
#' @importFrom Matrix t
#' @importFrom SingleCellExperiment normcounts logcounts sizeFactors
#' 
#' @export
#'
normaliseByRLE <- function(object) {
  # Functions to use
  calcNormFactor <- function(x, rowData = NULL) {
    x_geomeans <- cbind(x, rowData[names(x), "geoMean"])
    x_geomeans <- subset(x_geomeans, x_geomeans[,1 ] > 0)
    nonzero_median <- median(x_geomeans[,1] / x_geomeans[, 2])
    return(nonzero_median)
  }
  
  # Functions to use
  calcGeoMeans <- function(x){
    x <- x[x > 0]
    x<- exp(mean(log(x)))
    return(x)
  }
  
  if (!is.null(progressLog(object)$NormalisationMethod)){
    stop("This data is already normalised.")
  }
  
  # Calculate geometric means, then use to calculate normFactors
  expression_matrix <- SingleCellExperiment::counts(object)
  
  # Drop nonzero rows
  expression_matrix <- expression_matrix[Matrix::rowSums(expression_matrix) > 0, ]
  
  print("Calculating geometric means...")
  geo_means <- apply(expression_matrix, 1, calcGeoMeans)
  
  rowData <- SummarizedExperiment::rowData(object)
  gene_ids <- rowData[, 1]
  geomean_list <- list()
  geomean_list[[colnames(rowData)[1]]] <- gene_ids
  geomean_list$geoMean <- 0
  geomean_df <- S4Vectors::DataFrame(geomean_list)
  geomean_df$geoMean[which(geomean_df[,1] %in% names(geo_means))] <- as.vector(geo_means)
  rownames(geomean_df) <- geomean_df[, 1]
  rowData <- S4Vectors::merge(rowData, geomean_df, by = 1)
  rownames(rowData) <- rowData[, 1]
  
  expression_matrix <- SingleCellExperiment::counts(object)
  norm_factor <- apply(expression_matrix, 2, calcNormFactor, rowData = rowData)
  
  # Add to sizeFactors
  SingleCellExperiment::sizeFactors(object, "RLE") <- norm_factor
  
  # Scale counts
  print("Scaling counts...")
  norm_matrix <- Matrix::t(Matrix::t(expression_matrix/norm_factor))
  remove(expression_matrix)
  log2_matrix <- log2(norm_matrix + 1)
  
  print("Storing normalised counts")
  SingleCellExperiment::normcounts(object) <- norm_matrix[rownames(object), colnames(object)]
  SingleCellExperiment::logcounts(object) <- log2_matrix[rownames(object), colnames(object)]
  SummarizedExperiment::rowData(object) <- rowData[rownames(object), ]
  
  log <- progressLog(object)
  log$NormalisationMethod <- "RLE"
  progressLog(object) <- log
  return(object)
}
