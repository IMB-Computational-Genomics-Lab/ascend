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
#' # Load example EMSet
#' em_set <- ascend::raw_set
#' 
#' # Normalise batches
#' norm_set <- normaliseBatches(em_set)
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom stats median
#' @importFrom S4Vectors merge
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
  
  # Extract components we require
  counts <- counts(object)
  col_info <- colInfo(object)
  col_data <- colData(object)
  
  # Merge into cell_info
  cell_info <- merge(col_info, col_data, by = "cell_barcode")
  
  # We want library sizes per batch
  batches <- sort(unique(cell_info$batch))
  info_df <- cell_info[, c("cell_barcode", "batch", "qc_libsize")]
  
  # Calculate between-batch sizeFactor
  batchSizeFactors <- sapply(batches, function(x) stats::median(info_df[which(info_df$batch == x), "qc_libsize"]))
  betweenBatchSizeFactor <- stats::median(batchSizeFactors)
  
  # For each batch, calculate batch sizeFactor
  scaleBatch <- function(x, info_df = NULL, counts = NULL, betweenBatchSizeFactor = NULL){
    median_libsize <- stats::median(info_df[which(info_df$batch == x), "qc_libsize"])
    batch_sizeFactor <- median_libsize / betweenBatchSizeFactor
    scaled_matrix <- (counts[, info_df[which(info_df$batch == x), "cell_barcode"]])/batch_sizeFactor
    return(list(sizeFactor = batch_sizeFactor, scaled_matrix = scaled_matrix))
  }
  
  # Get number of workers
  workers <- BiocParallel::bpnworkers(BiocParallel::bpparam())
  
  # Divide among workers if there are more
  # Otherwise no
  if (length(batches) >= workers){
    scaled_data <- BiocParallel::bplapply(batches, scaleBatch, 
                                          info_df = info_df, 
                                          counts = counts, 
                                          betweenBatchSizeFactor = betweenBatchSizeFactor)
  } else{
    scaled_data <- lapply(batches, scaleBatch, 
                          info_df = info_df, 
                          counts = counts, 
                          betweenBatchSizeFactor = betweenBatchSizeFactor)
  }
  
  names(scaled_data) <- batches
  
  # Get matrix
  scaled_matrix <- do.call(cbind, sapply(names(scaled_data), function(x) scaled_data[[x]][["scaled_matrix"]]))
  scaled_matrix <- scaled_matrix[rownames(counts), colnames(counts)]
  size_factors <- sapply(names(scaled_data), function(x) scaled_data[[x]][["sizeFactor"]])
  sizeFactors <- S4Vectors::DataFrame(cell_barcode = info_df$cell_barcode, sizeFactor = 0)
  
  for (batch in names(size_factors)){
    sizeFactors[which(info_df$batch == batch), "sizeFactor"] = size_factors[[batch]]
  }
  rownames(sizeFactors) <- sizeFactors$cell_barcode
  sizeFactors <- sizeFactors[colnames(counts), ]
  
  counts(object) <- scaled_matrix
  SingleCellExperiment::sizeFactors(object, "batch") <- sizeFactors$sizeFactor
  
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
#' # Load example EMSet
#' em_set <- ascend::raw_set
#' 
#' # Normalise batches
#' norm_set <- normaliseByRLE(em_set)
#' 
#' @importFrom BiocParallel bpvec
#' @importFrom Matrix t
#' @importFrom SingleCellExperiment normcounts logcounts sizeFactors
#' 
#' @export
#'
normaliseByRLE <- function(object) {
  print("Retrieving objects from EMSet")
  # Get data from object
  counts <- counts(object)
  rowInfo <- rowInfo(object)
  colInfo <- colInfo(object)
  rowData <- SummarizedExperiment::rowData(object)
  colData <- SummarizedExperiment::colData(object)
  gene_id_name <- colnames(rowData)[1]
  obj_types <- is(counts)
  
  # Drop nonzero counts
  counts <- counts[rowData[rowData$qc_ncells > 0, 1], ]
  
  # Calculate geometric means using Rcpp
  print("Calculating geometric means for nonzero genes...")
  if ("sparseMatrix" %in% obj_types){
    geo_means <- calcGeoMeansSparse(counts)
  } else{
    geo_means <- calcGeoMeansDense(counts)
  }
  
  # Calculate geoMeans
  print("Geometric means calculation complete!")
  geomean_list <- rownames(counts)
  geomean_df <- S4Vectors::DataFrame(gene_id = geomean_list, geoMean = geo_means)
  colnames(geomean_df) <- c(gene_id_name, "geoMean")
  rowData <- S4Vectors::merge(rowData, geomean_df, by = 1, all = TRUE)
  rowData$geoMean[which(is.na(rowData$geoMean))] <- 0
  rownames(rowData) <- rowData[, 1]
  
  # Retrieve original counts and rearrange rowData to match
  counts <- SingleCellExperiment::counts(object)
  rowData <- rowData[rownames(counts), ]
  
  # Use Rcpp to calculate normFactor
  geomeans_list <- rowData[, "geoMean"]
  print("Calculating size factors...")
  if ("sparseMatrix" %in% obj_types){
    sizeFactors <- calcNormFactorSparse(counts, geomeans_list)
  } else{
    sizeFactors <- calcNormFactorDense(counts, geomeans_list)
  }
  
  # Add sizeFactors to object
  print("Size factors calculated.")
  SingleCellExperiment::sizeFactors(object, "RLE") <- sizeFactors
  
  # Reorganise to save on memory because it is a memory hog for 22K cells
  # Garbage collection!
  # Apply sizeFactors
  print("Applying size factors to counts...")
  normcounts <- Matrix::t(Matrix::t(counts)/sizeFactors)
  normcounts <- normcounts[rownames(object), colnames(object)]
  normcounts <- as(normcounts, obj_types[1])
  SingleCellExperiment::normcounts(object) <- normcounts 
  print("Normalised counts stored in EMSet.")
  
  print("Converting to logcounts...")  
  remove(counts)
  logcounts <- log2(normcounts + 1)
  remove(normcounts)
  logcounts <- logcounts[rownames(object), colnames(object)]
  logcounts <- as(logcounts, obj_types[1])
  SingleCellExperiment::logcounts(object) <- logcounts
  # Get ready to store
  print("Log-normalised counts complete. Returning EMSet...")
  SummarizedExperiment::rowData(object) <- rowData[rownames(object), ]
  log <- progressLog(object)
  log$NormalisationMethod <- "RLE"
  progressLog(object) <- log
  return(object)
}
