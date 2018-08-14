################################################################################
#
# ascend_normalisation.R
# description: Functions related to the normalisation of data.
#
################################################################################

#' normWithinBatch
#'
#' Normalise cells within each batch.
#' 
#' @param object An \linkS4class{EMSet}.
#' @return An EMSet with values normalised within a single batch.
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
#' em_set <- newEMSet(assays = list(counts = count_matrix), colInfo = col_info)
#' 
#' # Normalise within batches
#' norm_set <- normWithinBatch(em_set)
#' 
#' @importFrom SingleCellExperiment counts sizeFactors
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors merge
#' @importFrom BiocParallel bplapply
#' @importFrom Matrix t rowSums
#'
#' @export
#' 
normWithinBatch <- function(object) {
  if (!is.null(progressLog(object)$normaliseWithinBatch)) {
    stop("This data is already normalised.")
  }
  # Retrieve variables from EMSet object
  expression_matrix <- as.matrix(SingleCellExperiment::counts(object), sparse = TRUE)
  col_info <- colInfo(object)
  col_data <- SummarizedExperiment::colData(object)
  cell_info <- S4Vectors::merge(col_info, col_data, by = "cell_barcode")
  batch_list <- unique(as.vector(cell_info[, "batch"]))
  
  func_normWithinBatch <- function(x, 
                                   expression_matrix = NULL,
                                   cell_info = NULL){
    
    # Function called by NormaliseBatches
    barcodes <- cell_info[which(cell_info[, "batch"] == x), "cell_barcode"]
    batch_matrix <- expression_matrix[, barcodes]
    
    # Collapse all cells into one cell
    collapsed_matrix <- Matrix::rowSums(batch_matrix)
    norm_factor <- sum(collapsed_matrix)
    
    ## Scale sub-matrix to normalisation factor
    norm_matrix <- Matrix::t(Matrix::t(batch_matrix)/norm_factor)
    
    # Output to list for further calculations
    norm_factor_list <- rep(norm_factor, length(barcodes))
    names(norm_factor_list) <- barcodes
    return(list(norm_factors = norm_factor_list, norm_matrix = as.matrix(norm_matrix)))
  }
  
  print("Normalising each batch...")
  norm_objs <- BiocParallel::bplapply(batch_list, func_normWithinBatch,
                                      expression_matrix = expression_matrix,
                                      cell_info = cell_info)
  names(norm_objs) <- batch_list
  norm_matrices <- sapply(batch_list, function(x) norm_objs[[x]]$norm_matrix)
  norm_matrix <- do.call("cbind", norm_matrices)
  norm_matrix <- norm_matrix[rownames(object), colnames(object)]
  norm_size_factors <- sapply(batch_list, function(x) norm_objs[[x]]$norm_factors)
  norm_size_factors <- unlist(norm_size_factors)[colnames(object)]
  
  SingleCellExperiment::counts(object) <- norm_matrix
  SingleCellExperiment::sizeFactors(object, "withinBatch") <- norm_size_factors
  
  print("Re-calculating QC metrics...")
  object <- calculateQC(object)
  progress_log <- progressLog(object)
  progress_log$normaliseWithinBatch <- TRUE
  progressLog(object) <- progress_log
  print("Batch normalisation complete! Returning object...")
  return(object)
}

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
#' em_set <- newEMSet(assays = list(counts = count_matrix), colInfo = col_info)
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
  
  # Retrieve variables from EMSet object
  expression_matrix <- Matrix::Matrix(SingleCellExperiment::counts(object), sparse = TRUE)
  col_info <- colInfo(object)
  col_data <- SummarizedExperiment::colData(object)
  cell_info <- S4Vectors::merge(col_info, col_data, by = "cell_barcode")
  
  # Generate list of batches
  batch_list <- as.numeric(sort(unique(cell_info[, "batch"])))
  calcNormFactor <- function(x, cell_info = NULL, expression_matrix = NULL){
    barcodes <- cell_info[which(cell_info$batch == x), "cell_barcode"]
    batch_matrix <- expression_matrix[, barcodes]
    y_sums <- Matrix::rowSums(batch_matrix)
    norm_factor <- stats::median(y_sums)
    return(norm_factor)
  }
  
  print("Calculating size factors...")
  batch_size_factors <- BiocParallel::bplapply(batch_list, calcNormFactor, 
                                               cell_info = cell_info, 
                                               expression_matrix = expression_matrix)
  
  # Calculate between-batch size factor
  between_batch_size_factor <- stats::median(as.vector(unlist(batch_size_factors)))
  
  # Scale counts
  print("Scaling counts...")
  scaled_matrix <- Matrix::t(Matrix::t(expression_matrix) / between_batch_size_factor)
 
  # Record size factor
  print("Storing data in EMSet...")
  batch_size_factor <- rep(between_batch_size_factor, ncol(expression_matrix))
  SingleCellExperiment::counts(object) <- as.matrix(scaled_matrix)
  SingleCellExperiment::sizeFactors(object, "batch") <- batch_size_factor 
  print("Re-calculating QC metrics...")
  progress_log <- progressLog(object)
  progress_log$normaliseBatches <- TRUE
  progressLog(object) <- progress_log
  object <- calculateQC(object)
  print("Batch normalisation complete! Returning object...")
  return(object)
}


#' normaliseByRLE
#'
#' Normalisation of expression between cells, by scaling to relative log 
#' expression (RLE). This method assumes all genes express a pseudo value higher 
#' than 0, and also assumes most genes are not differentially expressed. Only 
#' counts that are greater than zero are considered in this normalisation method.
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
#' em_set <- newEMSet(assays = list(counts = count_matrix), colInfo = col_info)
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
  if (!is.null(progressLog(object)$NormalisationMethod)){
    stop("This data is already normalised.")
  }
  
  calcNormFactor <- function(x, geo_means = NULL) {
    norm_factors <- apply(x, 2, function(y){
      x_geomeans <- cbind(y, geo_means)
      x_geomeans <- x_geomeans[(x_geomeans[, 1] > 0), ]
      nonzero_median <- median(apply(x_geomeans, 1, function(z) {
        z <- as.vector(z)
        z[1]/z[2]
      }))
      return(nonzero_median)    
    })
    return(norm_factors)
  }
  
  calcGeoMeans <- function(x){
    geo_means <- apply(x, 1, function(y){
      y <- y[y > 0]
      y <- exp(mean(log(y)))
      return(y)
    })
    return(geo_means)
  }
  
  # Calculate geometric means, then use to calculate normFactors
  expression_matrix <- as(SingleCellExperiment::counts(object), "sparseMatrix")
  
  print("Chunking matrix for calculations")
  chunked_genes <- split(rownames(expression_matrix), 1:BiocParallel::bpworkers())
  chunked_cells <- split(colnames(expression_matrix), 1:BiocParallel::bpworkers())
  chunked_genecounts <- BiocParallel::bplapply(chunked_genes, function(x) expression_matrix[x, ])
  chunked_cellcounts <- BiocParallel::bplapply(chunked_cells, function(x) expression_matrix[, x])
  
  print("Calculating geometric means...")
  # Calculate geometric means - memory intensive step
  geo_means <- BiocParallel::bplapply(chunked_genecounts, calcGeoMeans)
  geo_means <- unlist(geo_means, recursive = FALSE, use.names = FALSE)
  names(geo_means) <- as.vector(unlist(chunked_genes, use.names = FALSE))
  remove(chunked_genecounts)
  # Replace Nan with 0
  geo_means[is.nan(geo_means)] <- 0
  
  print("Calculating normalisation factor...")
  norm_factor <- BiocParallel::bplapply(chunked_cellcounts, calcNormFactor, geo_means = geo_means)
  norm_factor <- unlist(norm_factor, recursive = FALSE, use.names = FALSE)
  names(norm_factor) <- as.vector(unlist(chunked_cells, use.names = FALSE))
  remove(chunked_cellcounts)
  
  # Store size factors
  print("Scaling counts...")
  norm_matrix <- Matrix::t(Matrix::t(expression_matrix/norm_factor))
  remove(expression_matrix)
  log2_matrix <- log2(norm_matrix + 1)
  
  print("Storing normalised counts...")
  SingleCellExperiment::normcounts(object) <- norm_matrix
  SingleCellExperiment::logcounts(object) <- log2_matrix
  SingleCellExperiment::sizeFactors(object, "RLE") <- norm_factor
  SummarizedExperiment::rowData(object)$geo_mean <- geo_means
  
  # Update log and return object
  log <- progressLog(object)
  log$NormalisationMethod <- "RLE"
  progressLog(object) <- log
  return(object)
}
