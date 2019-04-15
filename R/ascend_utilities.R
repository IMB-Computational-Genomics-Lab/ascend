################################################################################
#
# ascend_utilities.R
# description: Functions that help other functions. 
#
################################################################################

convertToScone <- function(x){
  # Coercian to experiment loses rowData for some reason
  # Merge colInfo and colData
  metadata <- metadata(x)
  colInfo <- metadata$colInfo
  rowInfo <- metadata$rowInfo
  colData <- S4Vectors::merge(colInfo, SummarizedExperiment::colData(x), by = 1)
  rowData <- S4Vectors::merge(rowInfo, SummarizedExperiment::rowData(x), by = 1)
  colData$batch <- factor(colData$batch, levels = unique(colData$batch))
  
  # Convert to SummarizedExperiment
  x <- as(x, "SummarizedExperiment")
  SummarizedExperiment::colData(x) <- colData
  SummarizedExperiment::rowData(x) <- rowData
  rownames(x) <- rownames(rowData) 
  colnames(x) <- rownames(colData)
  
  qc_columns <- grep("^qc_", colnames(colData))
  batch <- grep("^batch$", colnames(colData))
  
  # Add to SconeExperiment
  out <- scone::SconeExperiment(x, which_qc = qc_columns, which_batch = batch)
  return(out)
}

convertToSCE <- function(x){
  # Retrieve EMSet-specific slots
  col_info <- colInfo(x)
  row_info <- rowInfo(x)
  log <- progressLog(x)
  cluster_analysis <- clusterAnalysis(x)
  
  # Convert into SingleCellExperiment
  object <- as(x, "SingleCellExperiment")
  
  # Load everything EMSet-related into metatadata 
  S4Vectors::metadata(object) <- list(colInfo = col_info,
                                      rowInfo = row_info,
                                      log = log,
                                      cluster_analysis)
  return(object)
}

convertSCEtoEMSet <- function(x){
  # Set ascend slots to retrieve
  ascend_slots <- c("colInfo", "rowInfo", "log", "clusterAnalysis")
  
  # Get metadata where these ascend slots are stored
  sce_metadata <- S4Vectors::metadata(x)
  
  # Retreive ascend metadata
  ascend_elements <- sce_metadata[names(sce_metadata) %in% ascend_slots]
  
  # Identify non-ascend entries
  non_ascend_indices <- which(!(names(sce_metadata) %in% ascend_slots))
  
  
  # If they are not ascend entries, retain for storage
  if (length(non_ascend_indices) > 0){
    metadata <- sce_metadata[non_ascend_indices]
  } else{
    metadata <- c()
  }
  
  # Cooerce into EMSet
  object <- as(x, "EMSet")
  
  # Restore old metadata to the object
  S4Vectors::metadata(object) <- metadata
  
  # Replace slots  
  colInfo(object) <- ascend_elements$colInfo[BiocGenerics::colnames(object), ]
  rowInfo(object) <- ascend_elements$rowInfo[BiocGenerics::rownames(object), ]
  progressLog(object) <- ascend_elements$log
  
  # If clustering has been done, replace the analysis
  if ("clusterAnalysis" %in% names(ascend_elements)){
    clusterAnalysis(object) <- ascend_elements$clusterAnalysis  
  }
  
  # Run QC
  object <- calculateQC(object)
  
  # Return to user
  return(object)
}


convertEMSetToSeurat <- function(x){
  if (is(x, "EMSet")){
    x <- convertToSCE(x)   
  } else{
    stop("Supplied object is not an EMSet.")
  }
  
  # Check if logcounts is present - required by Seurat
  if (!("logcounts" %in% SummarizedExperiment::assayNames(x))){
    if ("normcounts" %in% SummarizedExperiment::assayNames(x)){
      normcounts <- SingleCellExperiment::normcounts(x)
      logcounts <- log2(normcounts + 1)
    } else if ("counts" %in% SummarizedExperiment::assayNames(x)){
      counts <- SingleCellExperiment::counts(x)
      logcounts <- log2(counts + 1)
    } else{
      stop("Please arrange your assays into lists.")
    }
    SingleCellExperiment::logcounts(x) <- logcounts
  }
  
  # Convert to Seurat
  object <- Seurat::Convert(from = x, to = "seurat")
  return(object)
}

convertSeuratToEMSet <- function(x){
  x <- Seurat::Convert(x, to = "sce")
  colData <- SummarizedExperiment::colData(x)
  colData <- colData[, !(colnames(colData) %in% grep("^qc_", colnames(colData), value = TRUE))]
  SummarizedExperiment::colData(x) <- colData
  object <- EMSet(x)
  return(object)  
}

#' convert
#' 
#' Conversion function for working with the `ascend`` package. This conversion
#' function will convert an EMSet to the following classes: SingleCellExperiment,
#' Seurat, and SconeExperiment. This function will also convert 
#' SingleCellExperiment and Seurat objects to an EMSet. 
#' To convert between Seurat and SingleCellExperiment, it is recommended you use 
#' the instructions for Seurat [here](https://satijalab.org/seurat/conversion_vignette.html#converting-tofrom-singlecellexperiment).
#' 
#' @param x Object to convert
#' @param to sce, seurat, scone, EMSet
#' @return Object in specified format
#' 
#' @export
#' 
convert <- function(x, to = c("sce", "seurat", "scone", "EMSet")){
  # Check if EMSet first - as it inherits from SingleCellExperiment
  # Check if SingleCellExperiment
  # Check if Seurat
  # Otherwise dump
  if (is(x, "EMSet")){
    if (to == "sce"){
      object <- convertToSCE(x)  
    }
    if (to == "seurat"){
      object <- convertEMSetToSeurat(x)
    }
    if (to == "scone"){
      x <- convert(x, to = "sce")
      object <- convertToScone(x)
    }
  } else if(is(x, "SingleCellExperiment")){
    if (to == "EMSet"){
      object <- convertSCEtoEMSet(x)
    }
  } else if(is(x, "seurat")){
    if (to == "EMSet"){
      object <- convertSeuratToEMSet(x)
    }
  }else{
    stop("Supplied object is not recognised.")
  }
  return(object)
}




unLog2Matrix <- function(x){
  # Convert to matrix
  x <- as(x, "matrix")
  
  # Unlog the matrix
  unlogged_matrix <- 2^x
  
  # Subtract pseudocount of 1
  unlogged_matrix_sub_1 <- unlogged_matrix - 1

  # Make negative values 0
  unlogged_matrix_sub_1[unlogged_matrix_sub_1 < 0] <- 0
  
  # Make infinite values 0
  unlogged_matrix_sub_1[!is.finite(unlogged_matrix_sub_1)] <- 0
  
  return(unlogged_matrix_sub_1)
}

joinPaths <- function(x){
  if (length(x) > 1){
    x <- gsub("/$", "", x)
    path <- do.call("file.path", as.list(x))
    return(path)
  } else{
    return(x)
  }
}

fileCheck <- function(x) {
  if (!(file.exists(x))) {
    stop(sprintf("%s is missing", x))
  } else {
    return(FALSE)
  }
}

mergeDF <- function(x, y, z){
  # Check if column is in both data frames
  if (!(all((z %in% colnames(x)) & (z %in% colnames(y))))){
    stop("Specified column merge does not exist in all data frames.")
  }
  
  # We want to keep the first column so discard it from the combined_colnames
  replace_cols <- colnames(x) %in% colnames(y) & colnames(x) != z
  
  # Remove columns to be replaced from old_df
  if (length(which(!(replace_cols))) > 1){
    x <- x[, !(replace_cols)]    
  } else{
    merge_vector <- list()
    merge_vector[[z]] <- x[, 1]
    x <- S4Vectors::DataFrame(merge_vector, row.names = rownames(x))
  }

  # Merge
  merged_df <- S4Vectors::merge(x, y, by = z)
  return(merged_df)
} 