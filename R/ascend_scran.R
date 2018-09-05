################################################################################
#
# ascend_scran.R
# description: Wrappers for scran
#
################################################################################

#' scranCellCycle
#'
#' Wrapper for \pkg{scran}'s cell cycle functions. Please ensure you are using
#' *ensembl_id* as your rownames in this dataset. Also ensure you are using
#' mitochondrial and ribosomal genes as controls.
#'
#' @param object An \code{\linkS4class{EMSet}} object.
#' @param training.set A training dataset containing pairs of marker genes.
#' @return \code{object} with cell cycle information loaded into colInfo.
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment counts
#' @importFrom BiocParallel bpparam
#' @export
#'
scranCellCycle <- function(object, training_set = NULL) {
  
  if(is.null(training_set)){
    stop("Please specify a cyclone training set.")
  }
  
  expression_matrix <- SingleCellExperiment::counts(object)
  row_info <- rowInfo(object)
  
  # Run cyclone with cell cycle assignments as a vector
  print("Running cyclone from scran...")
  cc_assignments <- scran::cyclone(expression_matrix, pairs = training_set, rownames(expression_matrix), BPPARAM = BiocParallel::bpparam())
  col_info <- as.data.frame(colInfo(object))
  col_data <- as.data.frame(SummarizedExperiment::colData(object))
  cc_assignments <- as.data.frame(cc_assignments)
  col_info$phase <- cc_assignments$phases
  col_data <- S4Vectors::DataFrame(cbind(col_data), cc_assignments[ ,2:ncol(cc_assignments)])
  colInfo(object) <- S4Vectors::DataFrame(col_info)
  SummarizedExperiment::colData(object) <- col_data
  return(object)
}

#' scranNormalise
#'
#' Normalise an \code{\linkS4class{EMSet}} with \pkg{scran}'s deconvolution 
#' method by Lun et al. 2016.
#'
#' @details Users may choose to run computeSumFactors using either preset 
#' group sizes (40,60, 80 and 100) or quickCluster. For datasets with over
#' 20,000 cells - it is recommended you use quickCluster.
#' 
#' @param object An \code{\linkS4class{EMSet}} that has not undergone 
#' normalisation.
#' @param quickCluster Use scran's quickCluster method (TRUE) or  use randomly-
#' assigned groups (FALSE, Default).
#' @param min.mean Threshold for average counts. This argument is for the 
#' \code{\link[scran]{computeSumFactors}} function from \pkg{scran}. The value 
#' of 1 is recommended for read count data, while the default value of 1e-5 is 
#' best for UMI data. This argument is only used for newer versions of 
#' \pkg{scran}.
#' @return An \code{\linkS4class{EMSet}} with an expression matrix with counts 
#' normalised by \code{\link[scater]{normalize}} function.
#' 
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment isSpike spikeNames sizeFactors counts normcounts logcounts
#' @export
#'
scranNormalise <- function(object, quickCluster = FALSE, min.mean = 1e-5){
    if (!is.null(progressLog(object)$NormalisationMethod)) {
      stop("This data is already normalised.")
    }
    # Remove controls
    # Check if mitochondria and ribosomes are in the set
    log <- progressLog(object)
    gene_id_name <- colnames(rowInfo(object))[1]
    
    # Remove controls prior to normalisation
    if (any(names(log$set_controls) %in% c("Mt", "Rb"))){
      if ("Mt" %in% names(log$set_controls)){
        object <- excludeControl(object, control = "Mt")    
      }
      if ("Rb" %in% names(log$set_controls)){
        object <- excludeControl(object, control = "Rb")    
      }
    }
    
    col_info <- colInfo(object)
    row_info <- rowInfo(object)
    col_data <- SummarizedExperiment::colData(object)
    row_data <- SummarizedExperiment::rowData(object)
    log <- progressLog(object)
    cluster_analysis <- clusterAnalysis(object)
    old_size_factors <- SingleCellExperiment::sizeFactors(object)
    counts <- SingleCellExperiment::counts(object)
    
    # Save column names for extracted data frames so they can be retrieved
    col_info_headers <- BiocGenerics::colnames(col_info)
    row_info_headers <- BiocGenerics::colnames(row_info)
    col_data_headers <- BiocGenerics::colnames(col_data)
    row_data_headers <- BiocGenerics::colnames(row_data)
    sce_obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    
    # Remove spike-ins if present
    if (!is.null(SingleCellExperiment::isSpike(object))){
      spike_ins <- SingleCellExperiment::spikeNames(sce_obj)
      exclude_list <- sapply(spike_ins, function(x) which(SingleCellExperiment::isSpike(sce_obj, x)))
      exclude_list <- unique(exclude_list)
      sce_obj <- sce_obj[-exclude_list, ]
    }
    
    # Non-quickcluster method
    ncells <- BiocGenerics::ncol(sce_obj)
    
    if (quickCluster){
      print(sprintf("%i cells detected. Running computeSumFactors with quickCluster...", ncells))
      quick_cluster <- scran::quickCluster(sce_obj, method = "hclust", min.mean = min.mean)
      sce_obj <- scran::computeSumFactors(sce_obj, clusters = quick_cluster, positive = TRUE, min.mean = min.mean)
    } else{
      print(sprintf("%i cells detected. Running computeSumFactors with preset sizes of 40, 60, 80, 100...", ncells))
      preset_sizes <- c(40, 60, 80, 100)
      sce_obj <- scran::computeSumFactors(sce_obj, sizes = preset_sizes, positive = TRUE, min.mean = min.mean)
    }
    
    print("scran's computeSumFactors complete. Adjusting zero sum factors...")
    size_factors <- SingleCellExperiment::sizeFactors(sce_obj)
    zero_size_factors <- which(size_factors == 0)
    
    # Adjust size factor to use minimum 
    if (length(zero_size_factors) > 0){
      min_size_factor <- min(size_factors[-zero_size_factors])
      size_factors[zero_size_factors] <- min_size_factor
      SingleCellExperiment::sizeFactors(sce_obj) <- size_factors
    }
    
    print("Running scater's normalize method...")
    norm_obj <- scater::normalize(sce_obj, exprs_values = "counts")
    print("Normalisation complete. Converting SingleCellExperiment back to EMSet...")
    counts <- SingleCellExperiment::counts(norm_obj)
    logcounts <- SingleCellExperiment::logcounts(norm_obj)
    normcounts <- unLog2Matrix(logcounts)
    size_factors <- SingleCellExperiment::sizeFactors(norm_obj)
    
    # Retrieve cells that are still in the expression matrix
    cell_list <- colnames(normcounts)
    gene_list <- rownames(normcounts)
    
    # Trim so all the data is there
    colInfo <- col_info[which(col_info$cell_barcode %in% cell_list), ]
    rowInfo <- row_info[which(row_info[,1] %in% gene_list), ]
    keep_colData <- which(!(col_data_headers %in% grep("^qc_", col_data_headers, value = TRUE)))
    keep_rowData <- which(!(row_data_headers %in% grep("^qc_", row_data_headers, value = TRUE)))
    
    if (length(keep_colData) > 1){
      colData <- col_data[, keep_colData]
      colData <- colData[match(cell_list, colData$cell_barcode), ]
    } else{
      colData <- S4Vectors::DataFrame(cell_barcode = cell_list, row.names = cell_list) 
    }
    
    if (length(keep_rowData) > 1){
      rowData <- row_data[, keep_rowData]
      rowData <- rowData[match(gene_list, rowData[, gene_id_name]), ]
    } else{
      rowData <- S4Vectors::DataFrame(gene_list, row.names = gene_list)
      colnames(rowData) <- gene_id_name
    }
    
    norm_emset <- newEMSet(assays = list(counts = counts,
                                         normcounts = normcounts,
                                         logcounts = logcounts),
                           colInfo = colInfo,
                           rowInfo = rowInfo,
                           colData = colData,
                           rowData = rowData
    )
    
    SingleCellExperiment::sizeFactors(norm_emset, "norm_deconv") <- size_factors
    
    # Update log
    log$NormalisationMethod <- "Deconvolution"
    progressLog(norm_emset) <- log
    norm_emset <- calculateQC(norm_emset)
    return(norm_emset)
}
