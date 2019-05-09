################################################################################
#
# ascend_filtering.R
# description: Functions related to the filtering of data
#
################################################################################

#' filterLowAbundanceGenes
#'
#' Removes genes if they are expressed in less than a set percentage of cells . 
#' This step is usually done after the other filtering steps and prior to 
#' normalisation. As this step may remove rare transcripts, this filtering step
#' is optional.
#' 
#' @param object An \linkS4class{EMSet} that has been filtered by 
#' \code{\link{filterByOutliers}} and \code{\link{filterByControl}}.
#' @param pct.threshold Percentage threshold as a whole number. Default: 1.
#' @return An \linkS4class{EMSet} with low abundance genes removed from the 
#' dataset.
#' 
#' @examples
#' # Load EMSet
#' raw_emset <- ascend::raw_set
#' 
#' # Filter low abundance genes expressed in less than 1% of cells
#' filtered_emset <- filterLowAbundanceGenes(raw_emset, pct.threshold = 0.1)
#' 
#' @importFrom SummarizedExperiment colData rowData 
#' @importFrom Matrix rowSums
#' @export
#'
filterLowAbundanceGenes <- function(object, pct.threshold = 1){
  # Merge data together to make it easier to work with
  cell_info <- S4Vectors::merge(colInfo(object), 
                                SummarizedExperiment::colData(object), 
                                by = "cell_barcode")
  gene_id_name <- colnames(rowInfo(object))[1]
  gene_info <- S4Vectors::merge(rowInfo(object),
                                SummarizedExperiment::rowData(object),
                                by = gene_id_name)
  # Reorder to match expression matrix
  cell_info <- cell_info[match(colnames(object), cell_info$cell_barcode), ]
  gene_info <- gene_info[match(rownames(object), gene_info[, gene_id_name]), ]
  
  # Generate list of genes and cells to keep
  cells_per_gene <- gene_info$qc_ncells
  threshold <- ncol(object) * (pct.threshold/100)
  remove_genes_indices <- which(cells_per_gene < threshold)
  
  if (length(remove_genes_indices) > 0){
    removed_genes <- gene_info[remove_genes_indices, 1]
    filtered_object <- object[which(!rownames(object) %in% removed_genes), ]
  } else{
    removed_genes <- list()
    filtered_object <- object
  }
  
  # Updating the log
  current_log <- progressLog(filtered_object)
  updated_log <- current_log
  
  if (!is.null(updated_log$filterLowAbundanceGenes)){
    updated_log$RemovedLowAbundanceGenes <- c(current_log$RemovedLowAbundanceGenes, removed_genes)
  } else{
    updated_log$RemovedLowAbundanceGenes <- removed_genes
  }
  
  # Update the data frame
  if (is.null(updated_log$FilteringLog)){
    filtered_df <- data.frame(FilteredLowAbundanceGenes = length(removed_genes))
  } else{
    filtered_df <- updated_log$FilteringLog
    filtered_df$FilteredLowAbundanceGenes <- length(updated_log$RemovedLowAbundanceGenes)
  }
  updated_log$FilteringLog <- filtered_df
  progressLog(filtered_object) <- updated_log
  # Validate controls
  
  if (!is.null(controls(filtered_object))){
    controls <- controls(filtered_object)
    updated_controls <- lapply(names(controls), function(x){
      if (any(!(controls[[x]] %in% rownames(filtered_object)))){
        return(controls[[x]][which(controls[[x]] %in% rownames(filtered_object))])
      } else{
        return(controls[[x]])
      }
    })
    names(updated_controls) <- names(controls)
    updated_controls[sapply(updated_controls, is.null)] <- NULL
    controls(filtered_object) <- updated_controls
  }
  
  return(filtered_object)
}


#' filterByControl
#'
#' Filter cells in an expression matrix based on the expression levels of a
#' specific control.
#' This function should be used AFTER the cells have undergone general filtering
#' with the \code{\link{filterByOutliers}} function.
#' 
#' @param object An \linkS4class{EMSet}.
#' @param control Name of the control group, as used in the named list 
#' supplied to the EMSet object.
#' @param pct.threshold Percentage threshold to filter cells by, as a whole 
#' number. Default: 20.
#' @return An \linkS4class{EMSet} with filtered controls.
#' 
#' @examples
#' # Load EMSet
#' raw_emset <- ascend::raw_set
#' 
#' # Filter by outliers with default settings
#' filtered_emset <- filterByControl(raw_emset, control = "Mt", 
#' pct.threshold = 20)
#'
#' @importFrom S4Vectors merge
#' @importFrom SummarizedExperiment colData rowData
#' 
#' @export
#'
filterByControl <- function(object, control = NULL, pct.threshold = 20){
    # Filter by control
    # Check in case user hasn't defined any controls.
    if(!progressLog(object)$controls){
      stop("Please define controls before attempting to filter this dataset.")
    }
    
    if(is.null(control)){
      stop("Please specify a control name before using this function.")
    } else{
      if (!control %in% names(progressLog(object)$set_controls)){
        stop("Please specify a control group present in this data.")
      }
    }
    
    filterControl <- function(control_group, pct.threshold = NULL, cell_info = NULL){
      # Get transcript counts for the controls
      control_transcript_counts <- cell_info[, paste0("qc_", control_group, "_pct_counts")]
      discard_indices <- which(control_transcript_counts > pct.threshold)
      discard_barcodes <- as.vector(cell_info$cell_barcode[discard_indices])
      return(discard_barcodes)
    }
    
    # Extract data from object
    # Merge data together to make it easier to work with
    gene_id_name <- colnames(rowInfo(object))[1]
    cell_info <- S4Vectors::merge(colInfo(object), 
                                  SummarizedExperiment::colData(object), 
                                  by = "cell_barcode")
    gene_info <- S4Vectors::merge(rowInfo(object),
                                  SummarizedExperiment::rowData(object),
                                  by = gene_id_name)
    
    # Reorder by matrix
    cell_info <- cell_info[match(colnames(object), cell_info[, "cell_barcode"]), ]
    gene_info <- gene_info[match(rownames(object), gene_info[, gene_id_name]), ]
    
    # Remove barcodes from the matrix
    discard_barcodes <- filterControl(control, pct.threshold = pct.threshold, cell_info = cell_info)
    
    if (length(discard_barcodes) > 0){
      filtered_object <- object[, which(!colnames(object) %in% discard_barcodes)]
    } else{
      discard_barcodes <- list()
      filtered_object <- object
    }
    
    # Update the log
    log <- progressLog(filtered_object)
    log_entry <- list()
    log_entry[[control]] <- discard_barcodes
    
    if (is.null(log$filterByControl)){
      control_log <- list()
    } else{
      control_log <- log$filterByControl
    }
    
    if (length(control_log[[control]] > 0)){
      control_log[[control]] <- c(control_log[[control]], log_entry[[control]])
    } else{
      control_log <- c(control_log, log_entry)
    }
    
    log$filterByControl <- control_log
    
    # Update the dable
    log_colname<- paste0("CellsFilteredBy", control, "Pct")
    column <- list()
    column[[log_colname]] <- length(control_log[[control]])
    
    # Get Existing Table
    if (is.null(log$FilteringLog)){
      filtering_log <- data.frame(column)
    } else{
      filtering_log <- log$FilteringLog
      filtering_log[[log_colname]] <- length(control_log[[control]])
    }
    
    log$FilteringLog <- filtering_log
    progressLog(filtered_object) <- log
    return(filtered_object)
}

#' filterByOutliers
#' 
#' Automatically filter cells based on expression levels.
#' These values are then used to filter out cells based on the following criteria:
#' \itemize{
#' \item{Low overall gene expression.}
#' \item{Low number of expressed genes.}
#' \item{Expression of control genes beyond set threshold.}
#' }
#'
#' This function then loads the filtered expression matrix into the \linkS4class{EMSet}.
#'
#' @param object An \linkS4class{EMSet}.
#' @param cell.threshold  Mean Absolute Deviation (MAD) value to filter cells by 
#' library size. Default: 3.
#' @param control.threshold  Mean Absolute Deviation (MAD) value to filter cells 
#' by proportion of control genes. Default: 3.
#' @param gene.threshold Mean Absolute Deviation (MAD) value to filter cells by
#' number of expressed genes. Default: 3
#' @return An \linkS4class{EMSet} with outlier cells filtered out.
#' 
#' @examples
#' # Load EMSet
#' raw_emset <- ascend::raw_set
#' 
#' # Filter by outliers with default settings
#' filtered_emset <- filterByOutliers(raw_emset, cell.threshold = 3, 
#' gene.threshold = 1.5, control.threshold = 3)
#' 
#' @importFrom dplyr semi_join
#' @importFrom BiocParallel bplapply
#' @importFrom S4Vectors merge
#' @importFrom SummarizedExperiment colData rowData
#' 
#' @export
#'
filterByOutliers <- function(object, 
                             cell.threshold = 3, 
                             gene.threshold = 3,
                             control.threshold = 3){
  #Input check
  if (!is.numeric(cell.threshold)){
    stop("Please set your Cell Threshold (NMAD value) to a valid integer.")
  }
  if(!is.numeric(control.threshold)){
    stop("Please set your Control Threshold (NMAD value) to a valid integer.")
  }
 
  # Stop if the user hasn't set any controls
  if(!progressLog(object)$controls){
    stop("Please define controls before filtering this dataset.")
  }
  
  # Gene ID
  gene_id_name <- colnames(rowInfo(object))[1]
  
  # Merge data together to make it easier to work with
  cell_info <- S4Vectors::merge(colInfo(object), 
                                SummarizedExperiment::colData(object), 
                                by = "cell_barcode")
  gene_info <- S4Vectors::merge(rowInfo(object),
                                SummarizedExperiment::rowData(object),
                                by = gene_id_name)
  
  # Reorder
  cell_info <- cell_info[match(colnames(object), cell_info$cell_barcode), ]
  gene_info <- gene_info[match(rownames(object), gene_info[, gene_id_name]), ]
  
  # Check values
  control_list <- progressLog(object)$set_controls
  
  total_counts <- cell_info$qc_libsize
  log10_total_counts <- log10(total_counts)
  total_features_counts_per_cell <- cell_info$qc_nfeaturecounts
  log10_features_counts_per_cell <- log10(total_features_counts_per_cell)
  pct_total_counts <- lapply(names(control_list), function(x) cell_info[, sprintf("qc_%s_pct_counts", x)])
  names(pct_total_counts) <- names(control_list)
  
  findOutliers <- function(values, 
                           nmads = 3, 
                           type = c("both", "lower", "upper"), 
                           na.rm = FALSE) {
    med_val <- stats::median(values, na.rm = na.rm)
    mad_val <- stats::mad(values, center = med_val, na.rm = na.rm)
    upper_limit <- med_val + nmads * mad_val
    lower_limit <- med_val - nmads * mad_val
    if (type == "lower"){
      upper_limit <- Inf
    } else if (type == "upper") {
      lower_limit <- -Inf
    }
    return(values < lower_limit | upper_limit < values)
  }
  
  ## Start identifying cells by barcodes
  # Returns boolean - do not merge yet
  cells_libsize_upper <- findOutliers(log10_total_counts, nmads=cell.threshold, type="lower") ## Remove cells with low expression
  cells_libsize_lower <- findOutliers(log10_total_counts, nmads=cell.threshold, type="upper") ## Remove cells with low expression
  cells_feature <- findOutliers(log10_features_counts_per_cell, nmads=control.threshold, type="lower") ## Remove cells with low number of genes
  cells_ngenes <- findOutliers(cell_info$qc_ngenes, nmads = gene.threshold, type = "lower")
  
  # Filter cells by MAD
  drop_barcodes_libsize <- c()
  if (any(cells_libsize_upper)){
    drop_barcodes_libsize <- c(drop_barcodes_libsize, which(cells_libsize_upper))
  }
  if (any(cells_libsize_lower)){
    drop_barcodes_libsize <- c(drop_barcodes_libsize, which(cells_libsize_lower))
  }
  
  if (length(drop_barcodes_libsize) > 0){
    drop_barcodes_libsize <- unique(drop_barcodes_libsize)
    cells_libsize <- rep(FALSE, length(cells_libsize_upper))
    cells_libsize[drop_barcodes_libsize] <- TRUE
  }
  
  # Filter features by MAD
  if (any(cells_feature)){
    drop_barcodes_feature <- which(cells_feature)
  } else {
    drop_barcodes_feature <- list()
  }
  
  if (any(cells_ngenes)){
    drop_barcodes_ngenes <- which(cells_ngenes)
  } else{
    drop_barcodes_ngenes <- list()
  }
  
  ## Identify cells to remove based on proportion of expression
  print("Identifying outliers...")
  control_counts <- BiocParallel::bplapply(pct_total_counts, findOutliers, nmads=control.threshold, type="higher") ## Use nmad to identify outliers
  drop_barcodes_controls <- BiocParallel::bplapply(control_counts, which) ## Identify cell barcodes to remove
  drop_barcodes_control_list <- BiocParallel::bplapply(names(drop_barcodes_controls), function(x) cell_info$cell_barcode[drop_barcodes_controls[[x]]])
  names(drop_barcodes_control_list) <- paste0("CellsFilteredBy", names(drop_barcodes_controls))
  
  ### Barcode master list of cells to remove
  remove_indices <- unique(c(as.vector(unlist(drop_barcodes_libsize)), 
                             as.vector(unlist(drop_barcodes_feature)),
                             as.vector(unlist(drop_barcodes_controls)),
                             as.vector(unlist(drop_barcodes_ngenes))))
  
  remove_barcodes <- colnames(object)[remove_indices]
  
  outlier_libsize_barcodes <- cell_info$cell_barcode[as.vector(unlist(drop_barcodes_libsize))]
  outlier_feature_barcodes <- cell_info$cell_barcode[as.vector(unlist(drop_barcodes_feature))]
  
  filtering_log <- list(CellsFilteredByLibSize = outlier_libsize_barcodes,
                        CellsFilteredByLowExpression = outlier_feature_barcodes)
  filtering_log <- c(filtering_log, drop_barcodes_control_list)
  
  # Remove cells from object
  filtered_object <- object[ , !(colnames(object) %in% remove_barcodes)]
  
  # Add records to log
  log <- progressLog(filtered_object)
  
  ### Update old log if it is still there
  if (!is.null(log$filterByOutliers)){
    old_filtering_log <- log$filterByOutliers
    keys <- unique(c(names(old_filtering_log), names(filtering_log)))
    filtering_log <- setNames(mapply(c, old_filtering_log[keys], filtering_log[keys]), keys)
  }
  
  filtering_df <- as.data.frame(t(as.matrix(sapply(names(filtering_log), function(x) length(filtering_log[[x]])))))
  
  # To go into the dataframe
  if (!is.null(log$FilteringLog)){
    old_filtering_df <- log$FilteringLog
    filtering_df <- dplyr::semi_join(filtering_log, old_filtering_df)
  }
  
  # Add to object
  log$filterByOutliers <- filtering_log
  log$FilteringLog <- filtering_df
  progressLog(filtered_object) <- log
  remove(object)
  return(filtered_object)
}
