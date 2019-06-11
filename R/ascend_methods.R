################################################################################
#
# ascend_methods.R
# description: Methods related to the updating and operation of EMSets
#
################################################################################

#' addGeneLabel
#' 
#' Labels cells if they express gene(s) as a condition.
#' 
#' @examples
#' em_set <- ascend::raw_set
#' em_set <- addGeneLabel(em_set, gene = c("MALAT1"))
#' 
#' @param x \linkS4class{EMSet}
#' @param gene List of gene markers
#' @return An EMSet with the expression of the specified genes noted in colInfo
#' 
#' @include ascend_objects.R
#' @importFrom SingleCellExperiment normcounts counts
#' @importClassesFrom S4Vectors DataFrame
#' @export
#' 
setGeneric("addGeneLabel", function(x, gene) 
  standardGeneric("addGeneLabel"))

#' @rdname addGeneLabel
setMethod("addGeneLabel", signature(x = "EMSet"), function(x, gene = c()){
  if (any(!(gene %in% rownames(x)))){
    stop("Please ensure all listed genes/transcripts are present in the matrix.")
  }
  
  # Use normcounts, otherwise use counts
  if ("normcounts" %in% names(SummarizedExperiment::assays(x))){
    counts <- SingleCellExperiment::normcounts(x)  
  } else(
    counts <- SingleCellExperiment::counts(x)
  )
  
  # Retrieve boolean of whether genes are present
  col_info <- colInfo(x)
  expressor_df <- as.matrix(counts[rownames(counts) %in% gene, ] > 0)
  
  # Flip it around if it's a vector of genes
  if (ncol(expressor_df) > 1){
    expressor_df <- t(expressor_df)
  }
  
  # Polish it to store in the colInfo slot
  expressor_df <- S4Vectors::DataFrame(expressor_df)
  colnames(expressor_df) <- gene
  # Add to col_info
  col_info <- cbind(col_info, expressor_df)
  colInfo(x) <- col_info 
  return(x)
})

calculateControlMetrics <- function(x, expression_matrix = NULL, gene_info = NULL, qc_cell_df = NULL){
  control_name <- x
  control_gene_info <- subset(gene_info, gene_info[, "control_group"] == control_name)
  control_counts <- expression_matrix[rownames(expression_matrix) %in% control_gene_info[ , 1], ]
  control_ntotal <- Matrix::colSums(control_counts)
  control_pct <- (control_ntotal / qc_cell_df$qc_libsize) * 100
  output <- list(control_ntotal = control_ntotal, control_pct = control_pct)
  return(output)
}

calcControlQC <- function(x, gene_info = NULL, qc_cell_df = NULL){
  # Group together - hopefully user hasn't annotated anything else in that column
  control_groups <- unique(gene_info[, "control_group"])
  control_groups <- control_groups[!is.na(control_groups)]
  qc_controls <- BiocParallel::bplapply(control_groups, 
                                        calculateControlMetrics, 
                                        expression_matrix = x,  
                                        gene_info = gene_info, 
                                        qc_cell_df = qc_cell_df)
  names(qc_controls) <- control_groups
  
  # Add controls to data frame
  for (control in control_groups){
    control_ntotal_name <- sprintf("qc_%s_ncounts", control)
    control_pct_name <- sprintf("qc_%s_pct_counts", control)
    ntotal <- qc_controls[[control]]$control_ntotal
    pct_total <- qc_controls[[control]]$control_pct
    qc_cell_df[[control_ntotal_name]] <- as.vector(ntotal)
    qc_cell_df[[control_pct_name]] <- as.vector(pct_total)
  }
  
  # Calculate proportion of non-control genes
  feature_matrix <- x[is.na(gene_info[,"control_group"]), ]
  n_featurecounts_per_cell <- Matrix::colSums(feature_matrix)
  pct_featurecounts_per_cell <- 100*(n_featurecounts_per_cell/Matrix::colSums(x))
  qc_cell_df$qc_nfeaturecounts <- n_featurecounts_per_cell
  qc_cell_df$qc_pctfeatures <- pct_featurecounts_per_cell
  
  return(qc_cell_df)
}

#' calculateQC
#' 
#' Calculates the following values for quality control:
#' 
#' \describe{
#'   \item{qc_libsize:}{Total number of counts per cell}
#'   \item{qc_ngenes}{Number of genes expressed by a cell}
#'   \item{qc_ncounts}{Total counts per gene} 
#'   \item{qc_ncells}{Number of cells expressing gene}
#'   \item{qc_meancounts}{Mean expression of the gene}
#'   \item{qc_topgeneranking}{Rank in most abundant gene expression}
#'   \item{qc_pct_total_expression}{Proportion of gene expression to total expression}
#' }
#' 
#' If users have defined controls, the following metrics are also calculated:
#' \describe{
#'   \item{qc_nfeaturecounts}{Number of reads from non-control genes (feature) for each cell}
#'   \item{qc_pctfeatures}{Percentage of features to total expression}
#' }
#' 
#' This function is called when changes are made to the count matrix.
#' 
#' @examples
#' # Load example dataset
#' object <- ascend::raw_set
#' 
#' # Recalculate quality control metrics
#' object <- calculateQC(object)
#' 
#' @param object An \linkS4class{EMSet}
#' @return QC values stored in an EMSet's rowData and colData slots
#' 
#' @include ascend_objects.R
#' @importFrom S4Vectors DataFrame merge
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @export
setGeneric("calculateQC", function(object) standardGeneric("calculateQC"))

#' @rdname calculateQC
setMethod("calculateQC", signature(object = "EMSet"), function(object){
  # Metrics - save trouble of retrieving object over and over
  expression_matrix <- SingleCellExperiment::counts(object)
  
  # Get old metadata
  old_rowData <- SummarizedExperiment::rowData(object)
  old_colData <- SummarizedExperiment::colData(object)
  old_colInfo <- colInfo(object)
  old_rowInfo <- rowInfo(object)
  gene_id_name <- colnames(old_rowInfo)[1]
  
  # Total expression for entire dataset
  qc_totalexpression <- sum(expression_matrix)
  
  # QC for cells
  # 1. qc_libsize: Total transcripts for each cell
  # 2. qc_ngenes: Number of genes expressed by each cell
  qc_libsize <- Matrix::colSums(expression_matrix)
  qc_ngenes <- Matrix::colSums(expression_matrix != 0)

  qc_cell_df <- data.frame(qc_libsize = qc_libsize, qc_ngenes = qc_ngenes)
  qc_cell_df$cell_barcode <- rownames(qc_cell_df)
  qc_cell_df <- qc_cell_df[match(colnames(expression_matrix), qc_cell_df$cell_barcode), c("cell_barcode", "qc_libsize", "qc_ngenes")]
  
  # QC for genes
  # 1. qc_genecounts: Total transcripts for each gene
  # 2. qc_cellspergene: Number of cells expressing each gene
  # 3. qc_meangenecounts: Average expression of each gene
  qc_genecounts <- Matrix::rowSums(expression_matrix)
  qc_cellspergene <- Matrix::rowSums(expression_matrix > 0)
  qc_meangenecounts <- Matrix::rowMeans(expression_matrix)
  qc_sortedcountspergene <- sort(qc_genecounts, decreasing = TRUE)
  
  qc_genecounts <- qc_genecounts[rownames(expression_matrix)]
  qc_cellspergene <- qc_cellspergene[rownames(expression_matrix)]
  qc_meangenecounts <- qc_meangenecounts[rownames(expression_matrix)]
  qc_generankings <- rank(-qc_genecounts, ties.method = "first")
  names(qc_generankings) <- rownames(expression_matrix)
  qc_pct_total_expression <- (qc_genecounts / qc_totalexpression) * 100
  
  qc_gene_df <- data.frame(qc_ncounts = qc_genecounts,
                           qc_ncells = qc_cellspergene,
                           qc_meancounts = qc_meangenecounts,
                           qc_topgeneranking = qc_generankings,
                           qc_pct_total_expression = qc_pct_total_expression)
  qc_gene_df[, gene_id_name] <- rownames(qc_gene_df)
  qc_gene_df <- qc_gene_df[match(rownames(expression_matrix), 
                                 qc_gene_df[ , gene_id_name]), 
                           c(gene_id_name, "qc_ncounts", "qc_ncells",
                             "qc_meancounts", "qc_topgeneranking", 
                             "qc_pct_total_expression")]
  
  # For control-related metrics
  if ("control_group" %in% colnames(old_rowInfo) && any(!is.na(old_rowInfo["control_group"]))){
    qc_cell_df <- calcControlQC(expression_matrix, gene_info = old_rowInfo, qc_cell_df = qc_cell_df)
  }
  
  # Sort again
  qc_cell_df <- qc_cell_df[match(colnames(expression_matrix), qc_cell_df$cell_barcode), ]
  qc_gene_df <- qc_gene_df[match(rownames(expression_matrix), as.vector(qc_gene_df[, gene_id_name])), ]

  # Convert to DataFrame
  qc_cell_df <- S4Vectors::DataFrame(qc_cell_df)
  qc_gene_df <- S4Vectors::DataFrame(qc_gene_df)
  
  # Merge with old data
  if (any(colnames(qc_cell_df)[2:length(qc_cell_df)] %in% colnames(old_colData))){
    new_colData <- mergeDF(old_colData, qc_cell_df, "cell_barcode")
  } else{
    new_colData <- S4Vectors::merge(old_colData, qc_cell_df, by = "cell_barcode")
  }
  
  if (any(colnames(qc_gene_df)[2:length(qc_gene_df)] %in% colnames(old_rowData))){
    new_rowData <- mergeDF(old_rowData, qc_gene_df, gene_id_name)
  } else{
    new_rowData <- S4Vectors::merge(old_rowData, qc_gene_df, by = gene_id_name)
  }

  # Re-order
  new_colData <- new_colData[match(colnames(expression_matrix), new_colData$cell_barcode), ]
  BiocGenerics::rownames(new_colData) <- new_colData[,1]
  new_rowData <- new_rowData[match(rownames(expression_matrix), as.vector(new_rowData[, gene_id_name])), ]
  BiocGenerics::rownames(new_rowData) <- new_rowData[ , gene_id_name]
  
  # Replace information
  SummarizedExperiment::colData(object) <- new_colData
  SummarizedExperiment::rowData(object) <- new_rowData
  return(object)
})

#' convertGeneID
#' 
#' Switches gene identifiers used in the EMSet to the identifiers stored in a 
#' column of rowInfo.
#' 
#' @examples
#' # Load data
#' x <- ascend::raw_set
#' 
#' # Switch identifiers
#' x <- convertGeneID(x, new.annotation = "ensembl_gene_id")
#' 
#' @param object \linkS4class{EMSet}
#' @param new.annotation Name of the column where the new set of identifiers 
#' are stored in rowInfo
#' 
#' @return An EMSet with genes swapped to identifiers stored in specified column
#' of colInfo.
#' 
#' @include ascend_objects.R
#' @include ascend_getters.R
#' @include ascend_setters.R
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @export
setGeneric("convertGeneID", function(object, new.annotation) 
  standardGeneric("convertGeneID"))

#' @rdname convertGeneID
setMethod("convertGeneID", signature(object = "EMSet"), function(object,
                                             new.annotation = NULL){
  # Get old information
  row_info <- rowInfo(object)
  row_data <- rowData(object)
  row_info_names <- colnames(row_info)
  row_data_names <- colnames(row_data)
  old.annotation <- row_info_names[1]
  
  # Check if the new annotation exists
  if(!new.annotation %in% row_info_names){
    stop("Your selected gene annotation is not present in the rowInfo dataframe.")
  } else{
    if (identical(old.annotation, new.annotation)){
      stop("Your selected gene annotation is already in use.")
    }
  }
  
  # Identify non-identifier columns
  other_columns_info <- row_info_names[which(!(row_info_names %in% 
                                                 c(old.annotation, new.annotation)))]
  other_columns_data <- row_data_names[which(!(row_data_names %in% 
                                                 c(old.annotation, new.annotation)))]
  
  # New order for rowInfo
  reordered_columns_info <- c(new.annotation, old.annotation, other_columns_info)
  reordered_rowInfo <- row_info[ , reordered_columns_info]
  new_annotations <- reordered_rowInfo[, 1]
  rownames(reordered_rowInfo) <- new_annotations
  
  # Replace annotation in rowData
  reordered_rowData <- row_data
  colnames(reordered_rowData)[1] <- new.annotation
  reordered_rowData[ ,1] <- new_annotations
  rownames(reordered_rowData) <- new_annotations
  
  # New identifiers
  renamed_set <- object
  
  # Check if any controls are defined
  # If there are controls, convert them to the new convention
  log <- progressLog(object)
  if (log$controls){
    control_list <- log$set_controls
    converted_controls <- lapply(names(control_list), function(x){
      old_gene_ids <- control_list[[x]]
      new_gene_ids <- reordered_rowInfo[which(reordered_rowInfo[ , old.annotation] 
                                              %in% old_gene_ids), new.annotation]
      return(new_gene_ids)
    })
    names(converted_controls) <- names(control_list)
    log$set_controls <- converted_controls
  }
  
  # Replace slots
  progressLog(renamed_set) <- log
  BiocGenerics::rownames(renamed_set) <- new_annotations
  rowData(renamed_set) <- DataFrame(reordered_rowData)
  rowInfo(renamed_set) <- DataFrame(reordered_rowInfo)
  
  # Recalculate QC
  renamed_set <- calculateQC(renamed_set)
  return(renamed_set)  
})

#' calculateCV
#' 
#' Calculates the Coefficient of Variation (CV) for each gene using normalised
#' values. Results are stored in the rowData slot.
#' 
#' @examples
#' # Load example dataset
#' x <- ascend::analyzed_set
#' 
#' # Calculate CV 
#' x <- calculateCV(x)
#' 
#' # Show CV results
#' rowData(x)
#' 
#' @param object An \linkS4class{EMSet}
#' @return An EMSet with calculated CV values.
#' 
#' @include ascend_objects.R
#' @include ascend_getters.R
#' @include ascend_setters.R
#' @importFrom SummarizedExperiment rowData
#' @export
#' 
setGeneric("calculateCV", function(object) standardGeneric("calculateCV"))

#' @rdname calculateCV
setMethod("calculateCV", signature(object = "EMSet"), function(object){
  # Check if data has been normalised
  if (!("normcounts" %in% SummarizedExperiment::assayNames(object))){
    stop("Please normalise your data before using this function.")
  }
  
  # Calculate CV
  expression_matrix <- SingleCellExperiment::normcounts(object)
  metrics <- as.data.frame(SummarizedExperiment::rowData(object))

  # Calculate standard deviation of gene
  sd_genes <- apply(expression_matrix, 1, stats::sd)
  
  # Coefficient of variation
  cv_genes <- (sd_genes/metrics$qc_meancounts)
  
  # Add information to rowData
  metrics$ascend_cv <- as.vector(unlist(cv_genes))
  metrics$ascend_cv_rank <- rank(-metrics$ascend_cv, ties.method = "first")
  metrics$ascend_sd <- sd_genes
  
  # Update rowData and return information
  SummarizedExperiment::rowData(object) <- S4Vectors::DataFrame(metrics)
  return(object)
})
