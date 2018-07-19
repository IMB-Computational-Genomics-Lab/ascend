################################################################################
#
# ascend_methods.R
# description: Methods related to the updating and operation of EMSets
#
################################################################################

#' @export
setGeneric("addControlInfo", function(x, ..., controls) standardGeneric("addControlInfo"))

#' addControlInfo
#' 
#' Adds and integrates control information to an \linkS4class{EMSet}.
#' 
#' @param x An \linkS4class{EMSet}.
#' @param controls A named list of control groups and controls.
#' @return An \linkS4class{EMSet} with controls integrated into \code{rowInfo}
#' and calculated QC metrics.
#'
#' @include ascend_objects.R
#' @include ascend_setters.R
#' @include ascend_getters.R
#' @export
#' @importFrom SummarizedExperiment rowData
setMethod("addControlInfo", "EMSet" , function(x, controls = NULL){
  # Get row data
  row_info <- rowInfo(x)
  
  # Set control group to NULL by default. This is for non-control genes
  row_info$control_group <- NA
  
  # For each control group...
  # Control group check in case user did not group the controls
  if ((length(names(controls)) > 1) & (length(names(controls)) != nrow(row_info))){
    for (control_name in names(controls)){
      gene_set <- controls[[control_name]]
      row_info$control_group[which(row_info[ ,1] %in% gene_set)] <- control_name
    }
  } else{
    row_info$control_group[which(row_info[,1] %in% unlist(controls))] <- "Control"
  }
  
  # Replace control information
  rowInfo(x) <- row_info
  
  # Get log information
  log <- progressLog(x)
  
  # Update log information with controls
  log$set_controls <- controls
  log$controls <- TRUE
  progressLog(x) <- log
  return(x)
})

#' @include ascend_objects.R
# Per control - add to cell information
#' @export
calculateControlMetrics <- function(x, expression_matrix = NULL, gene_info = NULL, qc_cell_df = NULL){
  control_name <- x
  control_gene_info <- subset(gene_info, gene_info[, "control_group"] == control_name)
  control_counts <- expression_matrix[which(rownames(expression_matrix) %in% control_gene_info[ , 1]), ]
  control_ntotal <- colSums(control_counts)
  control_pct <- (control_ntotal / qc_cell_df$qc_libsize) * 100
  output <- list(control_ntotal = control_ntotal, control_pct = control_pct)
  return(output)
}

#' @include ascend_objects.R
#' @export
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
  feature_matrix <- x[which(is.na(gene_info[,"control_group"])), ]
  n_featurecounts_per_cell <- colSums(feature_matrix)
  pct_featurecounts_per_cell <- 100*(n_featurecounts_per_cell/colSums(x))
  qc_cell_df$qc_nfeaturecounts <- n_featurecounts_per_cell
  qc_cell_df$qc_pctfeatures <- pct_featurecounts_per_cell
  
  return(qc_cell_df)
}


#' @export
setGeneric("calculateQC", function(object, ...) standardGeneric("calculateQC"))


#' @include ascend_objects.R
#' @importFrom S4Vectors DataFrame merge
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @export
setMethod("calculateQC", "EMSet", function(object){
  # Metrics - save trouble of retrieving object over and over
  expression_matrix <- SingleCellExperiment::counts(object)
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
  qc_libsize <- colSums(expression_matrix)
  qc_ngenes <- colSums(expression_matrix != 0)
  
  # QC for genes
  # 1. qc_genecounts: Total transcripts for each gene
  # 2. qc_cellspergene: Number of cells expressing each gene
  # 3. qc_meangenecounts: Average expression of each gene
  qc_genecounts <- rowSums(expression_matrix)
  qc_cellspergene <- rowSums(expression_matrix > 0)
  qc_meangenecounts <- rowMeans(expression_matrix)
  qc_sortedcountspergene <- sort(qc_genecounts, decreasing = TRUE)
  
  qc_genecounts <- qc_genecounts[rownames(expression_matrix)]
  qc_cellspergene <- qc_cellspergene[rownames(expression_matrix)]
  qc_meangenecounts <- qc_meangenecounts[rownames(expression_matrix)]
  qc_generankings <- match(rownames(expression_matrix), names(qc_sortedcountspergene))
  names(qc_generankings) <- rownames(expression_matrix)
  qc_pct_total_expression <- (qc_genecounts / qc_totalexpression) * 100
  
  qc_cell_df <- data.frame(qc_libsize = qc_libsize, qc_ngenes = qc_ngenes)
  qc_cell_df$cell_barcode <- rownames(qc_cell_df)
  qc_cell_df <- qc_cell_df[match(colnames(expression_matrix), qc_cell_df$cell_barcode), c("cell_barcode", "qc_libsize", "qc_ngenes")]
  
  qc_gene_df <- data.frame(qc_ncounts = qc_genecounts,
                           qc_ncells = qc_cellspergene,
                           qc_meancounts = qc_meangenecounts,
                           qc_topgeneranking = qc_generankings,
                           qc_pct_total_expression = qc_pct_total_expression)
  qc_gene_df[, gene_id_name] <- rownames(qc_gene_df)
  qc_gene_df <- qc_gene_df[match(rownames(expression_matrix), qc_gene_df[ , gene_id_name]), c(gene_id_name, "qc_ncounts", "qc_ncells",
                                                                                     "qc_meancounts", "qc_topgeneranking", "qc_pct_total_expression")]
  # For control-related metrics
  if ("control_group" %in% colnames(old_rowInfo) && any(!is.na(old_rowInfo["control_group"]))){
    qc_cell_df <- calcControlQC(expression_matrix, gene_info = old_rowInfo, qc_cell_df = qc_cell_df)
  }
  
  # Sort again
  qc_cell_df <- qc_cell_df[match(colnames(expression_matrix), qc_cell_df$cell_barcode), ]
  qc_gene_df <- qc_gene_df[match(rownames(expression_matrix), qc_gene_df[, gene_id_name]), ]

    # Convert to DataFrame
  qc_cell_df <- S4Vectors::DataFrame(qc_cell_df)
  qc_gene_df <- S4Vectors::DataFrame(qc_gene_df)
  
  # We'll just replace pre-existing colnames
  if (all(colnames(qc_cell_df) %in% colnames(old_colData))){
    new_colData <- old_colData
    new_colData$qc_libsize <- qc_cell_df$qc_libsize
    new_colData$qc_ngenes <- qc_cell_df$qc_ngenes
  } else{
    # Identify common columns
    common_columns <- intersect(colnames(old_colData), colnames(qc_cell_df))
    # Merge dataframes
    new_colData <- S4Vectors::merge(old_colData, qc_cell_df, by = common_columns)    
  }
  
  # Re-order
  new_colData <- new_colData[match(colnames(expression_matrix), new_colData[ ,1]), ]
  BiocGenerics::rownames(new_colData) <- new_colData[,1]
  
  if (all(colnames(qc_gene_df) %in% colnames(old_rowData))){
    new_rowData <- old_rowData
    new_rowData$qc_ncells <- qc_gene_df$qc_ncells
    new_rowData$qc_meancounts <- qc_gene_df$qc_meancounts
    new_rowData$qc_topgeneranking <- qc_gene_df$qc_topgeneranking
    new_rowData$qc_pct_total_expression <- new_rowData$qc_pct_total_expression
    
  } else{
    # Merge dataframes
    new_rowData <- S4Vectors::merge(old_rowData, qc_gene_df, by = 1)
  }

  # Re-order
  new_rowData <- new_rowData[match(rownames(expression_matrix), new_rowData[, gene_id_name]), ]
  BiocGenerics::rownames(new_rowData) <- new_rowData[ , gene_id_name]
  
  # Replace information
  SummarizedExperiment::colData(object) <- new_colData
  SummarizedExperiment::rowData(object) <- new_rowData
  
  return(object)
})
