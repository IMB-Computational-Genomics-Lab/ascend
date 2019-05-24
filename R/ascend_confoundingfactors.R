regressGene <- function(x, covariate_matrix = NULL, counts = NULL) {
  regress_matrix <- cbind(counts, covariate_matrix[ ,x])
  colnames(regress_matrix) <- c(gsub("\\-", "_", colnames(counts)), "GENE")
  fmla <- stats::as.formula(paste("GENE ", " ~ ", paste(gsub("\\-", "_", 
                                                             colnames(covariate_matrix)), 
                                                        collapse = "+"), sep = ""))
  # Fit values to gene
  resid_gene <- stats::lm(fmla, data = regress_matrix)
  residuals <- as.data.frame(stats::residuals(resid_gene))
  return(residuals)
}

#' regressConfoundingFactors
#'
#' This function generates a scaled regression matrix based on candidate genes
#' supplied by the user. This function should be used after normalisation.
#'  
#' @param object An \code{\linkS4class{EMSet}} that has been normalised.
#' @param candidate.genes A list of genes you wish to regress from the dataset.
#' Refer to the vignette on how to choose genes for regression.
#' 
#' @return An \code{\linkS4class{EMSet}} with confounding factors regressed 
#' from the expression values.
#' 
#' @examples
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Define genes to regress
#' genes <- c("CDK4","CCND1")
#' 
#' regressedSet <- regressConfoundingFactors(em_set, candidate.genes = genes)
#' 
#' @importFrom Matrix t
#' @importFrom BiocParallel bplapply
#' @importFrom data.table as.data.table
#' @export
#' 
regressConfoundingFactors <- function(object, candidate.genes = c()) {
  # Check if EMSet
  if (!is(object, "EMSet")){
    stop("Please supply an EMSet.")
  }
  
  # Check if genes are supplied
  if (!(length(candidate.genes) > 0)) {
    stop("Please supply a list of gene candidates.")
  }
  
  # Check if data has been normalised
  if (!("normcounts" %in% SummarizedExperiment::assayNames(object))){
    stop("Data has not been normalised. Please normalise your data first.")
  }
  
  # Check that candidate genes are present
  if (!(any(candidate.genes %in% rownames(object)))){
    stop("No genes supplied in candidate.genes are present.")
  } else{
    candidate.genes <- candidate.genes[which(candidate.genes %in% rownames(object))]
  }
  
  # Insert checks for data type
  counts <- logcounts(object)
  cell_barcodes <- colnames(counts)
  gene_ids <- rownames(counts)
  
  # Original type
  object_type <- is(counts)
  
  if (object_type[1] == "dgeMatrix"){
    object_type <- "dgCMatrix"
  }
  
  # Weird matrix handling
  if (is(counts, "dgeMatrix")){
    counts <- as(counts, "dgCMatrix")
    counts <- Matrix::t(counts)
    counts <- as.matrix(counts)
    counts <- as.data.frame(counts)
  }
  
  # If sparse matrix
  if (is(counts, "dgCMatrix")){
    counts <- Matrix::t(counts)
    counts <- as.matrix(counts)
    counts <- as.data.frame(counts)
  }
  
  # If matrix
  if (is.matrix(counts)){
    counts <- t(counts)
    counts <- as.data.frame(counts)
  }
  
  if (length(candidate.genes) > 1){
    covariate_matrix <- counts[, candidate.genes]    
  } else{
    covariate_matrix <- as.data.frame(counts[, candidate.genes])
    colnames(covariate_matrix) <- candidate.genes
  }

  
  # Calculate residuals
  residuals <- BiocParallel::bplapply(candidate.genes,
                                      regressGene,
                                      covariate_matrix = covariate_matrix,
                                      counts = counts)
  
  residual_df <- do.call("cbind", residuals)
  colnames(residual_df) <- candidate.genes
  
  # Scale data
  scaled_residuals <- scale(residual_df, center = TRUE, scale = TRUE)
  counts <- as.matrix(counts)
  confounder_mat <- counts[, candidate.genes] - scaled_residuals
  counts[rownames(confounder_mat), colnames(confounder_mat)] <- confounder_mat
  counts <- counts[colnames(object), rownames(object)]
  counts <- Matrix::t(counts)
  
  # Update object
  regcounts(object) <- counts
  log <- progressLog(object)
  log$regressConfoundingFactors <- list(candidate.genes = candidate.genes, regressed = TRUE)
  progressLog(object) <- log
  
  return(object)
}
