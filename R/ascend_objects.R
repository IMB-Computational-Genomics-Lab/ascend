################################################################################
#
# ascend_objects.R
# description: Code related to the creation and maintenance of the EMSet
#
################################################################################ 
#' validateEMSet
#'
#' Validation function for the EMSet.
#' 
#' @param object An EMSet that has been created, but needs validation.
#' @keywords internal
#' @export
validateEMSet <- function(object){
  # Store errors in an object. If length is greater than 0, indicates multiple
  # errors are present.
  errors <- c()
  
  # 1. Retrieve cell identifiers and gene names
  cell_list <- BiocGenerics::colnames(object)
  gene_list <- BiocGenerics::rownames(object)
  
  # 2. Retrieve data from colInfo and rowInfo slots
  col_info <- object@colInfo
  row_info <- object@rowInfo
  
  # 3. Check they match the rownames of the slotted items
  if (any(cell_list != rownames(col_info))){
    errors <- c(errors, "Cell identifiers in count matrix do not match rownames of colInfo dataframe.")
  }
  if (any(gene_list != rownames(row_info))){
    errors <- c(errors, "Gene identifiers in count matrix do not match rownames of rowInfo dataframe.")
  }
  
  # 4. Check if they match the first column of the slotted items
  if (any(cell_list != col_info[, 1])){
    errors <- c(errors, "Cell identifiers in count matrix do not match the first column of colInfo dataframe.")
  }
  if (any(gene_list != row_info[, 1])){
    errors <- c(errors, "Gene identifiers in count matrix do not match the first column of rowInfo dataframe.")
  }
  
  # If any messages were collected, return an error, otherwise return TRUE
  if (length(errors)){
    errors
  } else TRUE
}

#' Expression and Metadata Set (EMSet)
#'
#' An \linkS4class{EMSet} is a S4 class that stores data in a format `ascend`` can 
#' easily work with for analysis.
#' 
#' This iteration of the \linkS4class{EMSet} inherits from the 
#' \code{\linkS4class{SingleCellExperiment}} superclass for integration into
#' Bioconductor.  
#' 
#' `ascend`-specific slots are as follows:
#' @slot colInfo A data frame containing each cell identifier, its 
#' associated batch/sample and additional information such as conditions. This 
#' slot was called CellInformation in earlier versions of `ascend`.
#' @slot rowInfo A data frame containing information a set of gene 
#' identifiers, such as gene symbols or ENSEMBL transcript identifiers. This 
#' data frame also holds information on controls and any information provided by 
#' the user. This slot was called GeneInformation in earlier versions of `ascend`.
#' @slot clusterAnalysis Objects related to clustering - distance matrix and
#' hclust object
#' @slot log A record of functions used on an \linkS4class{EMSet}.
#' 
#' The following data is loaded into \linkS4class{SingleCellExperiment} slots:
#' @slot assays Raw expression matrix is loaded into the \code{counts} slot.
#' Normalised counts are loaded into the following based on their normalisation
#' method:
#'  \describe{
#'   \item{normcounts}{Counts normalised by `ascend`'s RLE method.}
#'   \item{logcounts}{Counts normalised by `scran`'s deconvolution method.}
#' }
#' @slot colData The following cell-based QC metrics are stored in this 
#' DataFrame:
#' \describe{
#' \item{qc_libsize}{Total number of transcripts per cell.}
#' \item{qc_ngenes}{Total number of expressed genes per cell.}
#' \item{qc_control_expression}{Percentage of control expression to total 
#' expression.}
#' }
#' @slot rowData The following gene-based QC metrics are stored in this
#' DataFrame:
#' \describe{
#' \item{qc_ncounts}{Total gene expression across dataset.}
#' \item{qc_ncells}{Number of cells the gene is expressed in.}
#' \item{qc_meancounts}{Mean expression level of the gene.}
#' \item{qc_topgeneranking}{Position of gene in top gene expression rankings.}
#' \item{qc_pct_total_expression}{Proportion of gene expression to total
#' expression.}
#' }
#' @name EMSet-class
#' @rdname EMSet-class
#' @export
#' @exportClass EMSet
#' 
.EMSet <- setClass("EMSet", 
                       slots = list(colInfo = "DataFrame", 
                                    rowInfo = "DataFrame", 
                                    clusterAnalysis = "list", 
                                    log = "list"), 
                       contains = "SingleCellExperiment", 
                       prototype = prototype(new("SingleCellExperiment")),
                       validity = validateEMSet)

#' @importFrom SummarizedExperiment assayNames rowData colData
#' @importFrom S4Vectors metadata
#' @importFrom BiocGenerics rownames colnames
#' @importFrom SingleCellExperiment reducedDimNames spikeNames
#' @aliases EMSet
#' @export
#' 
setMethod("show", "EMSet", function(object){
  # 1. Basic EMSet structure - all will have these features
  set_dims <- dim(object)
  metadata_names <- paste(names(S4Vectors::metadata(object)), collapse = ", ")
  assay_names <- paste(SummarizedExperiment::assayNames(object), collapse = ", ")
  rownames <- paste(BiocGenerics::rownames(object)[1:10], collapse = ", ")
  rowInfo <- paste(colnames(rowInfo(object)), collapse = ", ")
  rowData <- paste(colnames(SummarizedExperiment::rowData(object)), collapse =", ")
  colnames <- paste(BiocGenerics::colnames(object)[1:10], collapse = ", ")
  colInfo <- paste(colnames(colInfo(object)), collapse = ", ")
  colData <- paste(colnames(SummarizedExperiment::colData(object)), collapse = ", ")
  dimnames <- paste(SingleCellExperiment::reducedDimNames(object), collapse = ", ")
  spikenames <- paste(SingleCellExperiment::spikeNames(object), collapse = ", ")
  clusteranlaysis <- paste(names(clusterAnalysis(object)), collapse = ", ")
  
  line1 <- "class: EMSet"
  line2 <- sprintf("dim: %i %i", set_dims[1], set_dims[2])
  line3 <- print(paste0("metadata:", metadata_names))
  line4 <- paste("assays:", assay_names, collapse = " ")
  line5 <- paste("rownames:", rownames, "...", collapse = " ")
  line6 <- paste("rowInfo:", rowInfo, collapse = " ")
  line7 <- paste("rowData:", rowData, collapse = " ")
  line8 <- paste("colnames:", colnames, "...")
  line9 <- paste("colInfo:", colInfo, collapse = " ")
  line10 <- paste("colData:", colData, collapse = " ")
  line11 <- paste("reducedDimNames:", dimnames, collapse = " ")
  line12 <- paste("spikeNames:", spikenames, collapse = " ")
  line13 <- paste("clusterAnalysis:", clusteranlaysis, collapse = " ")
  cat(paste(line1, line2, line3, line4, line5, line6, line7,
              line8, line9, line10, line11, line12, line13, sep = "\n"))
})

#' newEMSet
#' 
#' \code{\link{newEMSet}} generates a \linkS4class{EMSet} object for use with 
#' the `ascend` package. This object contains an expression matrix, 
#' associated metadata, downstream analysis and a log documenting the actions 
#' used to shape the data in this object.
#' 
#' @param assays A list of matrices representing count data. You must have a 
#' count matrix stored under "counts" in this list.
#' @param colInfo A data frame containing cell-related metadata.
#' @param rowInfo A data frame containing gene-related metadata.
#' @param colData A data frame containing cell-related data.
#' @param rowData A data frame containing gene-related data.
#' @param controls A named list containing control genes grouped into named 
#' control groups.
#' 
#' @examples
#' # Randomly generate count matrix
#' count_matrix <- matrix(sample(0:1, 100, replace=TRUE),10,10)
#' 
#' # Generate cell barcodes
#' cell_barcodes <- paste0("Cell-", 1:10)
#' gene_ids <- paste0("Gene-", 1:10)
#' 
#' # Add to matrix
#' colnames(count_matrix) <- cell_barcodes
#' rownames(count_matrix) <- gene_ids
#' 
#' # Create an EMSet
#' em_set <- newEMSet(assays = list(counts = count_matrix))
#' 
#' 
#' @return An \linkS4class{EMSet}.
#' 
#' @export
newEMSet <- function(assays = NULL, 
                     colInfo = NULL,
                     rowInfo = NULL,
                     colData = NULL,
                     rowData = NULL,
                     controls = NULL){
  if (!("counts" %in% names(assays))){
    stop("Please ensure 'counts' are in your supplied assay list.")
  } else{
    counts <- assays$counts
  }
  
  
  if (is.null(colInfo)){
    colInfo <- S4Vectors::DataFrame(cell_barcode = colnames(counts), batch = rep(1, ncol(counts)), row.names = colnames(counts))  
  } else{
    if(!("cell_barcode" %in% colnames(colInfo)[1])){
      stop("Please specify the name of the first column in colInfo as 'cell_barcode'")
    }
  }
  
  if (is.null(rowInfo)){
    rowInfo <- S4Vectors::DataFrame(gene_id = rownames(counts), row.names = rownames(counts))
  }
  
  # Prime rowData and colData with identifiers
  if (is.null(colData)){
    colData <- S4Vectors::DataFrame(cell_barcode = colInfo[, 1], row.names = colInfo[,1])
  }
  
  if (is.null(rowData)){
    # Ensure we use the same header as other datasets
    gene_id_name <- colnames(rowInfo)[1]
    gene_id_list <- list()
    gene_id_list[[gene_id_name]] <- rowInfo[ ,1]
    rowData <- S4Vectors::DataFrame(gene_id_list, row.names = rowInfo[, 1])
  } else{
    # Check the columns match
    if (!identical(colnames(rowData)[1], colnames(rowInfo)[1])){
      stop("Please ensure the name of the first column of rowInfo is the
           same as the first column of rowData.")
    }
  }
  
  # Load metadata as required
  colInfo <- S4Vectors::DataFrame(colInfo, row.names = colInfo[, 1])
  rowInfo <- S4Vectors::DataFrame(rowInfo, row.names = rowInfo[, 1])
  
  # Load into object
  object <- new("EMSet", 
                SingleCellExperiment::SingleCellExperiment(assays = assays), 
                colInfo = colInfo, 
                rowInfo = rowInfo)
  
  # Append rowData and colData
  SummarizedExperiment::colData(object) <- colData
  SummarizedExperiment::rowData(object) <- rowData
  
  # Add controls if present
  if (!is.null(controls)){
    object <- addControlInfo(object, controls = controls)
  }
  
  # Calculate QC
  BiocGenerics::colnames(object) <- colnames(SingleCellExperiment::counts(object))
  object <- calculateQC(object)
  
  # Return to user
  return(object)
}

#' @export
setGeneric("EMSet2SCE", function(x, ...) standardGeneric("EMSet2SCE"))

#' EMSet2SCE
#' 
#' Converts an EMSet to a SingleCellExperiment object. This function preserves
#' the information stored in the colInfo and rowInfo steps by merging them with
#' colData and rowData, respectively. Users should extract information in
#' the log and clusterAnalysis slots prior to conversion.
#' 
#' @param x An EMSet
#' @return A SingleCellExperiment object.
#' 
#' @examples
#' # Randomly generate count matrix
#' count_matrix <- matrix(sample(0:1, 100, replace=TRUE),10,10)
#' 
#' # Generate cell barcodes
#' cell_barcodes <- paste0("Cell-", 1:10)
#' gene_ids <- paste0("Gene-", 1:10)
#' 
#' # Add to matrix
#' colnames(count_matrix) <- cell_barcodes
#' rownames(count_matrix) <- gene_ids
#' 
#' # Create an EMSet
#' em_set <- newEMSet(assays = list(counts = count_matrix))
#' 
#' # Convert to SingleCellExperiment
#' single_cell_experiment <- EMSet2SCE(em_set)
#' 
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors merge
#' @importFrom BiocGenerics rownames colnames
#' @export
setMethod("EMSet2SCE", "EMSet", function(x){
  # Convert to SingleCellExperiment
  col_info <- colInfo(x)
  row_info <- rowInfo(x)
  
  # Save the name so we can use it
  gene_id_name <- names(col_info)[1]
  log <- progressLog(x)
  cluster_analysis <- clusterAnalysis(x)
  
  # Coerce into single cell experiment
  single_cell_set <- as(x, "SingleCellExperiment")
  
  # Stash data into metadata slots
  metadata_list <- S4Vectors::metadata(single_cell_set)
  metadata_list$EMSet <- TRUE
  metadata_list$ascend_log <- log
  metadata_list$ascend_cluster_analysis <- cluster_analysis
  metadata_list$ascend_colinfo_headers <- colnames(col_info)
  metadata_list$ascend_rowinfo_headers <- colnames(row_info)
  metadata_list$ascend_coldata_headers <- colnames(SummarizedExperiment::colData(x))
  metadata_list$ascend_rowdata_headers <- colnames(SummarizedExperiment::rowData(x))
  S4Vectors::metadata(single_cell_set) <- metadata_list
  
  col_data <- SummarizedExperiment::colData(single_cell_set)
  row_data <- SummarizedExperiment::rowData(single_cell_set)
  colnames <- BiocGenerics::colnames(SingleCellExperiment::counts(single_cell_set))
  rownames <- BiocGenerics::rownames(SingleCellExperiment::counts(single_cell_set))
  updated_col_data <- S4Vectors::merge(col_info, col_data, by = 1)
  updated_row_data <- S4Vectors::merge(row_info, row_data, by = 1)
  updated_col_data <- updated_col_data[match(colnames, updated_col_data$cell_barcode), ]
  updated_row_data <- updated_row_data[match(rownames, updated_row_data[, 1]), ]
  BiocGenerics::rownames(updated_col_data) <- colnames
  BiocGenerics::rownames(updated_row_data) <- rownames
  SummarizedExperiment::colData(single_cell_set) <- updated_col_data
  SummarizedExperiment::rowData(single_cell_set) <- updated_row_data
  return(single_cell_set)
})

#' @export
setGeneric("SCE2EMSet", function(x, ...) standardGeneric("SCE2EMSet"))

#' SCE2EMSet
#' 
#' Converts an SingleCellExperiment object back into an EMSet. This function 
#' rebuilds the slots from data kept in the SingleCellExperiment. This should
#' not be used to initiate a new EMSet from a SingleCellExperiment.
#' 
#' @param x SingleCellExperiment object
#' @return An EMSet
#' 
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors merge
#' @importFrom BiocGenerics rownames colnames
#' @export
setMethod("SCE2EMSet", "SingleCellExperiment", function(x){
  # Extract data from SCE
  sce_col_data <- SummarizedExperiment::colData(x)
  sce_row_data <- SummarizedExperiment::rowData(x)
  metadata <- S4Vectors::metadata(x)
  
  if (is.null(metadata$EMSet)){
    stop("Please supply a SingleCellExperiment that was generated from an EMSet.")
  }
  
  # Extract parent data
  sce_col_info <- sce_col_data[ , colnames(sce_col_data) %in% metadata$ascend_colinfo_headers]
  sce_row_info <- sce_row_data[ , colnames(sce_row_data) %in% metadata$ascend_rowinfo_headers]
  sce_col_data <- sce_col_data[, colnames(sce_col_data) %in% metadata$ascend_coldata_headers]
  sce_row_data <- sce_row_data[, colnames(sce_row_data) %in% metadata$ascend_rowdata_headers]
  rownames(sce_col_info) <- sce_col_info$cell_barcode
  rownames(sce_row_info) <- sce_row_info[, 1]
  
  SummarizedExperiment::colData(x) <- sce_col_data
  SummarizedExperiment::rowData(x) <- sce_row_data
  
  # Create an EMSet
  em_set <- new("EMSet", x, colInfo = sce_col_info, rowInfo = sce_row_info, 
                log = metadata$ascend_log, clusterAnalysis = metadata$ascend_cluster_analysis)
  em_set <- calculateQC(em_set)
  return(em_set)
})

#' updateObject
#' 
#' Converts EMSets from ascend version < 0.5.0 to new EMSet that inherits from
#' SingleCellExperiment. Please note that normalised data will be used as the 
#' main count matrix, and loaded into the normcounts slot if a normalisation
#' method was logged. Quality control metrics will be re-calculated upon 
#' conversion.
#' 
#' @param object An old EMSet.
#' @examples
#' \dontrun{
#' # Make sure you replace your original object with the updated object
#' object <- updateObject(object)
#' }
#' @return An updated EMSet.
#' 
#' @export
setMethod("updateObject", "EMSet", function(object){
  if (.hasSlot(object, "ExpressionMatrix")){
    print("Old EMSet detected! Updating to new EMSet structure...")
    log <- object@Log
    
    counts <- object@ExpressionMatrix
    
    # Retrieve metadata  
    cell_info <- object@CellInformation
    gene_info <- object@GeneInformation
    
    rownames(cell_info) <- colnames(counts)
    rownames(gene_info) <- rownames(counts)
    
    # Create base EMSet
    updated_object <- newEMSet(assays = list(counts = as.matrix(counts)),
                               colInfo = S4Vectors::DataFrame(cell_info),
                               rowInfo = S4Vectors::DataFrame(gene_info))
    
    SummarizedExperiment::mcols(updated_object) <- S4Vectors::DataFrame(gene_info)
    
    if(length(object@Controls) > 0){
      updated_object <- addControlInfo(updated_object, controls = object@Controls)
    }
    
    progressLog(updated_object) <- log
    
    # Check if values are normalised
    if (!is.null(object@Log$NormalisationMethod)){
      normcounts <- object@ExpressionMatrix
      SingleCellExperiment::normcounts(updated_object) <- normcounts
    }
    
    if (length(object@PCA) > 0){
      pca_matrix <- object@PCA$PCA
      pct_variance <- object@PCA$PCAPercentVariance
      SingleCellExperiment::reducedDim(updated_object, "PCA") <- pca_matrix
      log <- progressLog(updated_object)
      log$PCAVariance <- pct_variance
      progressLog(updated_object) <- log
    }
    if (length(object@Clusters) > 0){
      # Retrieve old cluster list
      cluster_list <- object@Clusters
      
      # Reformat so it matches new cluster list
      distanceMatrix <- cluster_list$DistanceMatrix
      hClust <- cluster_list$Hclust
      putativeCluster <- cluster_list$PutativeClusters
      clusteringMatrix <- cluster_list$ClusteringMatrix # Need to change resolution labels
      colnames(clusteringMatrix)[1:ncol(clusteringMatrix) - 1] <- cluster_list$KeyStats$Height
      clusters <- cluster_list$Clusters
      nClusters <- cluster_list$NumberOfClusters
      optimalTreeHeight <- cluster_list$OptimalTreeHeight
      keyStats <- cluster_list$KeyStats
      
      # New list
      clusterAnalysis <- list(distanceMatrix = distanceMatrix,
                              hClust = hClust,
                              putativeCluster = putativeCluster,
                              clusteringMatrix = clusteringMatrix,
                              clusters = clusters,
                              nClusters = nClusters,
                              optimalTreeHeight = optimalTreeHeight,
                              keyStats = keyStats)
      clusterAnalysis(updated_object) <- clusterAnalysis
    }
    
    # Validate
    if (validateEMSet(updated_object)){
      # Calculate QC metrics
      object <- calculateQC(updated_object)
      print("Conversion complete! Returning object...")
      return(object)
    } else{
      stop("Converstion unsuccessful.")
    }
  }  
})
