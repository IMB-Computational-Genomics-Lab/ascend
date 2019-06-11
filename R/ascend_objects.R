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
#' 
#' @examples
#' # Load example EMSet
#' em_set <- ascend::raw_set
#' 
#' # Validate
#' validateEMSet(em_set)
#' 
#' @return Validated \linkS4class{EMSet}.
#' @export
validateEMSet <- function(object){
  # Store errors in an object. If length is greater than 0, indicates multiple
  # errors are present.
  errors <- c()
  
  # LET'S VALIDATE THIS OBJECT!
  # Get gene identifiers
  gene_list <- rownames(object)
  
  # Get cell identifiers
  cell_list <- colnames(object)
  
  # Now we check all the elements are the same
  ## Check cell identifiers for EMSet match those in colInfo
  ## Check first column of colInfo matches others
  colInfo <- colInfo(object)
  if (!(identical(cell_list, rownames(colInfo)))){
    errors <- c(errors, "colInfo rownames do not match EMSet colnames.")
  } else{
    if (!(identical(as.vector(colInfo[,1]), cell_list))){
      errors <- c(errors, "First column of colInfo does not match EMSet colnames and colInfo rownames.")
    }
  }
  
  rowInfo <- rowInfo(object)
  ## Check gene names for EMSet match those in rowInfo
  if (!(identical(gene_list, rownames(rowInfo)))){
    errors <- c(errors, "rowInfo rownames do not match EMSet rownames.")
  } else{
    if (!(identical(as.vector(rowInfo[,1]), gene_list))){
      errors <- c(errors, "First column of rowInfo does not match EMSet rownames and rowInfo rownames.")
    }
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
#' 
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
#' @examples 
#' # Load example EMSet from package
#' em_set <- ascend::raw_set
#' 
#' class(em_set)
#' 
#' @name EMSet-class
#' @rdname EMSet-class
#'
#' @import methods
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
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

#' @param object \linkS4class{EMSet}
#' @importFrom SummarizedExperiment assayNames rowData colData
#' @importFrom S4Vectors metadata
#' @importFrom BiocGenerics rownames colnames
#' @importFrom SingleCellExperiment reducedDimNames spikeNames
#' @rdname EMSet
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

loadCounts <- function(x){
  if (is.data.frame(x)){
    cell_barcodes <- colnames(x)
    x <- as.matrix(x)
    x <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = x))
    colnames(x) <- cell_barcodes
  }
  
  # If it is a matrix
  else if (is.matrix(x)){
    cell_barcodes <- colnames(x)
    x <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = x))
    colnames(x) <- cell_barcodes
  }
  
  # If it is a sparse matrix
  else if (is(x, "sparseMatrix")){
    cell_barcodes <- colnames(x)
    x <- SingleCellExperiment::SingleCellExperiment(assays = list(counts= x))
    colnames(x) <- cell_barcodes
  }
  
  # If it is a list
  else if (is.list(x)){
    if ("counts" %in% names(x)){
      cell_barcodes <- colnames(x$counts)
      x <- SingleCellExperiment::SingleCellExperiment(x)
      colnames(x) <- cell_barcodes
    } else{
      stop("Please supply a list with an object labelled counts.")
    }
  }
  
  # If it is a simple list
  else if (is(x, "SimpleList")){
    if ("counts" %in% names(x)){
      cell_barcodes <- colnames(x$counts)
      x <- SingleCellExperiment::SingleCellExperiment(x)
      colnames(x) <- cell_barcodes
    } else{
      stop("Please supply a list with an object labelled counts.")
    }
  }
  
  # Check that we have created a single cell experiment
  if (is(x, "SingleCellExperiment")){
    if (is.null(colnames(x))){
      cell_barcodes <- colnames(counts(x))
      colnames(x) <- cell_barcodes
    }
    return(x)
  } else{
    stop("Please supply counts in an accepted format.")
  }
}

#' EMSet
#' 
#' \code{\link{EMSet}} generates a \linkS4class{EMSet} object for use with 
#' the `ascend` package. This object contains an expression matrix, 
#' associated metadata, downstream analysis and a log documenting the actions 
#' used to shape the data in this object.
#' 
#' @param x An object containing count data. The counts can be stored in the 
#' following formats: SingleCellExperiment, "counts" in a list or SimpleList,
#' dense matrix or sparse matrix.
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
#' # EMSet from a list
#' em_set <- EMSet(list(counts = count_matrix))
#' 
#' # EMSet from a matrix
#' em_set <- EMSet(count_matrix)
#' 
#' @return An \linkS4class{EMSet} object.
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData rowData
#' @export
EMSet <- function(x,
                     colInfo = NULL,
                     colData = NULL,
                     rowInfo = NULL,
                     rowData = NULL,
                     controls = NULL){

      # Load count data into a SingleCellExperiment object if possible
      # Expect errors to be picked up by SingleCellExperiment methods
      counts <- loadCounts(x)
      
      # Get colData supplied from the SingleCellExperiment
      # OR check one user has supplied
      # OR create a new one
      if ((nrow(SummarizedExperiment::colData(counts)) > 0) & (ncol(SummarizedExperiment::colData(counts)) > 0)){
        # Set first column as cell_barcode
        # Validate cell barcode
        if ("cell_barcode" %in% colnames(colData(counts))){
          cell_barcode_idx <- which("cell_barcode" == colnames(colData(counts)))
          other_idx <- which("cell_barcode" != colnames(colData(counts)))
          colData <- colData(counts)
          colData <- colData[, c(cell_barcode_idx, other_idx)]
        } else{
          colData <- colData(counts)
          cell_barcode <- colnames(counts)
          colData$cell_barcode <- cell_barcode
          rownames(colData) <- cell_barcode
          colData <- colData[, c(ncol(colData), which(colnames(colData) != "cell_barcode"))]
        }
      } else{
        # Create colData if user hasn't supplied one
        if (missing(colData)){
          if (missing(colInfo)){
            colData <- DataFrame(cell_barcode = colnames(counts), row.names = colnames(counts))      
          } else{
            colData <- DataFrame(cell_barcode = colInfo[, 1], row.names = colInfo[, 1]) 
          }
        } else{
          # Check if the supplied colData matches cell identifiers in the count data
          if (!(all(colData[, 1] %in% colnames(counts)))){
            stop("Supplied colData DataFrame does not match supplied counts.")
          }
          
          # If supplied, check if colInfo matches colData
          if (!(missing(colInfo))){
            if (!(all(colData[, 1] %in% colInfo[, 1]))){
              stop("First column of colData does not match first column of colInfo.")
            }     
          }
        }
      }
      colData(counts) <- colData
      
      if (ncol(rowData(counts)) > 0){
        rowData <- rowData(counts)
        
        # Check gene ids are present
        gene_ids <- rownames(counts)
        gene_col_query <- apply(rowData, 2, function(y) all(y == gene_ids))
        
        # Retrieve index if present, make one if not
        if (any(gene_col_query)){
          gene_id_name <- names(which(gene_col_query))
        } else{
          rowData$gene_id <- gene_ids
          gene_id_name <- "gene_id"
          gene_idx <- ncol(rowData)
        }
        
        # Move column to the front
        if (which(colnames(rowData) == gene_id_name) != 1){
          other_idx <- 1:ncol(rowData)
          other_idx <- other_idx[which(colnames(rowData) != gene_id_name)]
          rowData <- rowData[, c(gene_idx, other_idx)]
        }
        
      } else{
        if (missing(rowData)){
          if (missing(rowInfo)){
            # Ensure we use the same header as other datasets
            gene_id_name <- "gene_id"
            rowData <- DataFrame(gene_id = rownames(counts), row.names = rownames(counts))      
          } else{
            gene_id_name <- colnames(rowInfo)[1]
            gene_ids <- list()
            gene_ids[[gene_id_name]] <- rowInfo[, 1]
            rowData <- DataFrame(gene_ids, row.names = rownames(counts))
          }
          
        } else{
          # Check the columns match
          gene_id_name <- colnames(rowData)[1]
          if (!(all(rowData[, gene_id_name] %in% rownames(counts)))){
            stop("Gene identifiers in rowData do not match identifiers used in supplied data.")
          }
          if (!(missing(rowInfo))){
            if (!(all(rowData[,1] %in% rowInfo[, 1]))){
              stop("First column of rowData does not match first column of rowInfo.")
            }
          }
        }
      }
      
      rowData(counts) <- rowData
      
      # Now check colInfo
      if (missing(colInfo)){
        colInfo <- DataFrame(cell_barcode = colnames(counts), 
                                        batch = rep(1, ncol(counts)), 
                                        row.names = colnames(counts)) 
      } else{
        if("cell_barcode" != colnames(colInfo)[1]){
          stop("Please specify the name of the first column in colInfo as 'cell_barcode'")
        }
        if (!("batch" %in% colnames(colInfo))){
          colInfo$batch <- 1
        }
      }
      
      # Check rowInfo
      if (missing(rowInfo)){
        gene_list <- list()
        gene_list[[gene_id_name]] <- rownames(counts)
        rowInfo <- DataFrame(gene_list, row.names = rownames(counts))
      }
      
      # Make sure rownames are column 1 in colInfo, rowInfo, colData, rowData
      colInfo <- DataFrame(colInfo, row.names = as.vector(colInfo[, 1]))
      rowInfo <- DataFrame(rowInfo, row.names = as.vector(rowInfo[, 1]))
      
      # Now that everything is in order, create an EMSet
      object <- new("EMSet",
                    counts,
                    colInfo = colInfo,
                    rowInfo = rowInfo)
      
      # Now add controls if they are present
      if (!(is.null(controls))){
        controls(object) <- controls
      }
      
  # Calculate QC metrics
  object <- calculateQC(object)
  
  # Return to user
  return (object)
}

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
    updated_object <- EMSet(list(counts = as.matrix(counts)),
                               colInfo = S4Vectors::DataFrame(cell_info),
                               rowInfo = S4Vectors::DataFrame(gene_info))
    
    SummarizedExperiment::mcols(updated_object) <- S4Vectors::DataFrame(gene_info)
    
    if(length(object@Controls) > 0){
      updated_object <- controls(updated_object, controls = object@Controls)
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
