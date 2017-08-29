# Validation checks
# Expression matrix check - COMPULSORY
CheckExpressionMatrix <- function(object, errors){
  # Check expression matrix using plyr's empty function.
  if (plyr::empty(object@ExpressionMatrix)) {
    msg <-
      "The expression matrix is empty. Please load an expression matrix."
    errors <- c(errors, msg)
  }

  # Check if there are any NA values in the expression matrix
  if (any(is.na(object@ExpressionMatrix))) {
    msg <-
      "There are NA values in the expression matrix. Please remove these values from the matrix."
    errors <- c(errors, msg)
  }

  # Check there are no duplicates in the colnames and rownames of the
  # expression matrix
  if (any(duplicated(colnames(object@ExpressionMatrix)))) {
    msg <-
      "Duplicate column names are present. Please check the expression matrix."
    errors <- c(errors, msg)
  }

  if (any(duplicated(rownames(object@ExpressionMatrix)))) {
    msg <-
      "Duplicate row names are present. Please check the expression matrix."
    errors <- c(errors, msg)
  }

  return(errors)
}

# Control check - COMPULSORY
CheckControls <- function(object, errors){
  # Control checks - verify presence of controls in the matrix.
  if (length(object@Controls) > 0) {
    ## Verify listed controls are present in the expression matrix
    check.controls <- unlist(object@Controls) %in% rownames(object@ExpressionMatrix)
    if (!any(check.controls)) {
      msg <- "Please make sure the gene identifiers you have listed in your control list matches those used in the expression matrix."
      errors <- c(errors, msg)
    }
  } else {
    msg <- "Please make sure that you have defined a control."
    errors <- c(errors, msg)
  }
  return(errors)
}

# Check Batch Information - OPTIONAL
CheckBatchInfo <- function(object, errors){
  ## Check if number of batch identifiers match the number of columns in the expression matrix
  if (length(object@BatchInformation) != length(colnames(object@ExpressionMatrix))) {
    msg <-
      "Please check your batch information. Number of batch labels does not match the number of cells listed in the expression matrix."
  } else {
    ## Check if the batch identifiers are linked to cell barcodes in the expression matrix
    if (!all(colnames(object@ExpressionMatrix) %in% names(object@BatchInformation))) {
      msg <- "Please make sure the list supplied to BatchInformation slot is formatted as follows: Cell Identifier, Batch Number"
    }
  }
  return(errors)
}

# Check Gene Annotation - OPTIONAL
CheckGeneAnnotation <- function(object, errors){
  check.annotations <- lapply(object@GeneAnnotation, function(x) CheckData(x, rownames(object@ExpressionMatrix)))
  if (!TRUE %in% check.annotations) {
    msg <- "Please check your annotation list. All genes in the expression matrix must be found in the gene annotation."
    errors <- c(errors, msg)
  }
  return(errors)
}

#' CheckAEMSet
#'
#' Validation function for the S4 class \linkS4class{AEMSet}.
#'
#' @details The checks are as follows:
#' \enumerate{
#'     \item Check if the expression matrix is empty.
#'     \item Check if there are any NA values in the expression matrix.
#'     \item Check if there are any duplicated rownames or colnames in the expression matrix.
#'     \item Check gene identifiers are present in the supplied control list.
#'     \item Check if gene identifiers listed in the supplied controls are present in the expression matrix.
#'     \item If supplied, verify batch information matches colnames in the expression matrix.
#'     \item If supplied, verify gene identifiers match rownames in the expression matrix.
#' }
#' @export
#'
CheckAEMSet <- function(object) {
  # Error message container
  # This variable is only populated if an error is encountered.
  errors <- character()

  # Compulsory checks
  errors <- CheckExpressionMatrix(object, errors)
  errors <- CheckControls(object, errors)

  # Optional checks
  # Batch label check This step is optional
  if (length(object@BatchInformation) > 0) {
    errors <- CheckBatchInfo(object, errors)
  }

  # Check if the annotation data matches rownames in the expression matrix
  if (!plyr::empty(object@GeneAnnotation)) {
    errors <- CheckGeneAnnotation(object, errors)
  }

  # If all is well, object is valid. If not, it's invalid.
  if (length(errors) == 0) {
    TRUE
  } else{
    errors
  }
}

# ASCEND Expression and Metadata Set Class Definition
#' ASCEND Expression and Metadata Set (AEMSet)
#'
#' An S4 class to contain data in a format ASCEND can work with for analysis.
#' @slot ExpressionMatrix Transcript counts stored as a sparse matrix, where rows are transcript/gene identifiers and columns are invididual cells.
#' @slot GeneAnnotation A data frame containing information on genes and their corresponding identifiers, such as gene names and ENSEMBL transcript identifiers.
#' @slot BatchInformation A named list containing each cell identifier and its associated batch/sample.
#' @slot Controls A named list featuring gene identifiers to use as controls. These gene identifiers must match the identifiers used in the expression matrix.
#' @slot Conditions A named list containing a list of conditions and what category each cell identifier falls into.
#' @slot PCA Objects related to dimension reduction, such as a PCA matrixand a list of percentage variance values per principle component (PC). Populated by \code{\link{RunPCA}}.
#' @slot Clusters Objects related to clustering, including a distance matrix, a hclust object, cell identifiers and their associated cluster. Populated by \code{\link{FindOptimalClusters}}.
#' @slot Metrics A list of values generated by the \code{\link{GenerateMetrics}} function.
#' @slot Log A record of functions used on an \linkS4class{AEMSet}.
#' @export
setClass(
  "AEMSet",
  representation(
    ExpressionMatrix = "Matrix",
    GeneAnnotation = "data.frame",
    BatchInformation = "list",
    Controls = "list",
    PCA = "list",
    Conditions = "list",
    Clusters = "list",
    Metrics = "list",
    Log = "list"
  ),
  prototype(
    ExpressionMatrix = Matrix::Matrix(
      0,
      nrow = 0,
      ncol = 0,
      sparse = TRUE
    ),
    GeneAnnotation = data.frame(matrix(nr = 0, nc = 0)),
    BatchInformation = list(),
    Controls = list(),
    Conditions = list(),
    PCA = list(),
    Clusters = list(),
    Metrics = list(),
    Log = list()
  ),
  validity = CheckAEMSet
)

# More methods for this class
setMethod("show", signature("AEMSet"), function(object) {
  # Get number of genes and cells from the expression matrix dimensions
  print("ASCEND Object - AEMSet")
  n.genes <- dim(object@ExpressionMatrix)[1]
  n.cells <- dim(object@ExpressionMatrix)[2]
  expression.matrix.str <-
    sprintf("Expression Matrix: %i cells and %i genes",
            n.cells, n.genes)
  print(expression.matrix.str)
  batches <- table(as.vector(unlist(object@BatchInformation)))
  print(sprintf("Batch Information:"))
  print(batches)
  print("Controls:")
  for (control.name in names(object@Controls)) {
    n.control.genes <- length(unlist(object@Controls[control.name]))
    control.str <-
      sprintf("%s: %i genes", control.name, n.control.genes)
    print(control.str)
  }

  if (length(object@Log) > 0) {
    print(object@Log)
  }

  if (length(object@Clusters) > 0) {
    if (length(object@Clusters$ClusterList) > 0) {
      clusters <- object@Clusters$ClusterList
      ncluster <- length(unique(clusters))
      cluster.str <- sprintf("Number of Clusters: %i", ncluster)
      print(cluster.str)
      print("Cluster sizes:")
      count.table <- table(clusters)
      print(count.table)
    }
  }
})

# Constructor function for AEMSet
#' NewAEMSet
#'
#' \code{\link{NewAEMSet}} generates a \linkS4class{AEMSet} object for use with the ASCEND package. This object contains an expression matrix, associated metadata, downstream analysis and a log documenting the actions taken to generate this object.
#' @param ExpressionMatrix An expression matrix in data.frame, dgCMatrix (sparse) or matrix format. Rows should represent a transcript and its counts, while columns should represent individual cells. This is usually the end point for Single Cell RNA-Seq pipelines such as Cell Ranger and DropSeq.
#' @param Controls A named list of controls, eg. mitochondrial genes, ribosomal genes and ERCC spike-ins. These genes must be specified using the identifier used in the expression matrix. This is required.
#' @param GeneAnnotation A data frame containing gene identifiers used in the expression matrix. Other columns can be used to represent corresponding identifiers in other formats (such as ENSEMBL transcript IDs). This is an optional field, and is best used if you need to convert between identifiers.
#' @param BatchInformation A named list of cell identifiers (usually barcodes) and an integer representing which batch they belong to. This is an optional field, and it best used for experiments that contain data from multiple samples.
#' @return This function generates an object belonging to the \linkS4class{AEMSet}.
#' @seealso \linkS4class{AEMSet}
#' @export
NewAEMSet <-
  function(ExpressionMatrix = NULL,
           GeneAnnotation = NULL,
           BatchInformation = NULL,
           Controls = list(Mt = list(), Rb = list())) {
    # Check that we have the essential arguments
    arg.check <-
      list(ExpressionMatrix = missing(ExpressionMatrix),
           Controls = missing(Controls))
    if (any(arg.check == TRUE)) {
      missing.args <- names(which(arg.check == TRUE))
      msg <-
        sprintf(
          "Please check that you have supplied the following arguments: %s\n",
          as.character(unlist(missing.args))
        )
      stop(msg)
    }

    # Automatically generate batch labels if it's not defined.
    if (is.null(BatchInformation)) {
      BatchInformation <-
        setNames(rep(1, ncol(ExpressionMatrix)), as.list(colnames(ExpressionMatrix)))
      BatchInformation <- as.list(BatchInformation)
    }

    # Automatically generate a gene annotation data frame if one hasn't
    # been supplied.
    if (is.null(GeneAnnotation)) {
      GeneAnnotation <- data.frame(gene_id = rownames(ExpressionMatrix))
    }

    # Convert the data.frame input into a sparse matrix
    if (is.data.frame(ExpressionMatrix)) {
      sparse.matrix <- LoadSparseMatrix(ExpressionMatrix)
    } else if (is.matrix(ExpressionMatrix)) {
      sparse.matrix <- Matrix::Matrix(ExpressionMatrix, sparse = TRUE)
    } else if (is(ExpressionMatrix, "sparseMatrix")) {
      sparse.matrix <- ExpressionMatrix
    } else {
      stop(
        "Please supply an expression matrix in one of the following formats: data.frame, dgCmatrix or matrix"
      )
    }

    # Create a new AEMSet object.
    aem.set <-
      new(
        "AEMSet",
        ExpressionMatrix = sparse.matrix,
        GeneAnnotation = GeneAnnotation,
        BatchInformation = BatchInformation,
        Controls = Controls
      )
    aem.set <- GenerateMetrics(aem.set)

    # Filter out identifiers that are not in the expression matrix
    filtered.control.list <- sapply(aem.set@Controls, function(x) x[x %in% rownames(aem.set@ExpressionMatrix)], USE.NAMES = TRUE, simplify = FALSE)
    aem.set@Controls <- filtered.control.list

    # Store column ID so we can keep track
    gene.list <- rownames(ExpressionMatrix)
    annotation.col <- sapply(names(GeneAnnotation), function(x){if(all(gene.list %in% GeneAnnotation[[x]])){return(TRUE)} else{return(FALSE)}}, USE.NAMES = TRUE)
    aem.set@Log$GeneAnnotation <- annotation.col

    return(aem.set)
  }


