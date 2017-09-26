# Validation checks
# Expression matrix check - COMPULSORY

#' CheckExpressionMatrix
#'
#' Check the expression matrix for issues.
#'
CheckExpressionMatrix <- function(object, errors){
  # Check expression matrix using plyr's empty function.
  if (any(nrow(object@ExpressionMatrix) == 0, ncol(object@ExpressionMatrix) == 0)) {
    msg <- "The expression matrix is empty. Please load an expression matrix."
    errors <- c(errors, msg)
  }

  # Check if there are any NA values in the expression matrix
  if (any(is.na(object@ExpressionMatrix))) {
    msg <- "There are NA values in the expression matrix. Please remove these values from the matrix."
    errors <- c(errors, msg)
  }

  # Check there are no duplicates in the colnames and rownames of the expression matrix
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

#' CheckControls
#' Check controls if they are present, for issues such as mismatches.
#'
CheckControls <- function(object, errors){
  # Control checks - verify presence of controls in the matrix.
  if (length(object@Controls) > 0) {
    ## Verify listed controls are present in the expression matrix
    check.controls <- unlist(object@Controls) %in% rownames(object@ExpressionMatrix)
    if (!any(check.controls)) {
      msg <- "Please make sure the gene identifiers you have listed in your control list matches those used in the expression matrix."
      errors <- c(errors, msg)
    }
  }
  return(errors)
}

#' CheckCellInformation
#'
#' Checks cell information supplied by the user.
#'
CheckCellInformation <- function(object, errors){
  ## Check if number of batch identifiers match the number of columns in the expression matrix
  if (nrow(object@CellInformation) != ncol(object@ExpressionMatrix)) {
    msg <- "Please check your batch information. Number of batch labels does not match the number of cells listed in the expression matrix."
  } else {
    ## Check if the cell identifiers are linked to cell barcodes in the expression matrix
    if (!all(colnames(object@ExpressionMatrix) %in% object@CellInformation[,1])) {
      msg <- "Please make sure cell identifiers in column 1 of the Cell Information data frame matches the cell identifiers used in the expression matrix."
    }
  }
  return(errors)
}

#' CheckGeneInformation
#'
#' Check Gene Annotation if supplied by the user.
#'
CheckGeneInformation <- function(object, errors){
  if( !all(rownames(object@ExpressionMatrix) %in% object@GeneInformation[,1]) ){
    msg <- "Please make sure gene identifiers in column 1 of the Gene Information data frame matches the gene identifiers used in the expression matrix."
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
  if (nrow(object@CellInformation) > 0) {
    errors <- CheckCellInformation(object, errors)
  }

  # Check if the annotation data matches rownames in the expression matrix
  if (nrow(object@GeneInformation) > 0) {
    errors <- CheckGeneInformation(object, errors)
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
#' @slot GeneInformation A data frame containing information a set of gene identifiers, such as gene symbols or ENSEMBL transcript identifiers. This data frame also holds information on controls and any information provided by the user.
#' @slot CellInformation A data frame containing each cell identifier, its associated batch/sample and additional information such as conditions.
#' @slot Controls A named list featuring gene identifiers to use as controls. These gene identifiers must match the identifiers used in the expression matrix.
#' @slot PCA Objects related to dimension reduction, such as a PCA matrixand a list of percentage variance values per principle component (PC). Populated by \code{\link{RunPCA}}.
#' @slot Clusters Objects related to clustering, including a distance matrix, a hclust object, cell identifiers and their associated cluster. Populated by \code{\link{FindOptimalClusters}}.
#' @slot Metrics A list of values generated by the \code{\link{GenerateMetrics}} function.
#' @slot Log A record of functions used on an \linkS4class{AEMSet}.
#' @export
setClass(
  "AEMSet",
  representation(
    ExpressionMatrix = "Matrix",
    GeneInformation = "data.frame",
    CellInformation = "data.frame",
    Controls = "list",
    PCA = "list",
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
    GeneInformation = data.frame(matrix(nr = 0, nc = 0)),
    CellInformation = data.frame(matrix(nr = 0, nc = 0)),
    Controls = list(),
    PCA = list(),
    Clusters = list(),
    Metrics = list(),
    Log = list()
  ),
  validity = CheckAEMSet
)

# More methods for this class
setMethod("show", signature("AEMSet"), function(object) {
  # Rethink this slot
  # Get number of genes and cells from the expression matrix dimensions
  print("ASCEND Object - AEMSet")
  n.genes <- nrow(object@ExpressionMatrix)
  n.cells <- ncol(object@ExpressionMatrix)
  print(sprintf("Expression Matrix: %i genes and %i cells", n.genes, n.cells))

  # If defined, show controls
  if (length(object@Controls) > 0){
    print("Controls:")
    print(object@Controls)
  }

  # If defined, show filtering log.
  if (!is.null(object@Log$FilteringLog)){
    print("Filtering Log:")
    print(object@Log$FilteringLog)
  }

  # If normalised, show normalisation
  if (!is.null(object@Log$NormalisationMethod)){
    print(sprintf("Normalisation method: %s", object@Log$NormalisationMethod))
  }

  # If clustered, show PCA
  if (!is.null(object@Log$PCA)){
    print(sprintf("PCA: %i dimensions", ncol(object@PCA$PCA)))
  }

  # If defined, show clustering information
  if (!is.null(object@Clusters$NumberOfClusters)){
    print(sprintf("Number of Clusters: %i", object@Clusters$NumberOfClusters))
  }
})

# Constructor function for AEMSet
#' NewAEMSet
#'
#' \code{\link{NewAEMSet}} generates a \linkS4class{AEMSet} object for use with the ASCEND package. This object contains an expression matrix, associated metadata, downstream analysis and a log documenting the actions taken to generate this object.
#' @param ExpressionMatrix An expression matrix in data.frame, dgCMatrix (sparse) or matrix format. Rows should represent a transcript and its counts, while columns should represent individual cells. This is usually the end point for Single Cell RNA-Seq pipelines such as Cell Ranger and DropSeq.
#' @param Controls A named list of controls, eg. mitochondrial genes, ribosomal genes and ERCC spike-ins. These genes must be specified using the identifier used in the expression matrix. This information is required for some functions.
#' @param GeneInformation A data frame containing gene identifiers used in the expression matrix. The first column should hold the cell identifiers you are using in the expression matrix. Other columns contain information about the genes, such as their corresponding ENSEMBL transcript identifiers, whether or not they are a control and any additional information supplied by the user. This is an optional field.
#' @param CellInformation A data frame containing cell identifiers (usually barcodes) and an integer representing which batch they belong to. This is an optional field, and it best used for experiments that contain data from multiple samples. This data frame can also hold additional information supplied by the user.
#' @return This function generates an object belonging to the \linkS4class{AEMSet}.
#' @seealso \linkS4class{AEMSet}
#' @export
NewAEMSet <-
  function(ExpressionMatrix = NULL,
           GeneInformation = NULL,
           CellInformation = NULL,
           Controls = list() ){
    # Check that we have the essential arguments - an expression matrix
    arg.check <- list(ExpressionMatrix = missing(ExpressionMatrix))
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
    if ( is.null(CellInformation) ){
      CellInformation <- data.frame(cell_barcode = colnames(ExpressionMatrix), batch = rep(1, ncol(ExpressionMatrix)))
    }

    # Automatically generate a gene annotation data frame if one hasn't
    # been supplied.
    if ( is.null(GeneInformation) ) {
      GeneInformation <- data.frame(gene_symbol = rownames(ExpressionMatrix))
    }

    # Convert the data.frame input into a sparse matrix
    if (is.data.frame(ExpressionMatrix) | is.matrix(ExpressionMatrix)) {
      sparse.matrix <- ConvertMatrix(ExpressionMatrix, format = "sparse.matrix")
    } else if (is(ExpressionMatrix, "sparseMatrix")) {
      sparse.matrix <- ExpressionMatrix
    } else {
      stop(
        "Please supply an expression matrix in one of the following formats: data.frame, dgCmatrix or matrix"
      )
    }

    # If the user has defined controls, add this information to the GeneInformation data frame.
    if( length(Controls) > 0){
      GeneInformation <- AddControlInfo(GeneInformation, Controls)
    }

    # Create a new AEMSet object.
    aem.set <-
      new(
        "AEMSet",
        ExpressionMatrix = sparse.matrix,
        GeneInformation = GeneInformation,
        CellInformation = CellInformation,
        Controls = Controls
      )

    # Generate metrics
    aem.set <- GenerateMetrics(aem.set)

    # Sync all the information
    aem.set <- SyncSlots(aem.set)

    # If controls
    if (length(Controls) > 0){
      aem.set@Log$Controls <- TRUE
    } else{
      aem.set@Log$Controls <- FALSE
    }

    # All clear, return the object
    return(aem.set)
  }

setGeneric(
  name = "SyncSlots",
  def = function(object) {
    standardGeneric("SyncSlots")
  }
)

#' SyncSlots
#'
#' Synchronises the data frames ExpressionMatrix, GeneInformation and CellInformation, to ensure everything is up to date.
#' This is important as the object will undergo a series of changes during the filtering and normalisation process.
#'
setMethod("SyncSlots", signature("AEMSet"), function(object){
  # Get data frames
  expression.matrix <- object@ExpressionMatrix
  gene.information <- object@GeneInformation
  cell.information <- object@CellInformation

  # Expression matrix serves as the basis for updating all the other information
  present.cells <- colnames(expression.matrix)
  present.genes <- rownames(expression.matrix)

  # Get data frames based on this info
  gene.information <- gene.information[gene.information[,1] %in% present.genes,]
  cell.information <- cell.information[cell.information[,1] %in% present.cells,]

  object@GeneInformation <- gene.information
  object@CellInformation <- cell.information

  return(object)
})


#' UpdateControls
#'
#' Replaces the control list in a \linkS4class{AEMSet} object with a new control list. This also recalculates the metrics associated with the \linkS4class{AEMSet} object.
#' For best results, define your controls before you attempt any filtering. You can also use this function to change an AEMSet into a control-less dataset.
#'
#' @param object An \linkS4class{AEMSet} object.
#' @param controls A named list containing the gene identifiers. These identifiers must match the identifiers used in the expression matrix. This list will replace any pre-existing controls.
#' @return This function returns an AEMSet object with the upated controls.
#' @export
setGeneric(
  name = "UpdateControls",
  def = function(object, controls) {
    standardGeneric("UpdateControls")
  }
)

setMethod("UpdateControls", signature("AEMSet"), function(object, controls) {
  errors <- character()

  # If user has supplied a list of controls
  if (length(controls) > 0) {
    check.controls <- unlist(controls) %in% rownames(object@ExpressionMatrix)
    if (!TRUE %in% check.controls) {
      msg <-
        "Please make sure the gene identifiers you have listed in your control list matches those used in the expression matrix."
      errors <- c(errors, msg)
    }
  }

  if (length(errors) > 0) {
    errors
  } else {
    object@Controls <- controls
    if (length(controls) > 0){
      GeneInformation <- object@GeneInformation
      GeneInformation <- AddControlInfo(GeneInformation, controls)
      object@GeneInformation <- GeneInformation
      object@Log$Controls <- TRUE
    } else{
      object@Log$Controls <- FALSE
    }
    object <- GenerateMetrics(object)
    return(object)
  }
})

#' ReplaceCellInfo
#'
#' Can be called by the user or by a filtering function. Updates Cell Information in a \linkS4class{AEMSet} object, replacing the old data frame with a new one.
#'
#' @param object An \linkS4class{AEMSet} object.
#' @param cell.info A data frame containing cell information. The first column must comprise of cell identifiers that match the cell identifiers used in the expression matrix. The second column, while optional - should contain batch information.
#' @include ASCEND_objects.R
#' @export
setGeneric(
  name = "ReplaceCellInfo",
  def = function(object, cell.info) {
    standardGeneric("ReplaceCellInfo")
  }
)

setMethod("ReplaceCellInfo", signature("AEMSet"), function(object, cell.info) {
  # Check replacement is valid.
  errors <- CheckCellInformation(cell.info)
  if (length(errors) > 0){
    stop("New Cell Information is invalid. Please check your data frame and try again.")
  } else{
    object@CellInformation <- cell.info
    return(object)
  }
})

#' ReplaceGeneInfo
#'
#' Can be called by the user or by a filtering function. Updates Gene Information in a \linkS4class{AEMSet} object, replacing the old data frame with a new one.
#'
#' @param object An \linkS4class{AEMSet} object.
#' @param gene.info A data frame containing cell information. The first column must comprise of gene identifiers that match the gene identifiers used in the expression matrix. The second column, while optional - should contain batch information.
#' @include ASCEND_objects.R
#' @export
setGeneric(
  name = "ReplaceGeneInfo",
  def = function(object, gene.info) {
    standardGeneric("ReplaceGeneInfo")
  }
)

setMethod("ReplaceGeneInfo", signature("AEMSet"), function(object, gene.info) {
  # Check replacement is valid.
  errors <- CheckGeneInformation(gene.info)
  if (length(errors) > 0){
    stop("New Gene Information is invalid. Please check your data frame and try again.")
  } else{
    object@GeneInformation <- gene.info
    return(object)
  }
})

#' ReplaceExpressionMatrix
#'
#' Replace the expression matrix in a \linkS4class{AEMSet} with a new expression matrix and re-calculate its metrics.
#'
#' @param object The \linkS4class{AEMSet} you would like to update.
#' @param expression.matrix Expression matrix in matrix form
#' @include ASCEND_objects.R
#' @export
setGeneric(
  name = "ReplaceExpressionMatrix",
  def = function(object, expression.matrix) {
    standardGeneric("ReplaceExpressionMatrix")
  }
)

setMethod("ReplaceExpressionMatrix", signature("AEMSet"), function(object, expression.matrix) {
  # Replace the matrix
  object@ExpressionMatrix <- ConvertMatrix(expression.matrix, format="sparse.matrix")

  # Check that it's okay before replacing.
  errors <- c()
  errors <- CheckExpressionMatrix(object, errors)
  if (length(errors) > 0){
    stop("Please check your expression matrix.")
  }else{
    object <- GenerateMetrics(object)
    return(object)
  }
})
