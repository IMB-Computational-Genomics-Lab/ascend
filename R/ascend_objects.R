# Validation checks Expression matrix check - COMPULSORY

#' CheckExpressionMatrix
#'
#' Check the expression matrix for issues.
#'
CheckExpressionMatrix <- function(object, errors) {
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
        msg <- "Duplicate column names are present. Please check the expression matrix."
        errors <- c(errors, msg)
    }

    if (any(duplicated(rownames(object@ExpressionMatrix)))) {
        msg <- "Duplicate row names are present. Please check the expression matrix."
        errors <- c(errors, msg)
    }

    return(errors)
}

#' CheckControls
#' Check controls if they are present, for issues such as mismatches.
#'
CheckControls <- function(object, errors) {
    # Control checks - verify presence of controls in the matrix.
    if (length(object@Controls) > 0) {
        ## Verify listed controls are present in the expression matrix
        check.controls <- unlist(object@Controls) %in% rownames(object@ExpressionMatrix)
        if (!any(check.controls)) {
            msg <- "Please make sure the gene identifiers you have listed in your
      control list matches those used in the expression matrix."
            errors <- c(errors, msg)
        }
    }
    return(errors)
}

#' CheckCellInformation
#'
#' Checks cell information supplied by the user.
#'
CheckCellInformation <- function(object, errors) {
    ## Check if number of batch identifiers match the number of columns in the expression matrix
    if (nrow(object@CellInformation) != ncol(object@ExpressionMatrix)) {
        msg <- "Please check your batch information. Number of batch labels does
    not match the number of cells listed in the expression matrix."
    } else {
        ## Check if the cell identifiers are linked to cell barcodes in the expression matrix
        if (!all(colnames(object@ExpressionMatrix) %in% object@CellInformation[, 1])) {
            msg <- "Please make sure cell identifiers in column 1 of the Cell
      Information data frame matches the cell identifiers used in the expression
      matrix."
        }
    }
    return(errors)
}

#' CheckGeneInformation
#'
#' Check Gene Annotation if supplied by the user.
#'
CheckGeneInformation <- function(object, errors) {
    if (!all(rownames(object@ExpressionMatrix) %in% object@GeneInformation[, 1])) {
        msg <- "Please make sure gene identifiers in column 1 of the Gene Information
    data frame matches the gene identifiers used in the expression matrix."
    }
    return(errors)
}

#' CheckEMSet
#'
#' Validation function for the S4 class \linkS4class{EMSet}.
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
CheckEMSet <- function(object) {
    # Error message container This variable is only populated if an error is encountered.
    errors <- character()

    # Compulsory checks
    errors <- CheckExpressionMatrix(object, errors)
    errors <- CheckControls(object, errors)

    # Optional checks Batch label check This step is optional
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
    } else {
        errors
    }
}

# Expression and Metadata Set Class Definition
#' Expression and Metadata Set (EMSet)
#'
#' An S4 class to contain data in a format ascend can work with for analysis.
#' @slot ExpressionMatrix Transcript counts stored as a sparse matrix, where rows are transcript/gene identifiers and columns are invididual cells.
#' @slot GeneInformation A data frame containing information a set of gene identifiers, such as gene symbols or ENSEMBL transcript identifiers. This data frame also holds information on controls and any information provided by the user.
#' @slot CellInformation A data frame containing each cell identifier, its associated batch/sample and additional information such as conditions.
#' @slot Controls A named list featuring gene identifiers to use as controls. These gene identifiers must match the identifiers used in the expression matrix.
#' @slot PCA Objects related to dimension reduction, such as a PCA matrixand a list of percentage variance values per principle component (PC). Populated by \code{\link{RunPCA}}.
#' @slot Clusters Objects related to clustering, including a distance matrix, a hclust object, cell identifiers and their associated cluster. Populated by \code{\link{RunCORE}}.
#' @slot Metrics A list of values generated by the \code{\link{GenerateMetrics}} function.
#' @slot Log A record of functions used on an \linkS4class{EMSet}.
#' @export
setClass("EMSet", representation(ExpressionMatrix = "Matrix", GeneInformation = "data.frame", CellInformation = "data.frame", Controls = "list", PCA = "list",
    Clusters = "list", Metrics = "list", Log = "list"), prototype(ExpressionMatrix = Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE), GeneInformation = data.frame(matrix(nr = 0,
    nc = 0)), CellInformation = data.frame(matrix(nr = 0, nc = 0)), Controls = list(), PCA = list(), Clusters = list(), Metrics = list(), Log = list()), validity = CheckEMSet)

# More methods for this class
setMethod("show", signature("EMSet"), function(object) {
    # Rethink this slot Get number of genes and cells from the expression matrix dimensions
    print("ascend Object - EMSet")
    n.genes <- nrow(object@ExpressionMatrix)
    n.cells <- ncol(object@ExpressionMatrix)
    print(sprintf("Expression Matrix: %i genes and %i cells", n.genes, n.cells))

    # If defined, show controls
    if (length(object@Controls) > 0) {
        print("Controls:")
        print(object@Controls)
    }

    # If defined, show filtering log.
    if (!is.null(object@Log$FilteringLog)) {
        print("Filtering Log:")
        print(object@Log$FilteringLog)
    }

    # If normalised, show normalisation
    if (!is.null(object@Log$NormalisationMethod)) {
        print(sprintf("Normalisation method: %s", object@Log$NormalisationMethod))
    }

    # If clustered, show PCA
    if (!is.null(object@Log$PCA)) {
        print(sprintf("PCA: %i dimensions", ncol(object@PCA$PCA)))
    }

    # If defined, show clustering information
    if (!is.null(object@Clusters$NumberOfClusters)) {
        print(sprintf("Optimal tree-cut height: %f", object@Clusters$OptimalTreeHeight))
        print(sprintf("Number of Clusters: %i", object@Clusters$NumberOfClusters))
    }

})

# Constructor function for EMSet
#' NewEMSet
#'
#' \code{\link{NewEMSet}} generates a \linkS4class{EMSet} object for use with the ascend package. This object contains an expression matrix, associated metadata, downstream analysis and a log documenting the actions taken to generate this object.
#' @param ExpressionMatrix An expression matrix in data.frame, dgCMatrix (sparse) or matrix format. Rows should represent a transcript and its counts, while columns should represent individual cells. This is usually the end point for Single Cell RNA-Seq pipelines such as Cell Ranger and DropSeq.
#' @param Controls A named list of controls, eg. mitochondrial genes, ribosomal genes and ERCC spike-ins. These genes must be specified using the identifier used in the expression matrix. This information is required for some functions.
#' @param GeneInformation A data frame containing gene identifiers used in the expression matrix. The first column should hold the cell identifiers you are using in the expression matrix. Other columns contain information about the genes, such as their corresponding ENSEMBL transcript identifiers, whether or not they are a control and any additional information supplied by the user. This is an optional field.
#' @param CellInformation A data frame containing cell identifiers (usually barcodes) and an integer representing which batch they belong to. This is an optional field, and it best used for experiments that contain data from multiple samples. This data frame can also hold additional information supplied by the user.
#' @return This function generates an object belonging to the \linkS4class{EMSet}.
#' @seealso \linkS4class{EMSet}
#' @export
NewEMSet <- function(ExpressionMatrix = NULL, GeneInformation = NULL, CellInformation = NULL, Controls = list()) {
    # Check that we have the essential arguments - an expression matrix
    arg.check <- list(ExpressionMatrix = missing(ExpressionMatrix))
    if (any(arg.check == TRUE)) {
        missing.args <- names(which(arg.check == TRUE))
        msg <- sprintf("Please check that you have supplied the following arguments: %s\n", as.character(unlist(missing.args)))
        stop(msg)
    }

    # Automatically generate batch labels if it's not defined.
    if (is.null(CellInformation)) {
        CellInformation <- data.frame(cell_barcode = colnames(ExpressionMatrix), batch = rep(1, ncol(ExpressionMatrix)))
    }

    # Automatically generate a gene annotation data frame if one hasn't been supplied.
    if (is.null(GeneInformation)) {
        GeneInformation <- data.frame(gene_symbol = rownames(ExpressionMatrix))
    }

    # Convert the data.frame input into a sparse matrix
    if (is.data.frame(ExpressionMatrix) | is.matrix(ExpressionMatrix)) {
        sparse.matrix <- ConvertMatrix(ExpressionMatrix, format = "sparseMatrix")
    } else if (is(ExpressionMatrix, "sparseMatrix")) {
        sparse.matrix <- ExpressionMatrix
    } else {
        stop("Please supply an expression matrix in one of the following formats: data.frame, dgCmatrix or matrix")
    }

    # If the user has defined controls, add this information to the GeneInformation data frame.
    if (length(Controls) > 0) {
        GeneInformation <- AddControlInfo(GeneInformation, Controls)
    }

    # Create a new EMSet object.
    em.set <- new("EMSet", ExpressionMatrix = sparse.matrix, GeneInformation = GeneInformation, CellInformation = CellInformation, Controls = Controls)

    # Generate metrics
    em.set <- GenerateMetrics(em.set)

    # Sync all the information
    em.set <- SyncSlots(em.set)

    # If controls
    if (length(Controls) > 0) {
        em.set@Log$Controls <- TRUE
    } else {
        em.set@Log$Controls <- FALSE
    }

    # All clear, return the object
    return(em.set)
}

setGeneric(name = "SyncSlots", def = function(object) {
    standardGeneric("SyncSlots")
})

#' SyncSlots
#'
#' Synchronises the data frames ExpressionMatrix, GeneInformation and CellInformation, to ensure everything is up to date.
#' This is important as the object will undergo a series of changes during the filtering and normalisation process.
#'
setMethod("SyncSlots", signature("EMSet"), function(object) {
    # Get data frames
    expression.matrix <- object@ExpressionMatrix
    gene.information <- object@GeneInformation
    cell.information <- object@CellInformation
    
    gene.cols <- colnames(gene.information)
    cell.cols <- colnames(cell.information)
    
    # Expression matrix serves as the basis for updating all the other information
    present.cells <- colnames(expression.matrix)
    present.genes <- rownames(expression.matrix)

    # Get data frames based on this info
    gene.information <- gene.information[gene.information[, 1] %in% present.genes, ]
    cell.information <- cell.information[cell.information[, 1] %in% present.cells, ]
    
    # Check in case it's not a dataframe
    if (!is.data.frame(gene.information)){
      gene.information <- data.frame(gene.information)
      colnames(gene.information) <- gene.cols
    }
    
    if (!is.data.frame(cell.information)){
      cell.information <- data.frame(cell.information)
      colnames(cell.information) <- cell.cols
    }
    
    # Put back into the slots
    object@GeneInformation <- gene.information
    object@CellInformation <- cell.information

    return(object)
})

#' UpdateControls
#'
#' Replaces the control list in a \linkS4class{EMSet} object with a new control list. This also recalculates the metrics associated with the \linkS4class{EMSet} object.
#' For best results, define your controls before you attempt any filtering. You can also use this function to change an EMSet into a control-less dataset.
#'
#' @param object An \linkS4class{EMSet} object.
#' @param controls A named list containing the gene identifiers. These identifiers must match the identifiers used in the expression matrix. This list will replace any pre-existing controls.
#' @return This function returns an EMSet object with the upated controls.
#' @export
setGeneric(name = "UpdateControls", def = function(object, controls) {
    standardGeneric("UpdateControls")
})

setMethod("UpdateControls", signature("EMSet"), function(object, controls) {
    errors <- character()

    # If user has supplied a list of controls
    if (length(controls) > 0) {
        check.controls <- unlist(controls) %in% rownames(object@ExpressionMatrix)
        if (!TRUE %in% check.controls) {
            msg <- "Please make sure the gene identifiers you have listed in your control list matches those used in the expression matrix."
            errors <- c(errors, msg)
        }
    }

    if (length(errors) > 0) {
        errors
    } else {
        object@Controls <- controls
        if (length(controls) > 0) {
            GeneInformation <- object@GeneInformation
            GeneInformation <- AddControlInfo(GeneInformation, controls)
            object@GeneInformation <- GeneInformation
            object@Log$Controls <- TRUE
        } else {
            object@Log$Controls <- FALSE
        }
        object <- GenerateMetrics(object)
        return(object)
    }
})

#' ReplaceCellInfo
#'
#' Can be called by the user or by a filtering function. Updates Cell Information in a \linkS4class{EMSet} object, replacing the old data frame with a new one.
#'
#' @param object An \linkS4class{EMSet} object.
#' @param cell.info A data frame containing cell information. The first column must comprise of cell identifiers that match the cell identifiers used in the expression matrix. The second column, while optional - should contain batch information.
#' @include ascend_objects.R
#' @export
setGeneric(name = "ReplaceCellInfo", def = function(object, cell.info) {
    standardGeneric("ReplaceCellInfo")
})

setMethod("ReplaceCellInfo", signature("EMSet"), function(object, cell.info) {
    # Check replacement is valid.
    replacement <- object
    replacement@CellInformation <- cell.info

    errors <- c()
    errors <- CheckCellInformation(replacement, errors)

    if (length(errors) > 0) {
        stop("New Cell Information is invalid. Please check your data frame and try again.")
    } else {
        return(replacement)
    }
})

#' ReplaceGeneInfo
#'
#' Can be called by the user or by a filtering function. Updates Gene Information in a \linkS4class{EMSet} object, replacing the old data frame with a new one.
#'
#' @param object An \linkS4class{EMSet} object.
#' @param gene.info A data frame containing cell information. The first column must comprise of gene identifiers that match the gene identifiers used in the expression matrix. The second column, while optional - should contain batch information.
#' @include ascend_objects.R
#' @export
setGeneric(name = "ReplaceGeneInfo", def = function(object, gene.info) {
    standardGeneric("ReplaceGeneInfo")
})

setMethod("ReplaceGeneInfo", signature("EMSet"), function(object, gene.info) {
    # Check replacement is valid.
    replacement <- object
    replacement@GeneInformation <- gene.info
    errors <- c()
    errors <- CheckGeneInformation(replacement, errors)
    if (length(errors) > 0) {
        stop("New Gene Information is invalid. Please check your data frame and try again.")
    } else {
        return(replacement)
    }
})

#' ReplaceExpressionMatrix
#'
#' Replace the expression matrix in a \linkS4class{EMSet} with a new expression matrix and re-calculate its metrics.
#'
#' @param object The \linkS4class{EMSet} you would like to update.
#' @param expression.matrix Expression matrix in matrix form
#' @include ascend_objects.R
#' @export
setGeneric(name = "ReplaceExpressionMatrix", def = function(object, expression.matrix) {
    standardGeneric("ReplaceExpressionMatrix")
})

setMethod("ReplaceExpressionMatrix", signature("EMSet"), function(object, expression.matrix) {
    # Replace the matrix
    object@ExpressionMatrix <- ConvertMatrix(expression.matrix, format = "sparseMatrix")

    # Check that it's okay before replacing.
    errors <- c()
    errors <- CheckExpressionMatrix(object, errors)

    if (length(errors) > 0) {
        stop("Please check your expression matrix.")
    } else {
        object <- GenerateMetrics(object)
        object <- SyncSlots(object)
        return(object)
    }
})
