#' scranCellCycle
#'
#' Wrapper for \pkg{scran}'s cell cycle functions. Please ensure you are using
#' *ensembl_id* as your rownames in this dataset. Also ensure you are using
#' mitochondrial and ribosomal genes as controls.
#'
#' @param object An \linkS4class{AEMSet} object.
#' @param training.set A training dataset containing pairs of marker genes.
#' @export
#'
scranCellCycle <- function(object, training.set) {
  # Cell Cycle Testing
  expression.matrix <- GetExpressionMatrix(object, format = "matrix")

  # Run cyclone with cell cycle assignments as a vector
  cc.assignments <- scran::cyclone(expression.matrix, pairs = training.set, rownames(expression.matrix), BPPARAM = BiocParallel::bpparam())

  cell.info <- GetCellInfo(object)
  cell.info$phase <- cc.assignments$phases

  object <- ReplaceCellInfo(object, cell.info)
  return(object)
}

#' ConvertToScater
#'
#' Convert a \linkS4class{AEMSet} object into a \linkS4class{SCESet} for use with
#' \pkg{scater} and \pkg{scran}. In order to use this function, you must have
#' mitochondrial and ribosomal genes in your expression data.
#'
#' @param object An \linkS4class{AEMSet} object.
#' @param control.list Optional - a named list containing mitochondrial and ribosomal genes.
#' @export
#'
ConvertToScater <- function(object, control.list = list()) {
    # Prepare control list
    control.names <- c("Mt", "Rb")

    # Retrieve controls from AEMSet object if user hasn't supplied a list
    if (length(control.list) == 0) {
        control.list <- object@Controls
    }

    if (length(control.list) == 0) {
        # Verify Mt and Rb are names in supplied list
        if (!all(control.names %in% names(control.list))) {
            stop("Please make sure you have supplied mitochondrial and ribosomal gene identifiers in a named list.")
        }

        # Verify genes are present in the list
        if (!all(sapply(control.list, function(x) length(x) > 0))) {
            stop("Please make sure you have supplied genes to the mitochondrial and ribosomal gene lists.")
        }

        expression.matrix <- GetExpressionMatrix(object, "data.frame")
        sce.obj <- scater::newSCESet(countData = expression.matrix)
        sce.obj <- scater::calculateQCMetrics(sce.obj, feature_controls = control.list)
    } else {
        expression.matrix <- GetExpressionMatrix(object, "data.frame")
        sce.obj <- scater::newSCESet(countData = expression.matrix)
    }


    # Convert AEMSet to SCESet
    return(sce.obj)
}

#' ConvertScaterToASCEND
#'
#' Loads data from a SCATER object to a pre-existing ASCEND object.
#' @param SCESet A \linkS4class{SCESet} from \pkg{scater}
#' @param AEMSet An \linkS4class{AEMSet}
#' @export
#'
ConvertScaterToASCEND <- function(SCESet, AEMSet) {
    # Retrieve counts from SCESet
    expression.matrix <- scater::counts(SCESet)

    # Convert to sparse, re-run metrics and add to slot
    AEMSet <- ReplaceExpressionMatrix(AEMSet, expression.matrix)

    # Sync up the object
    AEMSet <- SyncSlots(AEMSet)

    # Return updated object
    return(AEMSet)
}
