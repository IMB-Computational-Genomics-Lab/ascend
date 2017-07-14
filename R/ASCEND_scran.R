#' ConvertToScater
#'
#' Convert a \linkS4class{AEMSet} object into a \linkS4class{SCESet} for use with \pkg{scater} and \pkg{scran}.
#' In order to use this function, you must have mitochondrial and ribosomal genes in your expression data.
#'
#' @param object A \linkS4class{AEMSet} object
#' @param control.list Optional - a named list containing mitochondrial and ribosomal genes
#'
ConvertToScater <- function(object, control.list = list()){
  # Prepare control list
  control.names <- c("Mt", "Rb")

  # Retrieve controls from AEMSet object if user hasn't supplied a list
  if (length(control.list) == 0){
    control.list <- object@Controls
  }

  # Verify Mt and Rb are names in supplied list
  if ( !all(control.names %in% names(control.list)) ){
    stop("Please make sure you have supplied mitochondrial and ribosomal gene identifiers in a named list.")
  }

  # Verify genes are present in the list
  if ( !all(sapply(control.list, function(x) length(x) > 0)) ){
    stop("Please make sure you have supplied genes to the mitochondrial and ribosomal gene lists.")
  }

  # Convert AEMSet to SCESet
  expression.matrix <- GetExpressionMatrix(object)
  sce.obj <- scater::newSCESet(countData = expression.matrix)
  sce.obj <- scater::calculateQCMetrics(sce.obj, feature_controls = control.list)
  return(sce.obj)
}

#' ConvertScaterToASCEND
#'
#' Loads data from a SCATER object to a pre-existing ASCEND object.
#' @param SCESet A \linkS4class{SCESet} from \pkg{scater}
#' @param AEMSet An \linkS4class{AEMSet}
#'
ConvertScaterToASCEND <- function(SCESet, AEMSet){
  # Retrieve counts from SCESet
  expression.matrix <- scater::counts(SCESet)

  # Convert matrix to sparse
  sparse.matrix <- LoadSparseMatrix(expression.matrix)

  # Load it back into AEMSet
  AEMSet@ExpressionMatrix <- sparse.matrix

  # Run Metrics again
  AEMSet <- GenerateMetrics(AEMSet)

  # Return updated object
  return(AEMSet)
}

