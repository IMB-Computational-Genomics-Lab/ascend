################################################################################
#
# ascend_utilities.R
# description: Functions that help other functions. 
#
################################################################################

#' unLogMatrix
#' 
#' Convenience function - converts a matrix that has been converted to logcounts
#' to un-logged counts.
#' 
#' @param x Logged matrix to convert into unlogged matrix
#' @examples
#' # Randomly populate a matrix with 0s and 1s
#' normal_matrix <- matrix(sample(0:1,10*10, replace=TRUE), 10 , 10)
#' 
#' # Log2 + 1 values in normal_matrix
#' logged_matrix <- log2(normal_matrix + 1)
#' 
#' # Return values to normal
#' unlogged_matrix <- unLog2Matrix(logged_matrix)
#' 
unLog2Matrix <- function(x){
  # Convert to matrix
  x <- as(x, "matrix")
  
  # Unlog the matrix
  unlogged_matrix <- 2^x
  
  # Subtract pseudocount of 1
  unlogged_matrix_sub_1 <- unlogged_matrix - 1

  # Make negative values 0
  unlogged_matrix_sub_1[unlogged_matrix_sub_1 < 0] <- 0
  
  # Make infinite values 0
  unlogged_matrix_sub_1[!is.finite(unlogged_matrix_sub_1)] <- 0
  
  return(unlogged_matrix_sub_1)
}
#' joinPaths
#' 
#' Convenience function - joins a list of strings to form a full path. 
#' Disregards "/" at the end of a substring to ensure final full path is not
#' malformed.
#' 
#' @param x A list of substrings to combine into a full path.
#' @examples
#' sub_paths <- c("Path", "To/", "Me")
#' full_path <- joinPaths(sub_paths)
#' print(full_path)
#' 
joinPaths <- function(x){
  if (length(x) > 1){
    x <- gsub("/$", "", x)
    path <- do.call("file.path", as.list(x))
    return(path)
  } else{
    return(x)
  }
}

#' fileCheck 
#' This little snippet is designed in a way in which you pair it up with a call 
#' to process the file. This checks if the file exists.
#' 
#' @param x Name of the file you want to check.
#' @return Returns a message if the file doesn't exist.
#' 
fileCheck <- function(x) {
  if (!(file.exists(x))) {
    stop(sprintf("%s is missing", x))
  } else {
    return(FALSE)
  }
}
