context("Test ascend utilities")

test_that("Unlog function check", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids

  # Generate a log matrix, add 1 to mimic steps done in pipeline
  log.matrix <- log2(test.matrix) + 1
  
  # Unlog it with the Unlog function
  expect_true(is.matrix(UnLog2Matrix(log.matrix)))
})

test_that("Matrix chunking function check", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids

  # Chunk Matrix
  chunked.matrix <- ChunkMatrix(test.matrix, axis = 0, chunks = 10)
  
  # Rejoin matrix and sort back into original order
  combined.df <- do.call("rbind", chunked.matrix)
  combined.df <- combined.df[rownames(test.matrix), colnames(test.matrix)]
  
  # Expect they're equal
  expect_equal(combined.df, test.matrix)
})

test_that("JoinMatrices check", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  
  test.matrix.2 <- matrix(rnbinom(1000*100, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.2 <- sapply(1:ncol(test.matrix.2), function(x) paste0("Cell", x))
  colnames(test.matrix.2) <- cell.ids.2
  rownames(test.matrix.2) <- gene.ids
  
  joined.matrix <- JoinMatrices(list(test.matrix.1, test.matrix.2))
  
  # Check if dimensions make sense
  expect_equal(nrow(joined.matrix), 1000)
  expect_equal(ncol(joined.matrix), 300)
  
  # Check if all dimensions are in the expression matrix
  expect_true(all(rownames(test.matrix.1) %in% rownames(joined.matrix)))
  expect_true(all(colnames(test.matrix.1) %in% colnames(joined.matrix)))
  expect_true(all(rownames(test.matrix.2) %in% rownames(joined.matrix)))
  expect_true(all(colnames(test.matrix.2) %in% colnames(joined.matrix)))
})

test_that("ConvertMatrix function check - data.frame to matrix", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  test.matrix.1 <- as.data.frame(test.matrix.1)
  
  expect_true(is.matrix(ConvertMatrix(test.matrix.1, "matrix")))
})

test_that("ConvertMatrix function check - data.frame to sparse matrix", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  test.matrix.1 <- as.data.frame(test.matrix.1)
  
  expect_true(is(ConvertMatrix(test.matrix.1, "sparseMatrix"), "sparseMatrix"))
})

test_that("ConvertMatrix function check - matrix to data.frame", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  
  expect_true(is.data.frame(ConvertMatrix(test.matrix.1, "data.frame")))
})

test_that("ConvertMatrix function check - matrix to sparse matrix", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  
  expect_true(is(ConvertMatrix(test.matrix.1, "sparseMatrix"), "sparseMatrix"))
})

test_that("ConvertMatrix function check - sparse matrix to data.frame", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  
  # Convert to sparse
  test.matrix.1 <- Matrix::Matrix(test.matrix.1, sparse = TRUE)
  
  expect_true(is.data.frame(ConvertMatrix(test.matrix.1, "data.frame")))
})

test_that("ConvertMatrix function check - sparse matrix to matrix", {
  # Generate dummy matrices
  test.matrix.1 <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids.1 <- sapply(1:ncol(test.matrix.1), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix.1), function(x) paste0("Gene", x))
  colnames(test.matrix.1) <- cell.ids.1
  rownames(test.matrix.1) <- gene.ids
  
  # Convert to sparse
  test.matrix.1 <- Matrix::Matrix(test.matrix.1, sparse = TRUE)
  
  expect_true(is.matrix(ConvertMatrix(test.matrix.1, "matrix")))
})
