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
  
  # Rejoin matrix
  
})