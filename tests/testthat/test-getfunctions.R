context("EMSet Get Functions")

test_that("GetExpressionMatrix returns a data frame", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Generate EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Expect output
  expect_true(is.data.frame(GetExpressionMatrix(em.set, format = "data.frame")))
})

test_that("GetExpressionMatrix returns a matrix", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Generate EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Expect output
  expect_true(is.matrix(GetExpressionMatrix(em.set, format = "matrix")))
})

test_that("GetControls returns a list", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Generate EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Expect output
  expect_true(is.list(GetControls(em.set)))
})

test_that("GetCellInfo returns a data frame", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Generate EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Expect output
  expect_true(is.data.frame(GetCellInfo(em.set)))
})

test_that("GetGeneInfo returns a data frame", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Generate EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Expect output
  expect_true(is.data.frame(GetGeneInfo(em.set)))
})

test_that("GetBatchMatrix returns a matrix", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Generate EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Expect output
  expect_true(is.matrix(GetBatchMatrix(em.set, 1)))
  
  # Expect output
  expect_true(is.matrix(GetBatchMatrix(em.set, "1")))
  
})

