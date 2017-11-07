# Tests for EMSet creation and modification
context("EMSet functions")

# Test for NewEMSet outputs an EMSet
test_that("NewEMSet works - no controls", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Initiate test
  expect_match(class(NewEMSet(ExpressionMatrix = test.matrix)), "EMSet")
})

test_that("NewEMSet works - with controls", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Initiate test
  expect_match(class(NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)), "EMSet")  
})

# Tests that ensure EMSet can be updated correctly
test_that("Check if UpdateControls works", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix)
  
  # Generate dummy controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list(Control = control.genes)
  
  # Initiate test
  expect_match(class(UpdateControls(em.set, controls)), "EMSet")
})

test_that("Check if ReplaceCellInfo works", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix)
  
  # Generate Dummy Replacement
  old.cell.info <- em.set@CellInformation
  new.cell.info <- old.cell.info
  new.cell.info$test <- rep("Test", nrow(new.cell.info))
  
  # Initiate test
  expect_match(class(ReplaceCellInfo(em.set, new.cell.info)), "EMSet")
})

test_that("Check if ReplaceGeneInfo works", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix)
  
  # Generate Dummy Replacement
  old.gene.info <- em.set@GeneInformation
  new.gene.info <- old.gene.info
  new.gene.info$test <- rep("Test", nrow(new.gene.info))
  
  # Initiate test
  expect_match(class(ReplaceGeneInfo(em.set, new.gene.info)), "EMSet")
})

test_that("Check if ReplaceExpressionMatrix works", {
  # Generate dummy expression matrix
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Generate dummy EMSet
  em.set <- NewEMSet(ExpressionMatrix = test.matrix)
  
  # Generate a new expression matrix
  replacement.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  new.cell.ids <- sapply(1:ncol(replacement.matrix), function(x) paste0("Cell", x))
  new.gene.ids <- sapply(1:nrow(replacement.matrix), function(x) paste0("Gene", x))
  
  colnames(replacement.matrix) <- new.cell.ids
  rownames(replacement.matrix) <- new.gene.ids

  # Test function
  # Initiate test
  expect_match(class(ReplaceExpressionMatrix(em.set, replacement.matrix)), "EMSet")
})