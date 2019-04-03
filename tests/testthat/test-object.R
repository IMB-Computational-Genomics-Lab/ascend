# Tests for EMSet creation and modification
context("EMSet constructor functions")

# Test EMSet function, ensure it outputs an EMSet
test_that("EMSet works without controls.", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  assay_list <- list(counts = test_matrix)
  
  # Initiate test
  expect_match(class(EMSet(assay_list)), "EMSet")
  
})

# Test EMSet works with controls
test_that("EMSet works - with controls", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  # Generate dummy controls
  control_genes <- sample(gene_ids, 10, replace = FALSE)
  controls <- list(Control = control_genes)
  
  # Initiate test
  expect_match(class(EMSet(list(counts = test_matrix), controls = controls)), "EMSet")  
})

# Test if colInfo works
test_that("EMSet works - colInfo loading", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids

  colInfo <- S4Vectors::DataFrame(cell_barcode = cell_ids, batch = 1, condition = "Test")
    
  # Initiate test
  expect_match(class(EMSet(list(counts = test_matrix), colInfo = colInfo)), "EMSet")  
})

# Test if rowInfo works
test_that("EMSet works - rowInfo loading", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  rowInfo <- S4Vectors::DataFrame(gene_id = gene_ids, condition = "Test")
  
  # Initiate test
  expect_match(class(EMSet(list(counts = test_matrix), rowInfo = rowInfo)), "EMSet")  
})
