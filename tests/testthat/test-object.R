# Tests for EMSet creation and modification
context("EMSet constructor functions")

# Test newEMSet function, ensure it outputs an EMSet
test_that("newEMSet works without controls.", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  assay_list <- list(counts = test_matrix)
  
  # Initiate test
  expect_match(class(newEMSet(assays = assay_list)), "EMSet")
  
})

# Test newEMSet works with controls
test_that("NewEMSet works - with controls", {
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
  expect_match(class(newEMSet(assays = list(counts = test_matrix), controls = controls)), "EMSet")  
})

# Test if colInfo works
test_that("NewEMSet works - colInfo loading", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids

  colInfo <- data.frame(cell_barcode = cell_ids, condition = "Test")
    
  # Initiate test
  expect_match(class(newEMSet(assays = list(counts = test_matrix), colInfo = colInfo)), "EMSet")  
})

# Test if rowInfo works
test_that("NewEMSet works - rowInfo loading", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  rowInfo <- data.frame(gene_id = gene_ids, condition = "Test")
  
  # Initiate test
  expect_match(class(newEMSet(assays = list(counts = test_matrix), rowInfo = rowInfo)), "EMSet")  
})
