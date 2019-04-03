context("EMSet-specific slot retrieval functions.")

# Test colInfo method
test_that("colInfo retrieval", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  colInfo <- S4Vectors::DataFrame(cell_barcode = cell_ids, batch = 1, condition = "Test")

  # Create EMSet
  em_set <- EMSet(list(counts = test_matrix), colInfo = colInfo)
  rownames(colInfo) <- cell_ids
  
  # Initiate test
  expect_equal(colInfo(em_set), colInfo)
})

test_that("rowInfo retrieval", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  rowInfo <- S4Vectors::DataFrame(gene_id = gene_ids, condition = "Test")
  em_set <- EMSet(list(counts = test_matrix), rowInfo = rowInfo)
  
  # Initiate test
  rownames(rowInfo) <- gene_ids
  expect_equal(rowInfo(em_set), S4Vectors::DataFrame(rowInfo))
})

test_that("progressLog set", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  em_set <- EMSet(list(counts = test_matrix))
  
  # Initiate test
  expect_message(progressLog(em_set), NA)
})

test_that("progressLog retrieval", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  progressLog <- list(Test = "Successful!")
  em_set <- EMSet(list(counts = test_matrix))
  progressLog(em_set) <- progressLog
  
  expect_equal(progressLog(em_set), progressLog)
})

test_that("clusterAnalysis set", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  clusterAnalysis <- list(Test = "Successful!")
  em_set <- EMSet(list(counts = test_matrix))
  
  # Initiate test
  expect_message(clusterAnalysis(em_set), NA)
})

test_that("clusterAnalysis retrieval", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  clusterAnalysis <- list(Test = "Successful!")
  em_set <- EMSet(list(counts = test_matrix))
  clusterAnalysis(em_set) <- clusterAnalysis
  
  expect_equal(clusterAnalysis(em_set), clusterAnalysis)
})

test_that("addControls test", {
  # Generate a test EMSet
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  # Define Cell Information
  cell_information <- S4Vectors::DataFrame(cell_barcode = cell_ids, batch = rep(1, ncol(test_matrix)), condition = rep(FALSE, ncol(test_matrix)))
  targets <- which(cell_information$cell_barcode %in% sample(cell_information$cell_barcode, 20, replace = FALSE))
  cell_information[targets, "condition"] <- TRUE
  
  # Define Controls
  control_genes1 <- sample(gene_ids, 10, replace = FALSE)
  control_genes2 <- sample(gene_ids, 10, replace = FALSE)
  
  if(all(control_genes2 %in% control_genes1)){
    control_genes2 <- control_genes2[!(which(control_genes2 %in% control_genes1))]    
  }
  
  controls <- list(control1 = control_genes1, control2 = control_genes2)
  
  # Build an EMSet
  em_set <- EMSet(list(counts = test_matrix), colInfo = cell_information)
  controls(em_set) <- controls
  # SubsetCondition - expect number of targets to be equal
  expect_true(progressLog(em_set)$controls)
})
