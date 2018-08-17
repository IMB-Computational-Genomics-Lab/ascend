context("EMSet-specific slot retrieval functions.")

# Test colInfo method
test_that("colInfo retrieval", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  colInfo <- data.frame(cell_barcode = cell_ids, condition = "Test")
  em_set <- newEMSet(assays = list(counts = test_matrix), colInfo = colInfo)
  rownames(colInfo) <- cell_ids
  # Initiate test
  expect_equal(colInfo(em_set), S4Vectors::DataFrame(colInfo))
})

test_that("rowInfo retrieval", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  rowInfo <- data.frame(gene_id = gene_ids, condition = "Test")
  em_set <- newEMSet(assays = list(counts = test_matrix), rowInfo = rowInfo)
  rownames(rowInfo) <- gene_ids
  
  # Initiate test
  expect_equal(rowInfo(em_set), S4Vectors::DataFrame(rowInfo))
})

test_that("progressLog set", {
  # Generate dummy expression matrix
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  progressLog <- list(Test = "Successful!")
  em_set <- newEMSet(assays = list(counts = test_matrix))
  
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
  em_set <- newEMSet(assays = list(counts = test_matrix))
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
  em_set <- newEMSet(assays = list(counts = test_matrix))
  
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
  em_set <- newEMSet(assays = list(counts = test_matrix))
  clusterAnalysis(em_set) <- clusterAnalysis
  
  expect_equal(clusterAnalysis(em_set), clusterAnalysis)
})

test_that("addControls test", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Define Cell Information
  cell.information <- data.frame(cell_barcode = cell.ids, batch = rep(1, ncol(test.matrix)), condition = rep(FALSE, ncol(test.matrix)))
  targets <- which(cell.information$cell_barcode %in% sample(cell.information$cell_barcode, 20, replace = FALSE))
  cell.information[targets, "condition"] <- TRUE
  
  # Define Controls
  control.genes1 <- sample(gene.ids, 10, replace = FALSE)
  control.genes2 <- sample(gene.ids, 10, replace = FALSE)
  if(all(control.genes2 %in% control.genes1)){
    control.genes2 <- control.genes2[!(which(control.genes2 %in% control.genes1))]    
  }
  controls <- list(control1 = control.genes1, control2 = control.genes2)
  em.set <- newEMSet(assays = list(counts = test.matrix), colInfo = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_true(progressLog(addControlInfo(em.set, controls = controls))$controls)
})
