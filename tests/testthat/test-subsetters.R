context("EMSet substter functions.")

test_that("subsetCondition returns an EMSet with the correct number of cells", {
  # Generate a test EMSet
  test_matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell_ids <- sapply(1:ncol(test_matrix), function(x) paste0("Cell", x))
  gene_ids <- sapply(1:nrow(test_matrix), function(x) paste0("Gene", x))
  colnames(test_matrix) <- cell_ids
  rownames(test_matrix) <- gene_ids
  
  # Define Cell Information
  cell_information <- S4Vectors::DataFrame(cell_barcode = cell_ids, batch = rep(1, ncol(test_matrix)), condition = rep(FALSE, ncol(test_matrix)), row.names = cell_ids)
  targets <- which(cell_information$cell_barcode %in% sample(cell_information$cell_barcode, 20, replace = FALSE))
  cell_information[targets, "condition"] <- TRUE
  
  # Define Controls
  control_genes <- sample(gene_ids, 10, replace = FALSE)
  controls <- list("Control" = control_genes)
  em_set <- newEMSet(list(counts = test_matrix), controls = controls, colInfo = cell_information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_equal(ncol((subsetCondition(em_set, by = "condition", conditions = list(condition = "TRUE")))), length(targets))
})

test_that("Check if excludeControls subsets the right number of genes", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Define Cell Information
  cell.information <- S4Vectors::DataFrame(cell_barcode = cell.ids, batch = rep(1, ncol(test.matrix)), condition = rep(FALSE, ncol(test.matrix)), row.names = cell.ids)
  targets <- which(cell.information$cell_barcode %in% sample(cell.information$cell_barcode, 20, replace = FALSE))
  cell.information[targets, "condition"] <- TRUE
  
  # Define Controls
  control.genes1 <- sample(gene.ids, 10, replace = FALSE)
  control.genes2 <- sample(gene.ids, 10, replace = FALSE)
  control.genes2 <- control.genes2[!(which(control.genes2 %in% control.genes1))]
  controls <- list(control1 = control.genes1, control2 = control.genes2)
  em.set <- newEMSet(list(counts = test.matrix), controls = controls, colInfo = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_equal(nrow((excludeControl(em.set, control = "control1"))), nrow(em.set) - length(control.genes1))
})

test_that("Check if excludeControls of all controls returns control set to FALSE", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Define Cell Information
  cell.information <- S4Vectors::DataFrame(cell_barcode = cell.ids, batch = rep(1, ncol(test.matrix)), condition = rep(FALSE, ncol(test.matrix)), row.names = cell.ids)
  targets <- which(cell.information$cell_barcode %in% sample(cell.information$cell_barcode, 20, replace = FALSE))
  cell.information[targets, "condition"] <- TRUE
  
  # Define Controls
  control.genes1 <- sample(gene.ids, 10, replace = FALSE)
  control.genes2 <- sample(gene.ids, 10, replace = FALSE)
  if(all(control.genes2 %in% control.genes1)){
    control.genes2 <- control.genes2[!(which(control.genes2 %in% control.genes1))]    
  }
  controls <- list(control1 = control.genes1, control2 = control.genes2)
  em.set <- newEMSet(list(counts = test.matrix), controls = controls, colInfo = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_false(progressLog(excludeControl(em.set, control = c("control1", "control2")))$controls)
})
