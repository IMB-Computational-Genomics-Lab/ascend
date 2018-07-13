context("EMSet substter functions.")

test_that("subsetCondition returns an EMSet with the correct number of cells", {
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
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  em.set <- newEMSet(assays = list(counts = test.matrix), controls = controls, colInfo = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_equal(ncol((subsetCondition(em.set, by = "condition", conditions = list(condition = "TRUE")))), length(targets))
})

test_that("Check if excludeControls subsets the right number of genes", {
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
  control.genes2 <- control.genes2[!(which(control.genes2 %in% control.genes1))]
  controls <- list(control1 = control.genes1, control2 = control.genes2)
  em.set <- newEMSet(assays = list(counts = test.matrix), controls = controls, colInfo = cell.information)
  
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
  em.set <- newEMSet(assays = list(counts = test.matrix), controls = controls, colInfo = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_false(progressLog(excludeControl(em.set, control = c("control1", "control2")))$controls)
})
