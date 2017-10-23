context("EMSet subsetting methods")

test_that("Test if SubsetCondition returns an EMSet with the correct number of cells", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Define Cell Information
  cell.information <- data.frame(cell_barcode = cell.ids, batch = rep(1, ncol(test.matrix)), condition = rep(FALSE, ncol(test.matrix)))
  targets <- which(colSums(test.matrix[controls$Control,]) > 1000)
  cell.information[targets, "condition"] <- TRUE
  
  # Define Controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls, CellInformation = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_equal(ncol((SubsetCondition(em.set, "condition")@ExpressionMatrix)), length(targets))
})

test_that("Test if SubsetBatch returns an EMSet with the correct number of cells", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  batch.1.length <- ncol(test.matrix) - 70
  batch.2.length <- 70
  cell.ids.1 <- sapply(1:batch.1.length, function(x) paste0("Cell", x, "-1"))
  cell.ids.2 <- sapply(1:batch.2.length, function(x) paste0("Cell", x, "-2"))
  cell.ids <- c(cell.ids.1, cell.ids.2)
  batches <- unlist(as.numeric(lapply(strsplit(as.character(cell.ids), "-"), `[`, 2)))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  
  # Define Cell Information
  cell.information <- data.frame(cell_barcode = cell.ids, batch = batches, condition = rep(FALSE, ncol(test.matrix)))
  targets <- which(colSums(test.matrix[controls$Control,]) > 1000)
  cell.information[targets, "condition"] <- TRUE
  
  # Define Controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls, CellInformation = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_equal(ncol((SubsetBatch(em.set, 2)@ExpressionMatrix)), batch.2.length)
})

test_that("Check if SubsetCells subsets the right set of cells", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  batch.1.length <- ncol(test.matrix) - 70
  batch.2.length <- 70
  cell.ids.1 <- sapply(1:batch.1.length, function(x) paste0("Cell", x, "-1"))
  cell.ids.2 <- sapply(1:batch.2.length, function(x) paste0("Cell", x, "-2"))
  cell.ids <- c(cell.ids.1, cell.ids.2)
  batches <- unlist(as.numeric(lapply(strsplit(as.character(cell.ids), "-"), `[`, 2)))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  select.list <- sample(cell.ids, 10, replace = TRUE)
  
  # Define Cell Information
  cell.information <- data.frame(cell_barcode = cell.ids, batch = batches)
  
  # Define Controls
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls, CellInformation = cell.information)
  
  # SubsetCondition - expect number of targets to be equal
  expect_equal(colnames((SubsetCells(em.set, select.list))@ExpressionMatrix), select.list)
})

