context("EMSet methods")

# Check if GenerateMetrics works
test_that("GenerateControlMetrics function check", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  total.counts <- colSums(test.matrix)
  
  # Run test
  expect_true(is.list(GenerateControlMetrics("Control", expression.matrix = test.matrix, control.list = controls, total.counts = total.counts)))
})

# Check ConvertGeneAnnotation
test_that("ConverGeneAnnotation successfully switches", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  gene.ids.2 <- sapply(1:nrow(test.matrix), function(x) paste0("AlternativeGene", x))
  gene.information <- data.frame(gene1 = gene.ids, gene2 = gene.ids.2) 
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, GeneInformation = gene.information)
  
  # Switch annotations
  switched.em.set <- ConvertGeneAnnotation(em.set, "gene1", "gene2")
  
  # Match switched annotations to second column of gene information dataframe
  expect_true(identical(rownames(switched.em.set@ExpressionMatrix), as.vector(gene.information$gene2)))
})

test_that("ConvertGeneAnnotation function check", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  gene.ids.2 <- sapply(1:nrow(test.matrix), function(x) paste0("AlternativeGene", x))
  gene.information <- data.frame(gene1 = gene.ids, gene2 = gene.ids.2) 
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, GeneInformation = gene.information)
  
  # Initiate test
  expect_match(class(ConvertGeneAnnotation(em.set, "gene1", "gene2")), "EMSet")
})

# Check ExcludeControl
test_that("Check if ExcludeControl works", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Initiate test
  expect_equal(length((ExcludeControl(em.set, "Control"))@Controls), 0)
})

# Check DisplayLog
test_that("Check if DisplayLog works", {
  # Generate a test EMSet
  test.matrix <- matrix(rnbinom(1000*200, mu=2^runif(1000, 3, 10), size=2), nrow=1000)
  cell.ids <- sapply(1:ncol(test.matrix), function(x) paste0("Cell", x))
  gene.ids <- sapply(1:nrow(test.matrix), function(x) paste0("Gene", x))
  colnames(test.matrix) <- cell.ids
  rownames(test.matrix) <- gene.ids
  control.genes <- sample(gene.ids, 10, replace = FALSE)
  controls <- list("Control" = control.genes)
  em.set <- NewEMSet(ExpressionMatrix = test.matrix, Controls = controls)
  
  # Iniitiate test
  expect_equal(em.set@Log, list(Controls = TRUE))
})