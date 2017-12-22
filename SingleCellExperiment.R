devtools::load_all(pkg = "~/CodeRepositories/ascend-dev")
library(Matrix)
library(SingleCellExperiment)

# Defined by user
cellranger.path <- "~/Data/IPSCRetina_scRNA_Aggr_V2/"
genome.name <- "GRCh38p7"

# Called by LoadCellRanger
loaded.data <- LoadCellRangerData(cellranger.path, genome.name)
exprs.mtx <- loaded.data$ExpressionMatrix
genes <- loaded.data$GeneInformation
batch <- loaded.data$CellInformation
controls <- loaded.data$Controls

# Add control information to gene information
genes$control <- rep("Feature", nrow(genes))
control.indexes <- sapply(controls, function(x) which(genes[,1] %in% x), USE.NAMES = TRUE)
for (control in names(control.indexes)){
  indices <- control.indexes[[control]];
  genes$control[indices] <- control
}

# Generate a SingleCellExperiment object
## Load matrix into SingleCellExperiement
expression.matrix <- as.matrix(exprs.mtx)


sce <- SingleCellExperiment(assays = list(counts = expression.matrix))
colData(sce) <- as(batch, "DataFrame")
rowData(sce) <- as(genes, "DataFrame")

# Prototype EMSet definition
setClass("PEMSet",
         slots = list(log = "list", 
                      clustering = "list", 
                      controls = "list"),
         contains = "SingleCellExperiment")

pem.set <- new("PEMSet", sce)
pem.set@controls <- controls

# Populate with metrics
counts <- counts(pem.set)
cell.info <- colData(sce)
gene.info <- rowData(sce)

rownames(counts) <- gene.info[,1]
colnames(counts) <- cell.info[,1]

# Cell Metrics
total.counts <- colSums(counts)
genes.per.cell <- colSums(counts != 0)
total.expression <- sum(counts)

# Gene Metrics
counts.per.gene <- rowSums(counts)
average.counts <- rowMeans(counts)
cells.per.gene <- rowSums(counts != 0)
mean.gene.expression <- rowMeans(counts)

# Top gene expression
sorted.counts.per.gene <- sort(counts.per.gene, decreasing = TRUE)
top.gene.list <- names(sorted.counts.per.gene)
sorted.exprs.matrix <- counts[top.gene.list, ]
top.genes.percentage <- 100 * sorted.counts.per.gene/total.expression

# For controls
for (control in names(controls)){
  # Grab indices of genes that belong to this group
  indices <- which(gene.info$control == control)
  control.transcript.total.counts <- colSums(counts[indices,])
  print(control.transcript.total.counts)
}