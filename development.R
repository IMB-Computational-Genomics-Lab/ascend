library(SingleCellExperiment)
library(Matrix)
library(scater)
# Load data from ascend package
load(system.file("extdata", "RGC_scRNASeq.RData", package = "ascend"))

# User input, as per ascend
dense.matrix <- as.matrix(matrix)

# Cell Information Prep
batch.information <- unlist(as.numeric(lapply(
  strsplit(as.character(barcodes$V1), "-"), `[`, 2)))
colnames(barcodes) <- c("cell_barcode")
barcodes$batch <- as.numeric(batch.information)

# Gene Information Prep
colnames(genes) <- c("ensembl_id", "gene_symbol")
genes <- genes[,c("gene_symbol", "ensembl_id")]
gene.names <- make.unique(as.vector(genes$gene_symbol))
rownames(matrix) <- gene.names
genes$gene_symbol <- gene.names

# Define Controls
mito.genes <- rownames(matrix)[grep("^MT-", rownames(matrix), ignore.case = TRUE)]
ribo.genes <- rownames(matrix)[grep("^RPS|^RPL", rownames(matrix), ignore.case = TRUE)]
genes$control <- genes[,1] %in% c(mito.genes, ribo.genes)
genes$control_group <- rep("Feature", nrow(genes))
genes$control_group[which(genes[,1] %in% mito.genes)] <- "Mt"
genes$control_group[which(genes[,1] %in% ribo.genes)] <- "Rb"

sce <- SingleCellExperiment(assays = list(counts = dense.matrix))
rowData(sce) <- as(genes, "DataFrame")
colData(sce) <- as(barcodes, "DataFrame")
rownames(sce) <- genes$gene_symbol
colnames(sce) <- barcodes$cell_barcode

sce <- calculateQCMetrics(sce, feature_controls = list(Rb = grepl("^RPS|^RPL", rownames(matrix), ignore.case = TRUE), Mt = grepl("^MT-", rownames(matrix), ignore.case = TRUE)))








setClass("EMSetPrototype",  representation(clusterObjects = "list", log = "list") , contains = "SingleCellExperiment")

prototype <- new("EMSetPrototype", sce)

# Conver the matrix into other formats
expression.matrix <- counts(sce)
cpm(sce) <- apply(expression.matrix, 2, function(x)(x/sum(x))*1e6)
logcounts(sce) <- log2(cpm(sce) + 1)

# Calculate Control Metrics
total.counts <- colSums(expression.matrix)

# Get Controls, subset out features
gene.info <- rowData(sce)
cell.info <- colData(sce)

# Get Control Information
control.info <- gene.info[which(gene.info$control),]
control.group.names <- unique(gene.info$control_group)
control.group.names <- control.group.names[-(which("Feature" == control.group.names))]

cell.info$total_counts <- total.counts

calculateControlCounts <- function(control, gene.info = NULL, expression.matrix = NULL){
  control.genes <- gene.info[which(gene.info$control_group == control), 1]
  control.transcript.counts <- expression.matrix[control.genes, ]
  control.transcript.total.counts <- Matrix::colSums(control.transcript.counts)
  control.pt.matrix <- (control.transcript.total.counts/total.counts) * 100
  output <- list(ControlTranscriptCounts = control.transcript.total.counts, PercentageTotalCounts = control.pt.matrix)
  return(output)
}

qc_list <- bplapply(control.group.names, calculateControlCounts, gene.info = gene.info, expression.matrix = expression.matrix)
names(qc_list) <- control.group.names


# Just do a normal loop for now, see if you can pack this into a bpvec
#for (control in control.group.names){
# control.genes <- gene.info[which(gene.info$control_group == control), 1]
# control.transcript.counts <- expression.matrix[control.genes, ]
# control.transcript.total.counts <- Matrix::colSums(control.transcript.counts)
# control.pt.matrix <- (control.transcript.total.counts/total.counts) * 100
# return(list(ControlTranscriptCounts = control.transcript.total.counts, PercentageTotalCounts = control.pt.matrix))
#}