# Set up environment
devtools::load_all("/Users/a.senabouth/CodeRepositories/ASCEND")
library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar = TRUE), default = TRUE)

# Load 1 FlowSeq sample
file.dir <- "/Volumes/Anne's External HD/NeuroSpheres_scRNA_V1/NeuroSpheres_scRNA_Aggr_ExpressionMatrix_V1"
matrix <- read.csv(paste0(file.dir, "/", "NeuroSpheres_scRNA_ExpressionMatrix.csv"), header = TRUE, row.names = 1)
cell.barcodes <- colnames(matrix)
batch.info <- unlist(as.numeric(lapply(strsplit(as.character(cell.barcodes), "[.]"), `[`, 2)))
batch.info[1:5]

cell.info <- data.frame(cell_barcode = cell.barcodes, batch = batch.info)
gene.info <- data.frame(gene_symbol = rownames(matrix))
mito.genes <- rownames(matrix)[grep("^MT-", rownames(matrix), ignore.case = TRUE)]
ribo.genes <- rownames(matrix)[grep("^RPS|^RPL", rownames(matrix), ignore.case = TRUE)]
controls <- list(Mt = mito.genes, Rb = ribo.genes)

# Create AEMSet
aem.set <- NewAEMSet(ExpressionMatrix = matrix, CellInformation = cell.info, GeneInformation = gene.info, Controls = controls)
general.qc <- PlotGeneralQC(aem.set)

batch.normalised <- NormaliseBatches(aem.set)
batch.normalised.matrix <- GetExpressionMatrix(batch.normalised)
