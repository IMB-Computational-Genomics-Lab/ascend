## ----SetupNix, eval = FALSE----------------------------------------------
#  library(BiocParallel)
#  ncores <- parallel::detectCores() - 1
#  register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

## ----SetupWin, eval = FALSE----------------------------------------------
#  library(BiocParallel)
#  workers <- 3 # Number of cores on your machine - 1
#  register(SnowParam(workers = workers,
#                     type = "SOCK",
#                     progressbar = TRUE), default = TRUE)

## ----LoadDevAscend, echo = FALSE, message = FALSE------------------------
code_dir <- dirname(getwd())
devtools::load_all(code_dir)
library(BiocParallel)
# For single-core systems
ncores <- parallel::detectCores() - 1
if (ncores < 1){
  ncores <- 1
}
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)
load(url("https://github.com/IMB-Computational-Genomics-Lab/ascend-guide/blob/master/data/sample_data/RGC_scRNASeq.rdata?raw=true"))

## ----Loadascend, eval = FALSE--------------------------------------------
#  devtools::load_all("ascend/")

## ----LoadascendGit, eval = FALSE-----------------------------------------
#  devtools::install_github("IMB-Computational-Genomics-Lab/ascend")
#  library(ascend)

## ----echo = FALSE, results = "asis"--------------------------------------
library(knitr)
kable(data.table(cell_barcode = c("Cell1-1", "Cell2-1", "...", "...", "Cellx-x"), 
                 batch = c("1", "1", "...", "...", "x")))

## ----echo = FALSE, results = "asis"--------------------------------------
library(knitr)
kable(data.table(gene_identifier1 = c("GENE1", "GENE2", "...", "...", "GENEX"), 
                 gene_identifier2 = c("IDENTIFIER1", "IDENTIFIER2", "...", "...", 
                                      "IDENTIFIERX")))

## ----LoadData, echo = FALSE----------------------------------------------
# Quick prep of data to use as examples. Duplicate of what's done in 
# ascend_Tutorial.md
batch.information <- lapply(strsplit(as.character(barcodes$V1), "-"), `[`, 2)
colnames(barcodes) <- c("cell_barcode")
barcodes$batch <- as.numeric(batch.information)
colnames(matrix) <- barcodes[,1]
colnames(genes) <- c("ensembl_id", "gene_name")
genes <- genes[,c("gene_name", "ensembl_id")]
gene.names <- make.unique(as.vector(genes$gene_name))
rownames(matrix) <- gene.names
genes$gene_name <- gene.names
mito.genes <- rownames(matrix)[grep("^MT-", rownames(matrix), 
                                    ignore.case = TRUE)]
ribo.genes <- rownames(matrix)[grep("^RPS|^RPL", rownames(matrix), 
                                    ignore.case = TRUE)]
control.list <- list(Mt = mito.genes, Rb = ribo.genes)
expression.matrix <- matrix
gene.information <- genes
cell.information <- barcodes

## ----NewEMSet------------------------------------------------------------
em.set <- NewEMSet(ExpressionMatrix = expression.matrix, 
                   CellInformation = cell.information, 
                   GeneInformation = gene.information, 
                   Controls = control.list)
em.set

## ----ReplaceCellInfo, eval = FALSE---------------------------------------
#  updated.em.set <- ReplaceCellInfo(em.set, new.cell.info)

## ----UpdateControls, eval = FALSE----------------------------------------
#  old.controls <- GetControls(em.set)
#  new.controls <- c(old.controls, list(ERCC = c("ERCC-00031",
#                                                "ERCC-00017",
#                                                "ERCC-00024")))
#  updated.em.set <- UpdateControls(em.set, new.controls)

## ----ReplaceExpressionMatrix,  eval = FALSE------------------------------
#  updated.em.set <- ReplaceExpressionMatrix(em.set, new.matrix)

## ----SubsetBatch, eval = FALSE-------------------------------------------
#  subset.batch <- SubsetBatch(em.set, batches = c("1", "2"))

## ----SubsetCluster,  eval = FALSE----------------------------------------
#  subset.clusters <- SubsetCluster(em.set, clusters = c("2", "3"))

## ----SubsetCondition, eval = FALSE---------------------------------------
#  thy1.set <- SubsetCondition(em.set, condition = "phase", subconditions = c("G1", "G2M"))

## ----SubsetCells, eval = FALSE-------------------------------------------
#  # Retrieve cell information from an EMSet
#  cell.info <- GetCellInfo(em.set)
#  
#  # Randomly sample 100 cell barcodes to isolate.
#  cell.barcodes <- sample(cell.info$cell_barcode, 100, replace = FALSE)
#  hundred.cell.set <- SubsetCells(em.set, cell_barcodes = cell.barcodes)

## ----LoadCompleteData, echo = FALSE--------------------------------------
em.set <- readRDS(url("https://github.com/IMB-Computational-Genomics-Lab/ascend-guide/blob/master/data/sample_data/RGC_scRNASeq.rds?raw=true"))

## ----ColourByBatch-------------------------------------------------------
library(ggplot2)
# PCA PLOT
pca.plot <- PlotPCA(em.set, dim1 = 1, dim2 = 2, condition = "cluster")
pca.plot <- pca.plot + scale_color_manual(values=c("#bb5f4c", "#8e5db0", "#729b57")) 
pca.plot <- pca.plot + ggtitle("PCA Plot", subtitle = "Labelled by batch")

# MDS PLOT
mds.plot <- PlotMDS(em.set, PCA = FALSE, dim1 = 1, dim2 = 2, condition = "cluster") 
mds.plot <- mds.plot + scale_color_manual(values=c("#bb5f4c", "#8e5db0", "#729b57")) 
mds.plot <- mds.plot + ggtitle("MDS Plot", subtitle = "Labelled by batch")

# TSNE PLOT
tsne.plot <- PlotTSNE(em.set, 
                      PCA = TRUE, 
                      condition = "cluster", 
                      seed = 0, perplexity = 30, theta = 0.5) 
tsne.plot <- tsne.plot + scale_color_manual(values=c("#bb5f4c", "#8e5db0", "#729b57")) 
tsne.plot <- tsne.plot + ggtitle("TSNE Plot", subtitle = "Labelled by batch")

print(pca.plot)
print(mds.plot)
print(tsne.plot)

## ----ConvertToSingleCellExperiment---------------------------------------
controls <- GetControls(em.set)
sce.object <- ConvertToSCE(em.set, control.list = controls)
sce.object

## ----ConvertToSCESet, eval = FALSE---------------------------------------
#  controls <- GetControls(em.set)
#  sce.set <- ConvertToSCESet(em.set, control.list = controls)

