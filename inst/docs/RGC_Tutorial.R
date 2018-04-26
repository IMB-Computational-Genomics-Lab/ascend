## ----Loadascend, message = FALSE, echo = FALSE---------------------------
# Load code from parent directory
code_dir <- dirname(getwd())
devtools::load_all(code_dir)
library(BiocParallel)
ncores <- parallel::detectCores() - 1
if (ncores < 1){
  ncores <- 1
}
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)
load(url("https://github.com/IMB-Computational-Genomics-Lab/ascend-guide/blob/master/data/sample_data/RGC_scRNASeq.rdata?raw=true"))

## ----LoadData, eval = FALSE----------------------------------------------
#  load("RGC_scRNASeq.RData")

## ----ReadExpressionMatrix------------------------------------------------
matrix <- as.data.frame(as.matrix(matrix))
matrix[1:5,1:5]

## ----ReadBarcodes--------------------------------------------------------
barcodes[1:5,]

## ----ProcessBatchID------------------------------------------------------
batch.information <- unlist(as.numeric(lapply(strsplit(as.character(barcodes$V1), "-"), `[`, 2)))
batch.information[1:5]

## ----BuildCellInfo-------------------------------------------------------
colnames(barcodes) <- c("cell_barcode")
barcodes$batch <- as.numeric(batch.information)
barcodes[1:5,]

## ----NameColumns---------------------------------------------------------
colnames(matrix) <- barcodes[,1]
matrix[1:5, 1:5]

## ----ReadGenes-----------------------------------------------------------
colnames(genes) <- c("ensembl_id", "gene_symbol")
genes[1:5,]

## ----SetRownames---------------------------------------------------------
genes <- genes[,c("gene_symbol", "ensembl_id")]
gene.names <- make.unique(as.vector(genes$gene_symbol))
rownames(matrix) <- gene.names
matrix[1:5, 1:5]

## ----ReplaceGeneNames----------------------------------------------------
genes$gene_symbol <- gene.names
genes[1:15,]

## ----DefineControls------------------------------------------------------
mito.genes <- rownames(matrix)[grep("^MT-", rownames(matrix), 
                                    ignore.case = TRUE)]
ribo.genes <- rownames(matrix)[grep("^RPS|^RPL", 
                                    rownames(matrix), 
                                    ignore.case = TRUE)]
controls <- list(Mt = mito.genes, Rb = ribo.genes)
controls

## ----BuildEMSet----------------------------------------------------------
em.set <- NewEMSet(ExpressionMatrix = matrix, GeneInformation = genes, 
                   CellInformation = barcodes, Controls = controls)
em.set

## ----AutoLoading, eval = FALSE-------------------------------------------
#  em.set <- LoadCellRanger("RGC_scRNASeq/", "GRCh38")

## ----LabelTHY1Cells------------------------------------------------------
cell.info <- GetCellInfo(em.set)
thy1.expression <- cell.info$batch
thy1.expression <- thy1.expression == 1
cell.info$THY1 <- thy1.expression
cell.info[1:5, ]

## ----LabelPOU4F2---------------------------------------------------------
# Create a list of transcript names
brn3.transcripts <- c("POU4F1", "POU4F2", "POU4F3")

# Extract expression matrix from the em.set as a data frame
expression.matrix <- GetExpressionMatrix(em.set, format = "data.frame")

# Extract rows from matrix belonging to these transcripts
brn3.transcript.counts <- expression.matrix[brn3.transcripts, ]

# Identify cells (columns) where transcript counts are greater than one
brn3.cells <- colSums(brn3.transcript.counts) > 0

# Add new information to cell information
cell.info$BRN3 <- brn3.cells

# View cell.info
cell.info[1:5,]

## ----ReplaceCellInfo-----------------------------------------------------
em.set <- ReplaceCellInfo(em.set, cell.info)

## ---- PlotGeneralQCChunk, fig.show = "hide", results="hide", include = FALSE----
raw.qc.plots <- PlotGeneralQC(em.set)

## ----CellCycle-----------------------------------------------------------
# Convert the EMSet's gene annotation to ENSEMBL IDs stored in the ensembl_id 
# column of the GeneInformation dataframe
em.set <- ConvertGeneAnnotation(em.set, "gene_symbol", "ensembl_id")

# Load scran's training dataset
training.data <- readRDS(system.file("exdata", "human_cycle_markers.rds", 
                                     package = "scran"))

# Run scranCellCycle
em.set <- scranCellCycle(em.set, training.data)

# View cell information
cell.info <- GetCellInfo(em.set)
cell.info[1:5, ]

# Convert annotation back to gene_symbol
em.set <- ConvertGeneAnnotation(em.set, "ensembl_id", "gene_symbol")

## ---- fig.align="center", fig.width = 4, fig.show = "hold"---------------
print(raw.qc.plots$LibSize)
print(raw.qc.plots$FeatureCountsPerCell)
print(raw.qc.plots$ControlPercentageTotalCounts$Mt)
print(raw.qc.plots$ControlPercentageTotalCounts$Rb)

## ----FilterByOutliers----------------------------------------------------
em.set <- FilterByOutliers(em.set, cell.threshold = 3, control.threshold = 3)

## ----ControlPercentagePlots, fig.width=4, fig.height=5.25, fig.align="center", fig.show = "hold"----
print(raw.qc.plots$ControlPercentageSampleCounts$Mt)
print(raw.qc.plots$ControlPercentageSampleCounts$Rb)

## ----GetControls---------------------------------------------------------
print(GetControls(em.set))

## ----FilterByControl-----------------------------------------------------
# Filter by mitochondrial genes
em.set <- FilterByControl(control.name = "Mt", pct.threshold = 20, em.set)
# Filter by ribosomal genes
em.set <- FilterByControl(control.name = "Rb", pct.threshold = 50, em.set)

## ----AverageGeneCountPlots, fig.width=4, fig.align="center"--------------
print(raw.qc.plots$AverageGeneCount)

## ----AverageGeneCountLogPlots, fig.align="center", fig.show = "hold", fig.width = 4----
print(raw.qc.plots$Log2AverageGeneCount)
print(raw.qc.plots$Log10AverageGeneCount)

## ----FilterLowAbundanceGenes, eval = FALSE-------------------------------
#  em.set <- FilterLowAbundanceGenes(em.set, pct.value = 1)

## ----CheckGenesOfInterest------------------------------------------------
expression.matrix <- GetExpressionMatrix(em.set, "data.frame")
brn3.transcripts <- c("POU4F1", "POU4F2", "POU4F3")
expression.matrix[brn3.transcripts, 
                  which(colSums(expression.matrix[brn3.transcripts,]) > 0)]

## ----DisplayLog----------------------------------------------------------
DisplayLog(em.set)

## ----CheckFiltering, fig.show = "hide", results="hide", include = FALSE----
filtered.qc.plots <- PlotGeneralQC(em.set)
print(filtered.qc.plots$LibSize)
print(filtered.qc.plots$FeatureCountsPerCell)

## ----NormaliseByRLE, eval = FALSE----------------------------------------
#  norm.set <- NormaliseByRLE(em.set)

## ----GetControlsScran----------------------------------------------------
print(GetControls(em.set))

## ----scranNormalise------------------------------------------------------
norm.set <- scranNormalise(em.set, quickCluster = FALSE, min.mean = 1e-5)

## ----PlotNormalisationQC-------------------------------------------------
norm.qc <- PlotNormalisationQC(original = em.set, normalised = norm.set, 
                               gene.list = c("GAPDH", "MALAT1"))

## ----NormLibsizePlot, fig.show = "hold", fig.align="center", fig.width = 4----
print(norm.qc$Libsize$Original)
print(norm.qc$Libsize$Normalised)

## ----NormScatter, fig.show="hold", fig.align="center", fig.width = 4.5, warning = FALSE----
print(norm.qc$GeneScatterPlots$GAPDH$Original)
print(norm.qc$GeneScatterPlots$GAPDH$Normalised)
print(norm.qc$GeneScatterPlots$MALAT1$Original)
print(norm.qc$GeneScatterPlots$MALAT1$Normalised)

## ----NormGenes, fig.show="hold", fig.align = "center", fig.width = 4.5, warning = FALSE----
print(norm.qc$GeneExpressionBoxplot$Original)
print(norm.qc$GeneExpressionBoxplot$Normalised)

## ----ControlRemovalPlot1, fig.width = 5, fig.height = 8.5, warning = FALSE----
print(filtered.qc.plots$TopGenes)

## ----ControlRemovalPlot2, fig.width = 5, fig.height = 4, fig.align="center", warning = FALSE----
top.20.plot <- PlotTopGeneExpression(norm.set, n = 20, controls = FALSE)
print(top.20.plot)

## ----ExcludeControls, eval = FALSE---------------------------------------
#  norm.set <- ExcludeControl(norm.set, "Mt")
#  norm.set <- ExcludeControl(norm.set, "Rb")

## ----ConfoundingFactors, eval = FALSE------------------------------------
#  cell.cycle.genes <- c("CDK4","CCND1","NOC2L","ATAD3C", "CCNL2")
#  em.set <- RegressConfoundingFactors(em.set, candidate.genes = cell.cycle.genes)

## ----DimReduction1-------------------------------------------------------
pca.set <- RunPCA(norm.set)

## ----PlotPCAVariance, fig.show="hold", fig.align="center", fig.width = 5----
pca.variance <- PlotPCAVariance(pca.set, n = 50)
print(pca.variance)

## ----ReduceDimensions----------------------------------------------------
pca.set <- ReduceDimensions(pca.set, n = 20)

## ----RunCOREfunctions, echo = FALSE, results = "asis"--------------------
library(knitr)
kable(data.table(Argument = c("conservative", "nres", "remove_outlier"), 
                 Description = c("Use conservative (more stable) clustering result (TRUE or FALSE). Default: TRUE", "Number of resolutions to test Default: 40", "Remove cells that weren't assigned a cluster with dynamicTreeCut. This is indicative of outlier cells within the sample. Default: TRUE")))

## ----RunCORE-------------------------------------------------------------
clustered.set <- RunCORE(pca.set, 
                         conservative = TRUE, 
                         nres = 40, 
                         remove_outlier = FALSE )

## ----PlotStabilityDendro, fig.width=4, fig.height=5, fig.show="hold", fig.align="center"----
PlotStabilityDendro(clustered.set)

## ----PlotStability, fig.width=5, fig.height=4, fig.show="hold", fig.align="center"----
PlotStability(clustered.set)

## ----GetRandMatrix-------------------------------------------------------
rand.matrix <- GetRandMatrix(clustered.set)
rand.matrix

## ----PlotDendrogram, fig.width=5, fig.height=4, fig.show="hold", fig.align="center"----
PlotDendrogram(clustered.set)

## ----GetCellInfo---------------------------------------------------------
cell.info <- GetCellInfo(clustered.set)
cell.info[1:5,]

## ----THY1DE--------------------------------------------------------------
thy1.de.result <- RunDiffExpression(clustered.set, 
                                    condition.a = "TRUE", 
                                    condition.b = "FALSE", 
                                    conditions = "THY1", 
                                    fitType = "local", 
                                    method = "per-condition")
thy1.de.result[1:10,]

## ----THY1DEplot----------------------------------------------------------
thy1.volcano.plot <- PlotDEVolcano(thy1.de.result, labels = FALSE)
print(thy1.volcano.plot)

## ----RunDiffExpression---------------------------------------------------
cluster.de.result <- RunDiffExpression(clustered.set,
                                       condition.a = "1",
                                       condition.b = "2",
                                       condition = "cluster", 
                                       fitType = "local", 
                                       method = "per-condition")
cluster.de.result[1:10,]

## ----RemoveCluster, fig.width = 5, fig.height = 4------------------------
clean.set <- SubsetCluster(clustered.set, clusters = "1")
clean.pca <- RunPCA(clean.set)
clean.cluster <- RunCORE(clean.pca, conservative = TRUE)
PlotDendrogram(clean.cluster)

## ----CleanDE, fig.width = 4, fig.height = 5------------------------------
# List of clusters to compare
cluster.list <- c("1", "2", "3", "4")

# Create a custom function to call RunDiffExpession
customFunction <- function(x, clean.cluster){
  # This is a standard RunDiffExpression call; The only difference is "x" will
  # be inputted by the sapply function
  de.result <- RunDiffExpression(clean.cluster, 
                                 condition.a = x, 
                                 condition.b = "Others",
                                 conditions = "cluster")
  # This will output the differential expression result as a list of dataframes
  return (de.result)
}

clean.cluster.de.results <- lapply(cluster.list, function(x)
  customFunction(x, clean.cluster))

# Generate volcano plots
cluster.de.1 <- PlotDEVolcano(clean.cluster.de.results[[1]], labels = FALSE)
cluster.de.2 <- PlotDEVolcano(clean.cluster.de.results[[2]], labels = FALSE)
cluster.de.3 <- PlotDEVolcano(clean.cluster.de.results[[3]], labels = FALSE)
cluster.de.4 <- PlotDEVolcano(clean.cluster.de.results[[4]], labels = FALSE)
print(cluster.de.1)
print(cluster.de.2)
print(cluster.de.3)
print(cluster.de.4)

## ----RunPairedDE---------------------------------------------------------
# Run differential expression on pairs
c1c2.de.results <- RunDiffExpression(clean.cluster, condition.a = "1", 
                                     condition.b = "2", conditions = "cluster")
c1c3.de.results <- RunDiffExpression(clean.cluster, condition.a = "1", 
                                     condition.b = "3", conditions = "cluster")
c1c4.de.results <- RunDiffExpression(clean.cluster, condition.a = "1", 
                                     condition.b = "4", conditions = "cluster")
c2c3.de.results <- RunDiffExpression(clean.cluster, condition.a = "2", 
                                     condition.b = "3", conditions = "cluster")
c2c4.de.results <- RunDiffExpression(clean.cluster, condition.a = "2", 
                                     condition.b = "4", conditions = "cluster")
c3c4.de.results <- RunDiffExpression(clean.cluster, condition.a = "3", 
                                     condition.b = "4", conditions = "cluster")

# Plot differential expression results
c1c2.plot <- PlotDEVolcano(c1c2.de.results, labels = FALSE)
c1c3.plot <- PlotDEVolcano(c1c3.de.results, labels = FALSE)
c1c4.plot <- PlotDEVolcano(c1c4.de.results, labels = FALSE)
c2c3.plot <- PlotDEVolcano(c2c3.de.results, labels = FALSE)
c2c4.plot <- PlotDEVolcano(c2c4.de.results, labels = FALSE)
c3c4.plot <- PlotDEVolcano(c3c4.de.results, labels = FALSE)

print(c1c2.plot)
print(c1c3.plot)
print(c1c4.plot)
print(c2c3.plot)
print(c2c4.plot)
print(c3c4.plot)

