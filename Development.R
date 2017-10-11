# Set up environment
devtools::load_all("/Users/a.senabouth/CodeRepositories/ASCEND")
library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar = TRUE), default = TRUE)

aem.set <- CellRangerToASCEND("/Volumes/Anne's External HD/APC_E7_scRNA/APC_E7_scRNA_Aggr_V1", "mm10")
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
aem.set <- ConvertGeneAnnotation(aem.set, "gene_symbol", "ensembl_id")
aem.set <- scranCellCycle(aem.set, mm.pairs)
aem.set <- ConvertGeneAnnotation(aem.set, "ensembl_id", "gene_symbol")

raw.qc <- PlotGeneralQC(aem.set)

filtered.set <- FilterByOutliers(aem.set)
filtered.set <- FilterByCustomControl(control.name = "Mt", pct.threshold = 20, filtered.set)
filtered.set <- FilterByCustomControl(control.name = "Rb", pct.threshold = 50, filtered.set)
filtered.set <- FilterByExpressedGenesPerCell(filtered.set, 1)

# Normalisation
rle.obj <- NormaliseByRLE(filtered.set)
scran.obj <- scranNormalise(filtered.set)
