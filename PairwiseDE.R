devtools::load_all("~/CodeRepositories/ascend-dev")
library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

em.set <- readRDS(system.file("extdata", "RGC_scRNASeq.rds", package = "ascend"))
object <- em.set
conditions = "cluster"
condition.a = "1"
condition.b = "2"
fitType = "local"
method = "per-condition"

de.results <- RunDiffExpression(object, condition.a = condition.a, condition.b = "Others", conditions = conditions)
