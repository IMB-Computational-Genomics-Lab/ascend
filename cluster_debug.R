# Devtools
devtools::load_all(pkg = "~/CodeRepositories/ASCEND")
library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

object <- readRDS("~/Test.RDS")
object <- ReduceDimensions(object, n = 10)

clustered.obj <- RunCORE(object, conservative = TRUE)
