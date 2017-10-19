devtools::load_all("~/CodeRepositories/ASCEND-dev/")

library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

object <- readRDS("~/Downloads/testObject.RDS")

# Clustered Object
clustered.object <- RunCORE(object, conservative = TRUE)

