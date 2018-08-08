## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install_cran, eval = FALSE------------------------------------------
#  # List of packages to install
#  cran_packages <- c("reshape2", "fields", "ggbeeswarm", "gridExtra",
#                     "dynamicTreeCut", "dendextend", "RColorBrewer",
#                     "locfit", "KernSmooth")
#  
#  # Easy command to install all at once
#  install.packages(cran_packages)

## ----setup_bioconductor, eval = FALSE------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite()

## ----bioconductor_packages, eval = FALSE---------------------------------
#  bioconductor_packages <- c("Biobase", "BiocGenerics", "BiocParallel",
#                             "SingleCellExperiment", "GenomeInfoDb", "GenomeInfoDbData")
#  biocLite(bioconductor_packages)

## ----install_ascend------------------------------------------------------
# Load devtools package
library(devtools)

# Use devtools to install the package
install_github("IMB-Computational-Genomics-Lab/ascend", ref = "devel")

## ----load_ascend---------------------------------------------------------
# Load the package in R
library(ascend)

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

