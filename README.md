# ASCEND
## Analysis of Single Cell Expression, Normalisation and Differential expression
### Workflow Summary
![Quality control and analytical workflow with the ASCEND package](ASCEND_workflow.png)

### Installation
#### Required Packages
- **BiocParallel**
- **dplyr**
- **plyr**
- **data.table**
- **reshape2**
- **ggplot2**
- **RColorBrewer**

#### Recommended Packages
- **scater**
  - Required for the following functions:
    - scranNormalise
- **scran**
  - Required for the following functions:
    - scranNormalise
- **limSolve**
  - Required for the following functions:
    -scranNormalise (in quickCluster mode)
- **dynamicTreeCut**
  - Required for the following functions:
    - FindOptimalClusters
- **Rtsne**
  - Required for the following functions:
    - PlotTSNE
- **DESeq**
  - Required for the following functions:
    - RunClusterDiffExpression
    - RunDiffExpression
- **cellrangerRkit**
  - Required for the following functions:
    - MergeCellRanger

#### Installing ASCEND
This package is still undergoing development, so is best installed via devtools.
```R
devtools::load_all("/Path/To/ASCEND")
```

### Getting started
#### System Requirements
The amount of computational power and time required to process and analyse a single-cell RNASeq dataset is dependant on its size. Generally, ASCEND can comfortably analyse datasets that consist of up to 10,000 cells on a single machine. A single machine is defined as a desktop or laptop with 8GB of memory and four CPUs.

#### Configuring BiocParallel
This package makes extensive use of [BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html) for its functions. Before you begin, you should register a BiocParallel backend in order to get the best performance out of this package.

Here are some example configurations:

#### Single-machine multicore (Linux/Unix)
```R
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar = TRUE), default = TRUE)
```

#### Single-machine multicore (Windows, SNOW)
```R
workers <- 3
register(SnowParam(workers = workers, type = "SOCK", progressbar = TRUE), default = TRUE)
```

#### Cluster (PBSPro)
```R
example cluster code
```

### Preparing your data
#### Test Datasets
The vignette uses the public dataset ["Aggregate of t_3k and t_4k"](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.0.1/t_3k_4k_aggregate) by 10x Genomics. This dataset consists of 8093 pan T cells that are divided into two batches. The cells were isolated from a healthy donor.

##### Expression Matrix
The main source of input is an expression matrix, or a gene/cell matrix containing transcript counts. They are usually produced at the end of single cell RNA-seq processing pipelines such as Cell Ranger and DropSeq.

In an expression matrix, each row represents a gene and each column represents a cell. The names of rows and columns will subsequently be named accordingly.

ASCEND is able to use any row and column names in the expression matrix, provided they abide by the following criteria:

1. Names should not repeat. If you have a list with repeats, you can make the names unique by using R's 'make.unique' function.
2. You should be able to identify which genes you would like to select as controls. This is why gene symbols or ENSEMBL transcript IDs should be used.
3. Cells from different batches, samples or sequencing runs should be given a numeric identifier at the end. eg. BARCODE-1, BARCODE-2, BARCODE-3.

##### Combining expression matrices from different batches
You can concatenate multiple expression matrices with the function *MergeExprsMtx*. The tutorial will later explain how to normalise counts between these combined expression matrices.

If you would like to aggregate expression matrices from Cell Ranger, it is recommended you use Cell Ranger's aggr function or the *MergeCellRanger* function that comes with this package.

#### Batch Information
Batch Information is a named list containing cell identifiers and their associated batch. ASCEND will automatically generate batch information for an expression matrix if none are provided. However, it will make the assumption that there is only one batch of cells in the expression matrix.

#### Gene Annotation
The Gene Annotation slot holds a data frame that contains the gene identifiers used in the expression matrix, in addition to their corresponding identifiers in other systems. This information is optional.

#### Controls
You must provide a list of gene identifier controls. Controls are important as ....... Generally, mitochondrial and ribosomal genes are used as controls. Spike-ins are also used as controls if they were included in the study.

Controls should be organised into a named list, and identifiers used should be present in the expression matrix.

#### Loading data from Cell Ranger
ASCEND is able to load data directly from Cell Ranger with the function *CellRangerToASCEND*. This function will prepare the expression matrix, batch information, gene annotations and controls for you.

Please note that mitochondrial and ribosomal genes will be selected as controls when you use this function. If you would like to modify or replace the controls with genes of your own choosing, you may use the *EditControls* function.

Refer to the tutorial for more information on how to prepare your data for use with the ASCEND package.
