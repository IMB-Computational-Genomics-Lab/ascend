#' Retinal Ganglion Cells (100 THY-1 positive, 100 THY-1 negative)
#'
#' A dataset containing a count matrix, cell-related metadata and gene-
#' related metdata.
#'
#' @format A list consisting of three objects:
#' \describe{
#'   \item{counts}{Count matrix stored in sparseMatrix format}
#'   \item{col_info}{DataFrame containing cell barcode and batch information}
#'   \item{row_info}{DataFrame containing gene names and corresponding ENSEMBL 
#'   gene identifiers}
#'   }
#' @name raw
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/}
"raw"

#' Retinal Ganglion Cells (100 THY-1 positive, 100 THY-1 negative)
#'
#' An EMSet containing the data from "raw" loaded into an EMSet.
#'
#' @format A list consisting of three objects:
#' \describe{
#'   \item{counts}{Count matrix stored in sparseMatrix format}
#'   \item{colInfo}{DataFrame containing cell barcode and batch information}
#'   \item{rowInfo}{DataFrame containing gene names and corresponding ENSEMBL 
#'   gene identifiers}
#'   \item{controls}{A list containing mitochondrial-associated transcripts and
#'   ribosomal-associated transcripts}
#'   }
#' @name raw_set
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/}
"raw_set"

#' Retinal Ganglion Cells (100 THY-1 positive, 100 THY-1 negative)
#'
#' An EMSet containing the data from "raw" loaded into an EMSet. This data has
#' been filtered and analysed as follows:
#'
#' @examples
#' \dontrun{
#' # Load package
#' library(ascend)
#' # Read dataset from disk
#' em_set <- loadCellRanger(
#' "~/Downloads/filtered_gene_bc_matrices_mex/GRCh38p7/")
#' 
#' # Subset 100 cells from each batch
#' col_info <- colInfo(em_set)
#' subset_barcodes1 <- sample(col_info$cell_barcode[which(col_info$batch == 1)], 
#' 100, replace = FALSE)
#' subset_barcodes2 <- sample(col_info$cell_barcode[which(col_info$batch == 2)], 
#' 100, replace = FALSE)
#' barcode_list <- c(subset_barcodes1, subset_barcodes2)
#' raw_set <- em_set[, barcode_list]
#' # Get raw elements to use in vignettes
#' raw_counts <- counts(raw_set)
#' raw_col_info <- colInfo(raw_set)
#' raw_row_info <- rowInfo(raw_set)
#' # Quick ascend workflow for cells
#' working_set <- normaliseBatches(raw_set)
#' batch_norm_qc <- plotBatchNormQC(raw_object = raw_set, 
#' norm_object = working_set)
#' # Batch normalisation fine.
#' qc_plots <- plotGeneralQC(working_set)
#' # Perform QC
#' working_set <- filterByOutliers(working_set, cell.threshold = 3, 
#' gene.threshold = 3, control.threshold = 3)
#' working_set <- filterLowAbundanceGenes(working_set, pct.threshold = 0.1)
#' # Check cell numbers
#' col_info <- colInfo(working_set)
#' working_set <- normaliseByRLE(working_set)
#' working_set <- excludeControl(working_set, control = c("Mt", "Rb"))
#' working_set <- runPCA(working_set, scaling = TRUE, ngenes = 100)
#' working_set <- runTSNE(working_set, PCA = TRUE, dims = 2, seed = 1)
#' working_set <- runUMAP(working_set, method = "naive")
#' working_set <- runCORE(working_set, conservative = FALSE, nres = 40, 
#' dims = 20, remove.outliers = TRUE)
#' }
#' @name analyzed_set
#' @source \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6108/}
"analyzed_set"


