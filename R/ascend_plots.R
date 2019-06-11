################################################################################
#
# ascend_plots.R
# description: Functions related to the plotting of data
#
################################################################################

#' plotBatchNormQC
#' 
#' Creates a boxplot that depicts the library sizes of each batch in a dataset.
#' This function is for the comparison of pre- and post- batch-normalised 
#' datasets that have been loaded into an EMSet. Cells found in these datasets 
#' must have the same identifiers.
#'
#' @param raw_object Dataset prior to batch normalisation
#' @param norm_object Dataset after batch normalisation
#' 
#' @examples
#' # Load example EMSet 
#' raw_object <- ascend::raw_set
#' 
#' # Run batch normalisation
#' norm_object <- normaliseBatches(raw_object)
#' 
#' # Generate plot
#' plot <- plotBatchNormQC(raw_object = raw_object, norm_object = norm_object)
#' 
#' @return A boxplot created with ggplot2
#' @export
plotBatchNormQC <- function(raw_object = NULL, norm_object = NULL){
  # Silence R CMD CHECK
  col_info <- "Shhh"
  type <- "Shhh"
  cell_barcode <- "Shhh"
  batch <- "Shhh"
  
  # Check if data is in correct state
  if (all(!is.null(progressLog(raw_object)$normaliseBatches) & 
          is.null(progressLog(norm_object)$normaliseBatches))){
    stop("Please ensure normaliseBatches has only been run on norm_object.")
  }
  
  # Load data from pre- and post- batch normalised objects
  exprs_counts <- SingleCellExperiment::counts(raw_object)
  exprs_normcounts <- SingleCellExperiment::counts(norm_object)
  
  # Load object information
  colInfo1 <- colInfo(raw_object)
  colInfo2 <- colInfo(norm_object)
  colData1 <- SummarizedExperiment::colData(raw_object)
  colData2 <- SummarizedExperiment::colData(norm_object)
  
  # Marry data
  objInfo1 <- S4Vectors::merge(colInfo1, colData1, by = "cell_barcode")
  objInfo2 <- S4Vectors::merge(colInfo2, colData2, by = "cell_barcode")
  
  # Take rows that are present in both datasets
  common_barcodes <- intersect(objInfo1$cell_barcode, objInfo2$cell_barcode)
  objInfo1 <- objInfo1[which(objInfo1$cell_barcode %in% common_barcodes), ]
  objInfo2 <- objInfo2[which(objInfo2$cell_barcode %in% common_barcodes), ]
  
  # Calculate library sizes
  libsize <- objInfo1[ , c("cell_barcode", "batch", "qc_libsize")]
  norm_libsize <- objInfo2[, c("cell_barcode", "batch", "qc_libsize")]
  colnames(libsize) <- c("cell_barcode", "batch", "libsize")
  colnames(norm_libsize) <- c("cell_barcode", "batch", "norm_libsize")
  
  # Build plot dataframe
  cell_df <- S4Vectors::merge(libsize, norm_libsize, 
                              by = c("cell_barcode", "batch"))
  cell_df$batch <- factor(cell_df$batch, levels = sort(unique(cell_df$batch))) 
  cell_df <- as.data.frame(cell_df)
  tidy_df <- tidyr::gather(cell_df, type, libsize, -cell_barcode, -batch)
  tidy_df$type <- factor(tidy_df$type, levels = c("libsize", "norm_libsize"))
  
  comparison_plot <- ggplot2::ggplot(tidy_df, 
                                     ggplot2::aes(x = batch, 
                                                  y = libsize, 
                                                  fill = type)) + 
    ggplot2::geom_boxplot()
  comparison_plot <- comparison_plot + 
    ggplot2::ggtitle("Library sizes per batch before and after normalisation") 
  comparison_plot <- comparison_plot + 
    ggplot2::xlab("Batch") + ggplot2::ylab("Library size") 
  comparison_plot <- comparison_plot + 
    ggplot2::scale_colour_brewer(name = "Library size", 
                                 labels = c("Raw", "Normalised"), 
                                 aesthetics = "fill")
  
  return(comparison_plot)
}

#' plotVariableGenes
#' 
#' Generates a scatter plot to aid in the detection of variable genes. Scatter
#' plot depicts Correlation of Variance (CV) vs log10(mean gene expression). 
#' Please use the \code{\link[ascend:calculateCV]{calculateCV}} 
#' function before using this function.
#' 
#' @param object An \code{\linkS4class{EMSet}} that has had CV values calculated
#' @param ngenes Select n most variable genes to plot (Optional)
#' @param label.size Size of gene labels
#' @param point.size Size of scatter points
#' @param check.overlap Hide overlapping labels (Default: FALSE)
#' @return A scatter plot rendered by ggplot2's geom_point function.
#' 
#' @examples
#' em_set <- ascend::analyzed_set
#' em_set <- calculateCV(em_set)
#' variable_gene_plot <- plotVariableGenes(em_set, ngenes = 1500, 
#' label.size = 3, point.size = 0.5, check.overlap = TRUE)
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom ggplot2 ggplot aes geom_point geom_text ggtitle xlab ylab theme_bw
#' @export
plotVariableGenes <- function(object, 
                              ngenes = NULL, 
                              label.size = 3, 
                              point.size = 0.5,
                              check.overlap = FALSE){
  # To silence checks...
  gene_mean <-"Shhh"
  gene_cv <- "Shhh"
  gene_id <- "Shhh"
  
  
  # Check CV values have been generated
  if (! "ascend_cv" %in% colnames(SummarizedExperiment::rowData(object)) ){
    stop("Please run calculateCV before using this function.")
  }
  
  # If ngenes not specified, use all values
  if (is.null(ngenes)){
    ngenes <- nrow(object)
  }
  
  # Extract rowData
  row_data <- as.data.frame(SummarizedExperiment::rowData(object))
  gene_id_name <- colnames(row_data)[1]
  plot_metrics <- row_data[ , c(gene_id_name, "qc_meancounts", 
                                "ascend_sd", "ascend_cv", "ascend_cv_rank")]
  plot_metrics[, gene_id_name] <- factor(plot_metrics[ , gene_id_name], 
                                         levels = unique(plot_metrics[ , gene_id_name]))
  colnames(plot_metrics) <- c("gene_id", "gene_mean", "gene_sd", "gene_cv", "gene_rank")
  cv_range <- round((max(plot_metrics$gene_cv)-min(plot_metrics$gene_cv))/2)
  plot_metrics$label <- plot_metrics$gene_cv >= cv_range
  
  # Subset genes by ranking
  plot_metrics <- subset(plot_metrics, plot_metrics$gene_rank <= ngenes)
  
  # Generate plot
  plot <- ggplot2::ggplot(plot_metrics, ggplot2::aes(x = log10(gene_mean), y = gene_cv)) + 
    ggplot2::geom_point(size = point.size) + 
    ggplot2::geom_text(data = subset(plot_metrics, plot_metrics$label == TRUE), 
                       ggplot2::aes(log10(gene_mean), gene_cv, label = gene_id), 
                       size = label.size, hjust = -0.25, 
                       check_overlap = check.overlap) + ggplot2::ggtitle("Variable genes") + 
    ggplot2::xlab(log[10]~"Mean gene expression") + ggplot2::ylab("Coefficiant of variation") + ggplot2::theme_bw()
  return(plot)
}

#' plotVolcano
#'
#' Produces a volcano plot featuring differential expression results.
#'
#' @param x Differential expression results generated by 
#' \code{\link{runDESeq}}.
#' @param threshold Threshold to determine significant results by (Default: 5e-3).
#' @param l2fc Threshold to determine significant Log2FoldChange by (Default: 2).
#' @param labels Display names of significant results on plot (Default: FALSE).
#' @param label.size Size of label text (Default: 5).
#' @param check.overlap If labels = TRUE, remove overlapping labels
#' @return A scatter plot rendered by ggplot2's geom_point.
#' 
#' @examples
#' # Load EMSet that has undergone analysis
#' em_set <- ascend::analyzed_set
#' 
#' # Run differential expression analysis using combined LRT
#' de_result_df <- runDiffExpression(em_set, group = "cluster", 
#' condition.a = 1, condition.b = 2, ngenes = 1500, subsampling = FALSE)
#' 
#' # Use function to generate volcano plot
#' my_volcano_plot <- plotVolcano(de_result_df, threshold = 5e-3, l2fc = 2,
#' labels = FALSE, check.overlap = FALSE)
#' 
#' @importFrom ggplot2 ggplot geom_point scale_color_manual xlab ylab theme_bw geom_text aes
#' @importFrom dplyr id
#' @export
#'
plotVolcano <- function(x, threshold = 5e-3,
                        l2fc = 2,
                        labels = FALSE,
                        label.size = 5,
                        check.overlap = TRUE){
  # Silence check
  log2FoldChange <- "Shhh"
  padj <- "Shhh"
  group <- "Shhh"
  id <- "Shhh"
  
  required_columns <- c("id", "padj", "log2FoldChange")
  if (!is.data.frame(x)){
    stop("Please supply a data frame.")
  } else{
    if (!any(required_columns %in% colnames(x))){
      stop("Please ensure your dataframe has the following columns: id, padj and log2FoldChange")
    }
  }
  
  # Convert padj to log10
  x$padj <- -log10(x$padj)
  
  # Identify significant limits
  padj_limit <- -log10(threshold)
  fc_limit1 <- l2fc
  fc_limit2 <- -l2fc
  
  # Label groups
  x$group <- "Not significant"
  sig_upregulated <- x$log2FoldChange > fc_limit1 & x$padj > padj_limit
  sig_downregulated <- x$log2FoldChange < fc_limit2 & x$padj > padj_limit
  
  if (any(sig_upregulated)){
    x[sig_upregulated, "group"] <- "Significant"
  }
  
  if (any(sig_downregulated)){
    x[sig_downregulated, "group"] <- "Significant"
  }
  
  # Remove Inf values
  x <- x[is.finite(x$log2FoldChange),]
  
  x$group <- factor(x$group, levels = c("Not significant", "Significant"))
  
  volcano_plot <- ggplot2::ggplot(x, ggplot2::aes(log2FoldChange, padj)) + 
    ggplot2::geom_point(ggplot2::aes(colour = group)) + 
    ggplot2::scale_colour_manual(name = "", values = c(
      "Not significant" = "#474747", "Significant" = "#ff5000")) + 
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + 
    ggplot2::xlab('log2 Fold Change') +
    ggplot2::ylab('-log10 adjusted P-value') 
  
  if(labels){
    if (any(sig_upregulated)){
      volcano_plot <- volcano_plot + ggplot2::geom_text(data = 
                                                          subset(x, padj > padj_limit & log2FoldChange > fc_limit1), 
                                                        ggplot2::aes(log2FoldChange, padj, label = id), 
                                                        size = label.size, vjust = 0, nudge_y = 2.5, 
                                                        check_overlap = check.overlap)
    }
    if (any(sig_downregulated)){
      volcano_plot <- volcano_plot + ggplot2::geom_text(data = subset(x, padj > padj_limit & log2FoldChange < fc_limit2), 
                                                        ggplot2::aes(log2FoldChange, padj, label = id), 
                                                        size = label.size, vjust = 0, nudge_y = 2.5, 
                                                        check_overlap = check.overlap)
    }
  }
  
  return(volcano_plot)
}

#' plotDendrogram
#' 
#' Generates a colour-labelled dendrogram to the device, using a clustered
#' \code{\linkS4class{EMSet}}.
#' 
#' @param object An \code{\linkS4class{EMSet}} that has undergone clustering.
#' @return A dendrogram plotted by base R graphics
#' 
#' @examples
#' # Load example EMSet
#' em_set <- ascend::analyzed_set
#' plotDendrogram(em_set)
#' 
#' @importFrom stats as.dendrogram order.dendrogram
#' @importFrom dendextend branches_attr_by_clusters set get_leaves_branches_col sort_levels_values colored_bars get_leaves_branches_attr
#' @importFrom graphics legend
#' @importFrom dendextend %>%
#' @export
#'
plotDendrogram <- function(object){
  # Input Checks
  if(length(clusterAnalysis(object)) == 0){
    stop("Please cluster the data in your EMSet with runCORE before using this function.")
  }
  
  # Retrieve data
  cluster_analysis <- clusterAnalysis(object)
  # Extract required values from the object
  hclust_obj <- cluster_analysis$hClust
  optimal_height <- cluster_analysis$optimalTreeHeight
  nclusters <- cluster_analysis$nClusters
  cluster_list <- cluster_analysis$clusters
  
  # Count table
  cluster_df <- as.data.frame(table(cluster_list))
  dendro_obj <- as.dendrogram(hclust_obj)
  
  # Reorder count table to match dendrogram order
  cluster_order <- order.dendrogram(dendro_obj)
  ordered_clusters <- as.vector(cluster_list)[cluster_order]
  cluster_df <- cluster_df[match(unique(ordered_clusters), cluster_df$cluster_list), ]
  
  # Generate coloured plot
  coloured_dendro <- dendextend::branches_attr_by_clusters(dendro_obj, clusters = ordered_clusters, attr = 'col')
  coloured_dendro <- dendextend::set(coloured_dendro, "labels", "")
  graphics::plot(coloured_dendro)
  
  # Generate cluster bar underneath
  dendro_colours <- unique(dendextend::get_leaves_branches_col(coloured_dendro))
  coloured_order <- cluster_order
  sorted_levels <- dendextend::sort_levels_values(as.vector(cluster_list)[coloured_order])
  sorted_levels <- sorted_levels[match(seq_along(coloured_order), coloured_order)]
  dendextend::colored_bars(dendro_colours[sorted_levels], coloured_dendro, rowLabels = "Cluster")
  
  # Get branch colours
  colour_order <- unique(coloured_dendro %>% dendextend::get_leaves_branches_attr("col"))
  
  # Generate legend labels
  legend_labels <- sapply(1:nrow(cluster_df), function(x) paste0("Cluster ", cluster_df$cluster_list[x], ": ", cluster_df$Freq[x]))
  legend("topright", legend = c(legend_labels), fill = colour_order, border = colour_order, bty = "n", title = "Cluster Populations")
}

#' plotStability
#'
#' Plots Stability, Consecutive RI and Rand Index. This can be used to determine
#' the optimal resolution of the clustering results.
#'
#' @param object An \code{\linkS4class{EMSet}} object that has undergone clustering.
#' @return A line graph generated by ggplot2's geom_line function
#' 
#' @examples
#' # Load example EMSet that has undergone processing
#' em_set <- ascend::analyzed_set
#' 
#' # Use function to plot stability scores
#' stability_plot <- plotStability(em_set)
#' 
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_line theme_bw theme element_text aes xlab ylab
#' @export
#'
plotStability <- function(object){
  # To shut R Check up...
  Height <- "Shhh"
  value <- "Shhh"
  variable <- "Shhh"
  
  if(length(clusterAnalysis(object)) == 0){
    stop("Please cluster the data in your EMSet with runCORE before using this function.")
  }
  
  # Grab info
  log <- progressLog(object)
  cluster_analysis <- clusterAnalysis(object)
  key_stats_df <- cluster_analysis$keyStats
  nres <- log$clustering$nres
  
  key_stats_df$ClusterCount <- key_stats_df$ClusterCount/10
  colnames(key_stats_df) <-c('Height', 'Stability', 'RandIndex', 'ConsecutiveRI', 'ClusterCount/10')
  key_stats_df$Height <-as.character(key_stats_df$Height)
  key_stats_df <- reshape2::melt(key_stats_df, id='Height')
  key_stats_df$Height <-as.numeric(key_stats_df$Height)
  step <- signif(1/nres, digits = 2)
  diagnostic_plot <- ggplot2::ggplot(key_stats_df)
  diagnostic_plot <- diagnostic_plot + 
    ggplot2::geom_line(ggplot2::aes(x=Height, y=value,  colour=variable)) + 
    ggplot2::theme_bw() + ggplot2::theme(axis.text = ggplot2::element_text(size=11), axis.title = ggplot2::element_text(size=11)) + 
    ggplot2::theme(legend.text = ggplot2::element_text(size=11)) + 
    ggplot2::theme(legend.title = ggplot2::element_blank()) + 
    ggplot2::xlab(sprintf('Parameter from %g to 1', step)) + ggplot2::ylab('Scores')
  return(diagnostic_plot)
}

 
PlotClusterDendro <- function (dendro, colors, groupLabels = NULL, rowText = NULL,
                               rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL,
                               textPositions = NULL, setLayout = TRUE, autoColorHeight = TRUE,
                               colorHeight = 0.2, rowWidths = NULL, dendroLabels = NULL, guideAll = FALSE, guideHang = 0.2,
                               addTextGuide = FALSE, cex.colorLabels = 0.8, cex.dendroLabels = 0.9,
                               cex.rowText = 0.8, marAll = c(1, 5, 3, 1), saveMar = TRUE,
                               abHeight = NULL, abCol = "red", ...)
{
  dendro$labels <- rep('', length(dendro$labels))
  
  oldMar = graphics::par("mar")
  if (!is.null(dim(colors))) {
    nRows = dim(colors)[2]
  }
  else nRows = 1
  if (!is.null(rowText))
    nRows = nRows + if (is.null(textPositions))
      nRows
  else length(textPositions)
  if (autoColorHeight)
    colorHeight = 0.2 + 0.3 * (1 - exp(-(nRows - 1)/6))
  if (setLayout)
    graphics::layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight,
                                             colorHeight))
  graphics::par(mar = c(0, marAll[2], marAll[3], marAll[4]))
  graphics::plot(dendro, labels = dendroLabels, cex = cex.dendroLabels,
       ...)
  if (!is.null(abHeight))
    graphics::abline(h = abHeight, col = abCol)
  graphics::par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
  PlotConsensusBars(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels,
                    rowText = rowText, rowTextAlignment = rowTextAlignment,
                    rowTextIgnore = rowTextIgnore, textPositions = textPositions,
                    cex.rowText = cex.rowText, rowWidths = rowWidths, addTextGuide = addTextGuide)
  if (saveMar)
    graphics::par(mar = oldMar)
}


PlotConsensusBars <- function (dendro, colors, rowLabels = NULL, rowWidths = NULL,
                               rowText = NULL, rowTextAlignment = c("left", "center", "right"),
                               rowTextIgnore = NULL, textPositions = NULL, addTextGuide = TRUE,
                               cex.rowLabels = 1, cex.rowText = 0.8, ...) {
  PlotOrderedColors(dendro$order, colors = colors, rowLabels = rowLabels,
                    rowWidths = rowWidths, rowText = rowText, rowTextAlignment = rowTextAlignment,
                    rowTextIgnore = rowTextIgnore, textPositions = textPositions,
                    addTextGuide = addTextGuide, cex.rowLabels = cex.rowLabels,
                    cex.rowText = cex.rowText, startAt = 0, ...)
}

PlotOrderedColors <- function (order, colors, rowLabels = NULL, rowWidths = NULL,
                               rowText = NULL, rowTextAlignment = c("left", "center", "right"),
                               rowTextIgnore = NULL, textPositions = NULL, addTextGuide = TRUE,
                               cex.rowLabels = 1, cex.rowText = 0.8, startAt = 0, ...) {
  colors = as.matrix(colors)
  dimC = dim(colors)
  
  # Create a colour ramp
  gradient_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  grDevices::palette(gradient_palette(max(colors)))
  
  if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) ==
                            dimC[2]))
    rowLabels = colnames(colors)
  sAF = options("stringsAsFactors")
  options(stringsAsFactors = FALSE)
  on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)
  nColorRows = dimC[2]
  if (length(order) != dimC[1])
    stop("ERROR: length of colors vector not compatible with number of objects in 'order'.")
  C = colors[order, , drop = FALSE]
  step = 1/(dimC[1] - 1 + 2 * startAt)
  graphics::barplot(height = 1, col = "white", border = FALSE, space = 0,
                    axes = FALSE)
  charWidth = graphics::strwidth("W")/2
  if (!is.null(rowText)) {
    if (is.null(textPositions))
      textPositions = c(1:nColorRows)
    if (is.logical(textPositions))
      textPositions = c(1:nColorRows)[textPositions]
    nTextRows = length(textPositions)
  }
  else nTextRows = 0
  nRows = nColorRows + nTextRows
  ystep = 1/nRows
  if (is.null(rowWidths)) {
    rowWidths = rep(ystep, nColorRows + nTextRows)
  }
  else {
    if (length(rowWidths) != nRows)
      stop("plotOrderedColors: Length of 'rowWidths' must equal the total number of rows.")
    rowWidths = rowWidths/sum(rowWidths)
  }
  hasText = rep(0, nColorRows)
  hasText[textPositions] = 1
  csPosition = cumsum(c(0, hasText[-nColorRows]))
  colorRows = c(1:nColorRows) + csPosition
  rowType = rep(2, nRows)
  rowType[colorRows] = 1
  physicalTextRow = c(1:nRows)[rowType == 2]
  yBottom = c(0, cumsum(rowWidths[nRows:1]))
  yTop = cumsum(rowWidths[nRows:1])
  if (!is.null(rowText)) {
    rowTextAlignment = match.arg(rowTextAlignment)
    rowText = as.matrix(rowText)
    textPos = list()
    textPosY = list()
    textLevs = list()
    for (tr in 1:nTextRows) {
      charHeight = max(graphics::strheight(rowText[, tr], cex = cex.rowText))
      width1 = rowWidths[physicalTextRow[tr]]
      nCharFit = floor(width1/charHeight/1.7/graphics::par("lheight"))
      if (nCharFit < 1)
        stop("Rows are too narrow to fit text. Consider decreasing cex.rowText.")
      set = textPositions[tr]
      textLevs[[tr]] = sort(unique(rowText[, tr]))
      textLevs[[tr]] = textLevs[[tr]][!textLevs[[tr]] %in%
                                        rowTextIgnore]
      nLevs = length(textLevs[[tr]])
      textPos[[tr]] = rep(0, nLevs)
      orderedText = rowText[order, tr]
      for (cl in 1:nLevs) {
        ind = orderedText == textLevs[[tr]][cl]
        sind = ind[-1]
        ind1 = ind[-length(ind)]
        starts = c(if (ind[1]) 1 else NULL, which(!ind1 &
                                                    sind) + 1)
        ends = which(c(ind1 & !sind, ind[length(ind)]))
        if (length(starts) == 0)
          starts = 1
        if (length(ends) == 0)
          ends = length(ind)
        if (ends[1] < starts[1])
          starts = c(1, starts)
        if (ends[length(ends)] < starts[length(starts)])
          ends = c(ends, length(ind))
        lengths = ends - starts
        long = which.max(lengths)
        textPos[[tr]][cl] = switch(rowTextAlignment,
                                   left = starts[long], center = (starts[long] +
                                                                    ends[long])/2 + 0.5, right = ends[long] +
                                     1)
      }
      if (rowTextAlignment == "left") {
        yPos = seq(from = 1, to = nCharFit, by = 1)/(nCharFit +
                                                       1)
      }
      else {
        yPos = seq(from = nCharFit, to = 1, by = -1)/(nCharFit +
                                                        1)
      }
      textPosY[[tr]] = rep(yPos, ceiling(nLevs/nCharFit) +
                             5)[1:nLevs][rank(textPos[[tr]])]
    }
  }
  jIndex = nRows
  if (is.null(rowLabels))
    rowLabels = c(1:nColorRows)
  C[is.na(C)] = "grey"
  for (j in 1:nColorRows) {
    jj = jIndex
    ind = (1:dimC[1])
    xl = (ind - 1.5 + startAt) * step
    xr = (ind - 0.5 + startAt) * step
    yb = rep(yBottom[jj], dimC[1])
    yt = rep(yTop[jj], dimC[1])
    if (is.null(dim(C))) {
      graphics::rect(xl, yb, xr, yt, col = as.character(C), border = as.character(C))
    }
    else {
      graphics::rect(xl, yb, xr, yt, col = as.character(C[, j]),
           border = as.character(C[, j]))
    }
    graphics::text(rowLabels[j], pos = 2, x = -charWidth/2 + xl[1],
         y = (yBottom[jj] + yTop[jj])/2, cex = cex.rowLabels,
         xpd = TRUE)
    textRow = match(j, textPositions)
    if (is.finite(textRow)) {
      jIndex = jIndex - 1
      xt = (textPos[[textRow]] - 1.5) * step
      xt[xt < graphics::par("usr")[1]] = graphics::par("usr")[1]
      xt[xt > graphics::par("usr")[2]] = graphics::par("usr")[2]
      yt = yBottom[jIndex] + (yTop[jIndex] - yBottom[jIndex]) *
        (textPosY[[textRow]] + 1/(2 * nCharFit + 2))
      nt = length(textLevs[[textRow]])
      if (addTextGuide)
        for (l in 1:nt) graphics::lines(c(xt[l], xt[l]), c(yt[l],
                                                 yTop[jIndex]), col = "darkgrey", lty = 3)
      textAdj = c(0, 0.5, 1)[match(rowTextAlignment, c("left",
                                                       "center", "right"))]
      graphics::text(textLevs[[textRow]], x = xt, y = yt, adj = c(textAdj,
                                                        1), xpd = TRUE, cex = cex.rowText)
    }
    jIndex = jIndex - 1
  }
  for (j in 0:(nColorRows + nTextRows)) graphics::lines(x = c(0, 1),
                                              y = c(yBottom[j + 1], yBottom[j + 1]))
}

#' plotStabilityDendro
#'
#' Generate a dendrogram with coloured bars representing each clustering result
#' below. This function is derived from the plotDendroAndColors function found
#' in the \pkg{WGCNA} package.
#' 
#' @param object An \code{\linkS4class{EMSet}} that has undergone clustering.
#' @return A dendrogram with coloured bars generated by R base graphics
#' 
#' @examples
#' # Load clustered EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Plot stability dendrogram
#' plotStabilityDendro(em_set)
#' 
#' @export
#'
plotStabilityDendro <- function(object){
  # Check that the user has done the required steps.
  if (length(clusterAnalysis(object)) == 0){
    stop("Please run RunCORE on this object before using this function.")
  }
  
  # Get the variables
  cluster_analysis <- clusterAnalysis(object)
  dendro <- cluster_analysis$hClust
  colours <- cluster_analysis$clusteringMatrix
  colnames(colours) <- c(1:40, "REF")
  # Plotting function
  print(PlotClusterDendro(dendro, colours))
}


#' plotMDS
#'
#' Generates a 2D MDS plot that can use the following inputs stored in an
#' EMSet:
#' 
#' 1. Distance matrix generated by \code{\link[ascend:runCORE]{runCORE}}.
#' 2. PCA matrix generated by \code{\link[ascend:runPCA]{runPCA}}.
#' 3. Normalised, log-transformed count matrix stored in logcounts.
#'   
#' @param object An \code{\linkS4class{EMSet}}.
#' @param Dim1 Dimension to plot on the x-axis.
#' @param Dim2 Dimension to plot on the y-axis.
#' @param PCA Use pre-calculated PCA-reduced values (Default: TRUE)
#' @param group (Optional) Name of the column in colInfo that 
#' describe a set of conditions you would like to colour cells by.
#' @param density Vary alpha by density (Default: FALSE)
#' @return A ggplot glob that contains a scatter plot.
#'
#' @examples
#' # Load EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Generate MDS
#' mds_plot <- plotMDS(em_set, Dim1 = 1, Dim2 = 2, 
#' group = "cluster", PCA = TRUE)
#' 
#' @importFrom stats dist cmdscale
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom ggplot2 ggplot aes geom_point labs ggtitle scale_alpha theme_bw
#' @export
#'
plotMDS <- function(object, Dim1 = 1, Dim2 = 2, group = NULL, density = FALSE, PCA = TRUE){
  # Load data from EMSet
  col_info <- colInfo(object)
  
  # Check the column has been defined
  if (!is.null(group)){
    if (!(group %in% colnames(col_info))){
      stop("Please ensure your specified column exists.")
    } else{
      condition_list <- col_info[, group]
      names(condition_list) <- col_info$cell_barcode
    }
  } else{
    condition_list <- NULL
  }
  
  # If clustering has been run...
  if (length(clusterAnalysis(object)) > 0){
    print("EMSet has undergone clustering. Retrieving distance matrix...")
    cluster_analysis <- clusterAnalysis(object)
    distance_matrix <- cluster_analysis$distanceMatrix
  } else{
    # If user has specified PCA
    if(PCA){
      # Check if the user has run PCA
      if (!("PCA" %in% SingleCellExperiment::reducedDimNames(object))){
        stop("Please use runPCA on this object before using the PCA argument for this function.")
      } else{
        print("Using PCA-reduced matrix to generate distance matrix...")
        pca_matrix <- SingleCellExperiment::reducedDim(object, "PCA")
        pca_matrix <- as.matrix(pca_matrix[ ,1:20])
        distance_matrix <- stats::dist(pca_matrix[ , 1:20])
      }
    } else{
      print("Using expression matrix to generate distance matrix...")
      expression_matrix <- SingleCellExperiment::logcounts(object)
      
      # Identify most variable genes to decrease scope of distance matrix
      gene_variance <- calcVariance(expression_matrix, axis = "row")
      names(gene_variance) <- rownames(expression_matrix)
      sorted_gene_variance <- gene_variance[order(unlist(gene_variance), decreasing = TRUE)]
      top_genes <- sorted_gene_variance[1:1000]
      
      expression_matrix <- expression_matrix[names(top_genes), ]
      expression_matrix <- Matrix::t(expression_matrix)
      distance_matrix <- stats::dist(expression_matrix)
    }  
  }
  
  # We finally have a distance matrix
  print("Running cmdscale...")
  mds_matrix <-stats::cmdscale(distance_matrix, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  
  print("Cmdscale complete! Processing scaled data...")
  mds_matrix <- as.data.frame(as.matrix(mds_matrix))
  colnames(mds_matrix) <- c("Dim1", "Dim2")
  
  if (density){
    mds_matrix$density <- fields::interp.surface(MASS::kde2d(mds_matrix$Dim1, mds_matrix$Dim2), mds_matrix[, c("Dim1", "Dim2")])
    
    print("Generating MDS plot...")
    if(!is.null(condition_list) && length(condition_list) > 0){
      mds_matrix <- as.data.frame(cbind(mds_matrix, group = factor(condition_list, levels = unique(condition_list))))
      mds_plot <- ggplot2::ggplot(mds_matrix, ggplot2::aes(Dim1, Dim2, alpha = 1/density)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::scale_alpha(range = c(0.3, 1), guide ="none") + 
        ggplot2::ggtitle("Multi-Dimensional Scaling (MDS) Plot")
    } else{
      mds_plot <- ggplot2::ggplot(mds_matrix, ggplot2::aes(Dim1, Dim2, alpha = 1/density)) + 
        ggplot2::geom_point() + ggplot2::theme_bw() + 
        ggplot2::scale_alpha(range = c(0.3, 1), guide ="none") + 
        ggplot2::ggtitle("Multi-Dimensional Scaling (MDS) Plot")
    }
  } else{
    print("Generating MDS plot...")
    if(!is.null(condition_list) && length(condition_list) > 0){
      mds_matrix <- as.data.frame(cbind(mds_matrix, group = factor(condition_list, levels = unique(condition_list))))
      mds_plot <- ggplot2::ggplot(mds_matrix, ggplot2::aes(Dim1, Dim2)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::ggtitle("Multi-Dimensional Scaling (MDS) Plot")
    } else{
      mds_plot <- ggplot2::ggplot(mds_matrix, ggplot2::aes(Dim1, Dim2)) + 
        ggplot2::geom_point() + ggplot2::theme_bw() + 
        ggplot2::ggtitle("Multi-Dimensional Scaling (MDS) Plot")
    }
  }
  return(mds_plot)
}

#' plotUMAP
#'
#' Generates a 2D UMAP plot using a TSNE matrix generated by the 
#' \code{\link[ascend:runUMAP]{runUMAP}} function.
#' 
#' @param object An \code{\linkS4class{EMSet}} object that hasv values 
#' calculated by UMAP in the reducedDim slot.
#' @param Dim1 Dimension to plot on the x-axis.
#' @param Dim2 Dimension to plot on the y-axis.
#' @param group (Optional) Name of the column in colInfo that 
#' describe a set of conditions you would like to colour cells by.
#' @param density Vary alpha by density (Default: FALSE)
#' @return A ggplot glob that contains a scatter plot.
#' @examples
#' # Load pre-processed data
#' em_set <- ascend::analyzed_set
#' 
#' # Plot with t-SNE 
#' umap_plot <- plotUMAP(em_set, Dim1 = 1, Dim2 = 2, 
#' group = "cluster", density = FALSE)
#' 
#' @importFrom ggplot2 ggplot aes geom_point labs ggtitle scale_alpha theme_bw
#' @export
#'
plotUMAP <- function(object, Dim1 = 1, Dim2 = 2, group = NULL, density = FALSE){
  if(!("UMAP" %in% SingleCellExperiment::reducedDimNames(object))){
    stop("Please generate UMAP coordinates with runUMAP.")
  }
  
  # Load data from EMSet
  col_info <- colInfo(object)
  
  # Check the column has been defined
  if (!is.null(group)){
    if (!(group %in% colnames(col_info))){
      stop("Please ensure your specified column exists.")
    } else{
      condition_list <- col_info[, group]
      names(condition_list) <- col_info$cell_barcode
    }
  } else{
    condition_list <- NULL
  }
  
  # Extract dimensions to plot
  pca_matrix <- reducedDim(object, "UMAP")
  plot_df <- as.data.frame(as.matrix(pca_matrix[, c(Dim1, Dim2)]))
  colnames(plot_df) <- c("Dim1", "Dim2")
  plot_df < as.data.frame(plot_df)
  
  if (density){
    # Adjust alpha by adding bivariate density for each point
    plot_df$density <- fields::interp.surface(MASS::kde2d(plot_df$Dim1, plot_df$Dim2), plot_df[, c("Dim1", "Dim2")])
    
    if (!is.null(condition_list) && length(condition_list) > 0){
      plot_df <- cbind(plot_df, group = condition_list)
      plot_df <- as.data.frame(plot_df)
      umap_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2, alpha = 1/density)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::scale_alpha(range = c(0.4, 1), guide ="none") + 
        ggplot2::ggtitle("Uniform Manifold Approximation and Projection (UMAP) Plot")
    } else{
      umap_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2, alpha = 1/density)) + 
        ggplot2::geom_point() +
        ggplot2::scale_alpha(range = c(0.4, 1), guide ="none") + 
        ggplot2::ggtitle("Uniform Manifold Approximation and Projection (UMAP) Plot")
    }
  } else{
    if (!is.null(condition_list) && length(condition_list) > 0){
      plot_df <- cbind(plot_df, group = condition_list)
      plot_df <- as.data.frame(plot_df)
      umap_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::ggtitle("Uniform Manifold Approximation and Projection (UMAP) Plot")
    } else{
      umap_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2)) + 
        ggplot2::geom_point() +
        ggplot2::ggtitle("Uniform Manifold Approximation and Projection (UMAP) Plot")
    }
  }
  return(umap_plot)
}



#' plotTSNE
#'
#' Generates a 2D TSNE plot using a TSNE matrix generated by the 
#' \code{\link[ascend:runTSNE]{runTSNE}} function.
#' 
#' @param object An \code{\linkS4class{EMSet}} object that has undergone PCA.
#' @param Dim1 Dimension to plot on the x-axis.
#' @param Dim2 Dimension to plot on the y-axis.
#' @param group (Optional) Name of the column in colInfo that 
#' describe a set of conditions you would like to colour cells by.
#' @param density Vary alpha by density (Default: FALSE)
#' @return A ggplot glob that contains a scatter plot.
#' @examples
#' # Load analyzed EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Use function to generate plots
#' tsne_plot <- plotTSNE(em_set, Dim1 = 1, Dim2 = 2, group = "cluster")
#' 
#' @importFrom ggplot2 ggplot aes geom_point labs ggtitle scale_alpha theme_bw
#' @export
#'
plotTSNE <- function(object, Dim1 = 1, Dim2 = 2, group = NULL, density = FALSE){
  if(!("TSNE" %in% SingleCellExperiment::reducedDimNames(object))){
    stop("Please supply an object that has undergone t-SNE.")
  }
  
  # Load data from EMSet
  col_info <- colInfo(object)
  
  # Check the column has been defined
  if (!is.null(group)){
    if (!(group %in% colnames(col_info))){
      stop("Please ensure your specified column exists.")
    } else{
      condition_list <- col_info[, group]
      names(condition_list) <- col_info$cell_barcode
    }
  } else{
    condition_list <- NULL
  }
  
  # Extract dimensions to plot
  pca_matrix <- reducedDim(object, "TSNE")
  plot_df <- as.data.frame(as.matrix(pca_matrix[, c(Dim1, Dim2)]))
  colnames(plot_df) <- c("Dim1", "Dim2")
  plot_df < as.data.frame(plot_df)
  
  if (density){
    # Adjust alpha by adding bivariate density for each point
    plot_df$density <- fields::interp.surface(MASS::kde2d(plot_df$Dim1, plot_df$Dim2), plot_df[, c("Dim1", "Dim2")])
    
    if (!is.null(condition_list) && length(condition_list) > 0){
      plot_df <- cbind(plot_df, group = condition_list)
      plot_df <- as.data.frame(plot_df)
      tsne_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2, alpha = 1/density)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::scale_alpha(range = c(0.4, 1), guide ="none") + 
        ggplot2::ggtitle("t-Distributed Stochastic Neighbor Embedding (t-SNE) Plot")
    } else{
      tsne_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2, alpha = 1/density)) + 
        ggplot2::geom_point() +
        ggplot2::scale_alpha(range = c(0.4, 1), guide ="none") + 
        ggplot2::ggtitle("t-Distributed Stochastic Neighbor Embedding (t-SNE) Plot")
    }
  } else{
    if (!is.null(condition_list) && length(condition_list) > 0){
      plot_df <- cbind(plot_df, group = condition_list)
      plot_df <- as.data.frame(plot_df)
      tsne_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::ggtitle("t-Distributed Stochastic Neighbor Embedding (t-SNE) Plot")
    } else{
      tsne_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(Dim1, Dim2)) + 
        ggplot2::geom_point() +
        ggplot2::ggtitle("t-Distributed Stochastic Neighbor Embedding (t-SNE) Plot")
    }
  }
  return(tsne_plot)
}


#' plotPCA
#'
#' Plot two principal components (PCs) on a scatter plot. This plot corresponds
#' more closely to the distance between points and therefore is good to use
#' to review the effectiveness of clustering by the CORE algorithm.
#'
#' @param object An \code{\linkS4class{EMSet}} object that has undergone PCA.
#' @param PCX Principal component to plot on the x-axis.
#' @param PCY Principal component to plot on the y-axis.
#' @param group (Optional) Name of the column in colInfo that 
#' describe a set of conditions you would like to colour cells by.
#' @param density Vary alpha by density (Default: FALSE)
#' @return A ggplot glob that contains a scatter plot.
#' 
#' @examples
#' # Load EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Generate PCA plot
#' pca_plot <- plotPCA(em_set, PCX = 1, PCY = 2, 
#' group = "cluster", density = TRUE)
#' 
#' @importFrom ggplot2 ggplot aes geom_point labs ggtitle scale_alpha theme_bw
#' @importFrom fields interp.surface
#' @importFrom MASS kde2d
#' @export
#'
plotPCA <- function(object, PCX = 1, PCY = 2, group = NULL, density = FALSE){
  if(!("PCA" %in% SingleCellExperiment::reducedDimNames(object))){
    stop("Please supply an object that has undergone PCA reduction.")
  }
  
  # Load data from EMSet
  col_info <- colInfo(object)
  
  # Check the column has been defined
  if (!is.null(group)){
    if (!(group %in% colnames(col_info))){
      stop("Please ensure your specified column exists.")
    } else{
      condition_list <- col_info[, group]
      names(condition_list) <- col_info$cell_barcode
    }
  } else{
    condition_list <- NULL
  }
  
  # Extract dimensions to plot
  pca_matrix <- SingleCellExperiment::reducedDim(object, "PCA")
  plot_df <- pca_matrix[, c(PCX, PCY)]
  plot_df <- as.data.frame(as.matrix(plot_df))
  colnames(plot_df) <- c("PCX", "PCY")
  plot_df < as.data.frame(plot_df)
  
  if (density){
    # Adjust alpha by adding bivariate density for each point
    plot_df$density <- fields::interp.surface(MASS::kde2d(plot_df$PCX, plot_df$PCY), plot_df[, c("PCX", "PCY")])
    
    if (!is.null(condition_list) && length(condition_list) > 0){
      plot_df <- cbind(plot_df, group = condition_list)
      plot_df <- as.data.frame(plot_df)
      pca_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(PCX, PCY, alpha = 1/density)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::scale_alpha(range = c(0.4, 1), guide ="none") + 
        ggplot2::ggtitle("Principal Component Analysis (PCA) Plot")
    } else{
      pca_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(PCX, PCY, alpha = 1/density)) + 
        ggplot2::geom_point() +
        ggplot2::theme_bw() + ggplot2::scale_alpha(range = c(0.4, 1), guide ="none") +
        ggplot2::ggtitle("Principal Component Analysis (PCA) Plot")
    }    
  } else{
    if (!is.null(condition_list) && length(condition_list) > 0){
      plot_df <- cbind(plot_df, group = condition_list)
      plot_df <- as.data.frame(plot_df)
      pca_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(PCX, PCY)) + 
        ggplot2::geom_point(ggplot2::aes(colour = factor(group))) + 
        ggplot2::labs(colour = group) + ggplot2::theme_bw() + 
        ggplot2::ggtitle("Principal Component Analysis (PCA) Plot")
    } else{
      pca_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(PCX, PCY)) + 
        ggplot2::geom_point() + ggplot2::theme_bw() +
        ggplot2::ggtitle("Principal Component Analysis (PCA) Plot")
    }    
  }
  return(pca_plot)
}

#' plotPCAVariance
#'
#' Generates a scree plot of PCs generated by 
#' \code{\link[ascend:runPCA]{runPCA}}. 
#' The principal components (PCs) are sorted from largest percent variance to 
#' the smallest percent variance. This plot is for determining the optimal 
#' number of PCs to retain for further analysis.
#'
#' @param object An \code{\linkS4class{EMSet}} that has undergone PCA.
#' @param n Number of PCs to view on the plot.
#' @return A ggplot2 glob that contains a scatter qplot.
#'
#' @examples
#' # Load EMSet
#' em_set <- ascend::analyzed_set
#' 
#' # Plot PCA variance
#' pca_variance_plot <- plotPCAVariance(em_set, n = 10)
#' 
#' @importFrom ggplot2 aes geom_point geom_segment ggtitle xlab ylab
#' @export
#'
plotPCAVariance <- function(object, n = 50){
  PC <- "Shhh"
  Variance <- "Shhh"
  
  # Generate scree plot
  if(is.null(progressLog(object)$PCAVariance)){
    stop("Please supply an object that has undergone PCA reduction.")
  }
  
  if (n > ncol(SingleCellExperiment::reducedDim(object, "PCA"))){
    n <- ncol(SingleCellExperiment::reducedDim(object, "PCA"))
  }
  
  # Data frame
  pca_df <- data.frame(PC = 1:length(progressLog(object)$PCAVariance), 
                       Variance = progressLog(object)$PCAVariance)
  pca_df <- pca_df[1:n, ]
  pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC, y = Variance)) + 
    ggplot2::geom_point(size = 3) + ggplot2::geom_segment(ggplot2::aes(x = PC,
                                                                       xend = PC,
                                                                       y = 0,
                                                                       yend = Variance))
  pca_plot <- pca_plot + ggplot2::ggtitle("Percent Variance per PC")
  pca_plot <- pca_plot + ggplot2::xlab("Principal Component (PC)")
  pca_plot <- pca_plot + ggplot2::ylab("Variance") + ggplot2::theme_bw()
  
  return(pca_plot)
}


#' plotNormQC
#'
#' Generates a series of plots comparing an un-normalised and a normalised 
#' \code{\linkS4class{EMSet}}.
#'
#' @param object A normalised \code{\linkS4class{EMSet}}.
#' @param gene_list OPTIONAL: A list of genes to plot expression levels for. 
#' If not defined, \pkg{ascend} will select GAPDH and MALAT1, or choose two 
#' genes at random.
#' 
#' @return A list containing the following plots:
#' \itemize{
#' \item{Library size histograms.}
#' \item{Scatter boxplots for selected genes.}
#' \item{Scatter boxplots for each gene.}
#' }
#' 
#' @examples
#' # Load normalised EMSet
#' EMSet <- ascend::analyzed_set
#' 
#' # Plot normalisation
#' norm_qc <- plotNormQC(EMSet, gene_list = c("GAPDH", "MALAT1"))
#' 
#' @importFrom ggplot2 aes geom_histogram xlab ylab ggtitle ggplot_build qplot geom_boxplot scale_y_continuous labs theme element_blank
#' @importFrom tidyr gather
#' @export
#'
plotNormQC <- function(object, gene_list = list()){
  # Get check to shut up
  libsize_count <- "Shhh"
  libsize_normcount <- "Shhh"
  gene <- "Shhh"
  count <- "Shhh"
  cell_barcode <- "Shhh"
  gene_id <- "Shhhh"
  
  
  # Check if the data has been normalised.
  if (is.null(progressLog(object)$NormalisationMethod)){
    stop("Please specify a normalised EMSet.")
  }
  
  # Create a list to output results to
  output_list <- list()
  
  # Get information from EMSet
  col_info <- colInfo(object)
  col_data <- SummarizedExperiment::colData(object)
  cell_info <- S4Vectors::merge(col_info, col_data, by = "cell_barcode")
  count_matrix <- SingleCellExperiment::counts(object)
  norm_matrix <- SingleCellExperiment::normcounts(object)
  
  # 1. Libsize histograms
  qc_normcount <- Matrix::colSums(norm_matrix[, cell_info$cell_barcode])
  wide_df <- data.frame(cell_barcode = cell_info$cell_barcode, libsize_count = cell_info$qc_libsize, libsize_normcount = qc_normcount)
  wide_df$cell_barcode <- factor(wide_df$cell_barcode, levels = wide_df$cell_barcode)
  
  min_lim <- min(min(wide_df$libsize_count), min(wide_df$libsize_normcount))
  max_lim <- max(max(wide_df$libsize_count), max(wide_df$libsize_normcount))
  breaks <- seq(min_lim, max_lim, by = 10^(min(round(abs(log10(min_lim))), round(abs(log10(max_lim))))))
  
  # Build un-normalised histogram
  libsize_count_hist <- ggplot2::ggplot(wide_df, ggplot2::aes(libsize_count)) +
    ggplot2::geom_histogram(breaks = breaks, fill = "#000000")
  
  # Build normalised histogram
  libsize_normcount_hist <- ggplot2::ggplot(wide_df, ggplot2::aes(libsize_normcount)) +
    ggplot2::geom_histogram(breaks = breaks, fill = "#000000")
  
  # Find common max
  count_max <- max(ggplot2::ggplot_build(libsize_count_hist)$data[[1]]$count)
  normcount_max <- max(ggplot2::ggplot_build(libsize_normcount_hist)$data[[1]]$count)
  y_max <- max(count_max, normcount_max)
  
  # Add to plots
  libsize_count_hist <- libsize_count_hist + ggplot2::scale_y_continuous(limits = c(0, y_max*1.1))
  libsize_normcount_hist <- libsize_normcount_hist + ggplot2::scale_y_continuous(limits = c(0, y_max*1.1))
  
  # Add labels
  libsize_count_hist <- libsize_count_hist + ggplot2::ggtitle("Pre-normalised library sizes")
  libsize_count_hist <- libsize_count_hist + ggplot2::xlab("Library size") + ggplot2::ylab("Number of cells") + ggplot2::theme_bw()
  
  libsize_normcount_hist <- libsize_normcount_hist + ggplot2::ggtitle("Normalised library sizes", 
                                                                      subtitle = sprintf("Normalised via %s", progressLog(object)$NormalisationMethod))
  libsize_normcount_hist <- libsize_normcount_hist + ggplot2::xlab("Library size") + ggplot2::ylab("Number of cells") + ggplot2::theme_bw()
  output_list[["libsize_histograms"]] <- list(count = libsize_count_hist, normcount = libsize_normcount_hist)
  
  # Plot scatters of specific gene
  # Validate user-supplied genes
  if (length(gene_list) > 0){
    # Check if gene list is in the matrix
    if (any(gene_list %in% rownames(object))){
      gene_list <- gene_list[which(gene_list %in% rownames(object))]
    }
  }
  
  # If users haven't supplied one, use house-keeping genes instead
  if (length(gene_list) == 0){
    gene_list <- c(grep("GAPDH", rownames(object), ignore.case = TRUE, value = TRUE), grep("MALAT1", rownames(object), ignore.case = TRUE, value = TRUE))
    
    # If house-keeping genes aren't present, choose genes by random
    if (length(gene_list) == 1){
      gene_list <- gene_list[which(!(sapply(gene_list, is.null)))]
    } else if (length(gene_list) == 0){
      # Sample random counts from transcrips that are in at least 50% of cells
      row_data <- rowData(object)
      gene_list <- sample(row_data[which(row_data$qc_ncells  > ncol(object) * 0.5), 1], 2)
    }
  }
  
  # Extract counts to decrease load on plotting function
  extractCounts <- function(gene, count_matrix = NULL, norm_matrix = NULL){
    gene_count_df <- data.frame(counts = count_matrix[gene, ], normcounts = norm_matrix[gene, ])
    return(gene_count_df)
  }
  
    
  plotNormGeneCounts <- function(gene, gene_counts = NULL, normalisation_method = NULL){
    cell_barcode <- "Shhhh"
    
    print(sprintf("Plotting %s expression...", gene))
    gene_count_df <- gene_counts[[gene]]
    
    ylim_max <- max(c(gene_count_df$counts, gene_count_df$normcounts))
    gene_count_df$cell_barcode <- factor(rownames(gene_count_df), levels = rownames(gene_count_df))
    scatter_count <- ggplot2::ggplot(gene_count_df, ggplot2::aes(x = cell_barcode, y = counts)) + ggplot2::geom_point() + ggplot2::ylim(c(0, ylim_max*1.1))
    scatter_normcount <- ggplot2::ggplot(gene_count_df, ggplot2::aes(x = cell_barcode, y = normcounts)) + ggplot2::geom_point() + ggplot2::ylim(c(0, ylim_max*1.1))
    scatter_count <- scatter_count + ggplot2::ggtitle(sprintf("Pre-normalised counts for %s", gene)) + ggplot2::xlab("Cell") + ggplot2::ylab("Count") + ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank())
    scatter_normcount <- scatter_normcount + ggplot2::ggtitle(sprintf("Normalised counts for %s", gene), subtitle = sprintf("Normalised via %s", normalisation_method)) + ggplot2::xlab("Cell") + ggplot2::ylab("Normalised count") + ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank())
    return(list(counts = scatter_count, normcounts = scatter_normcount))
  }
  
  # Extract gene counts
  gene_counts <- BiocParallel::bplapply(gene_list, extractCounts, count_matrix = count_matrix, norm_matrix = norm_matrix)
  names(gene_counts) <- gene_list
  
  # Generate plots
  normalisation_method <- progressLog(object)$NormalisationMethod
  norm_genecounts <- BiocParallel::bplapply(gene_list, plotNormGeneCounts, gene_counts = gene_counts, normalisation_method = normalisation_method)
  names(norm_genecounts) <- gene_list
  output_list[["sampled_genes"]] <- norm_genecounts
  
  print("Plotting gene expression box plots...")
  ordered_counts <- count_matrix[order(Matrix::rowSums(count_matrix), decreasing=TRUE),]
  ordered_normcounts <- norm_matrix[order(Matrix::rowSums(norm_matrix), decreasing=TRUE),]
  
  if (is(ordered_counts, "sparseMatrix")){
    ordered_counts <- as.matrix(ordered_counts)
  }
  if (is(ordered_normcounts, "sparseMatrix")){
    ordered_normcounts <- as.matrix(ordered_normcounts)
  }
  
  ordered_counts <- as.data.frame(ordered_counts)
  ordered_normcounts <- as.data.frame(ordered_normcounts)
  ordered_counts <- cbind(gene_id = rownames(ordered_counts), ordered_counts)
  ordered_normcounts <- cbind(gene_id = rownames(ordered_normcounts), ordered_normcounts)
  
  ordered_counts <- tidyr::gather(ordered_counts[,1:101], gene, count, -gene_id)
  ordered_normcounts <- tidyr::gather(ordered_normcounts[,1:101], gene, count, -gene_id)
  
  ylim <- max(ordered_counts$count, ordered_normcounts$count)
  
  counts_boxplot <- ggplot2::ggplot(ordered_counts, ggplot2::aes(x=gene, y=count)) + ggplot2::geom_boxplot() + ggplot2::scale_y_continuous(limits = c(0, ylim*1.1))
  normcounts_boxplot <- ggplot2::ggplot(ordered_normcounts, ggplot2::aes(x=gene, y=count)) + ggplot2::geom_boxplot() + ggplot2::scale_y_continuous(limits = c(0, ylim*1.1))
  
  counts_boxplot <- counts_boxplot + ggplot2::labs(x='Gene', y='Counts') + ggplot2::ggtitle('Gene expression (Sampled from 100 cells)', subtitle = "Before normalisation")
  normcounts_boxplot <- normcounts_boxplot + ggplot2::labs( x='Gene', y='Counts') + ggplot2::ggtitle('Gene expression (Sampled from 100 cells)', subtitle = sprintf("After normalisation via %s", progressLog(object)$NormalisationMethod))
  
  counts_boxplot <- counts_boxplot + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  normcounts_boxplot <- normcounts_boxplot + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  output_list[["sampled_cell_gene_expression"]] <- list(counts = counts_boxplot, normcounts = normcounts_boxplot)
  return(output_list)
}

#' plotControlHist
#' 
#' Generates a histogram of percentage of control per samples.
#' 
#' @param object An \linkS4class{EMSet} object.
#' @param control Name of control group to plot.
#' @return A ggplot2 histogram
#' 
#' @examples 
#' # Load unprocessed EMSet
#' em_set <- ascend::raw_set
#' 
#' # Plot controls
#' control_histogram <- plotControlHist(em_set, control = "Mt")
#' 
#' @export
#' @importFrom dplyr left_join
#' 
plotControlHist <- function(object, control = NULL){
  # Silence R Check
  percentage <- "Shhh"
  
  # User must specify a control group
  if (is.null(control)){
    stop("Please specify a control group to plot.")
  }
  
  # Control group must be specified in the EMSet
  if (!(control %in% rowInfo(object)$control_group)){
    stop("Please check that your control group has been defined in the object.")
  }
  
  # Get metrics
  metrics_df <- as.data.frame(SummarizedExperiment::colData(object))
  
  # Get cell info
  cell_info <- as.data.frame(colInfo(object))
  
  # Combine data
  metrics_df <- dplyr::left_join(cell_info, metrics_df, by = "cell_barcode")
  metrics_df <- metrics_df[, c("cell_barcode", "batch", sprintf("qc_%s_pct_counts", control))]
  colnames(metrics_df) <- c("cell_barcode", "batch", "percentage")
  pct_hist <- ggplot2::ggplot(metrics_df, ggplot2::aes(percentage))
  pct_hist <- pct_hist + ggplot2::geom_histogram(binwidth = 2, fill = "#000000")
  pct_hist <- pct_hist + ggplot2::ggtitle(sprintf("Proportion of control: %s", control))
  pct_hist <- pct_hist + ggplot2::xlab("% control") + ggplot2::ylab("Number of cells") + ggplot2::theme_bw()
  return(pct_hist)
}  

#' plotFeatureHist
#' 
#' Generates a histogram of all cell feature counts.
#' 
#' @param object An \linkS4class{EMSet} object.
#' @examples 
#' # Load unprocessed EMSet
#' em_set <- ascend::raw_set
#' 
#' # Use function to plot features
#' feature_plot <- plotFeatureHist(em_set)
#' 
#' @export
#' @importFrom dplyr left_join
#' @return A ggplot2 glob containing the histogram.
#' 
plotFeatureHist <- function(object){
  # For R Check
  qc_nfeaturecounts <- "Shhhh"
  
  # Get metrics
  metrics_df <- as.data.frame(SummarizedExperiment::colData(object))
  
  # Get cell info
  cell_info <- as.data.frame(colInfo(object))
  
  # Combine data
  metrics_df <- dplyr::left_join(cell_info, metrics_df, by = "cell_barcode")
  metrics_df <- metrics_df[, c("cell_barcode", "batch", "qc_nfeaturecounts")]
  feature_hist <- ggplot2::ggplot(metrics_df, ggplot2::aes(qc_nfeaturecounts))
  feature_hist <- feature_hist + ggplot2::geom_histogram(breaks = seq(min(metrics_df$qc_nfeaturecounts), 
                                                                      max(metrics_df$qc_nfeaturecounts), 
                                                                      by = 1000), fill = "#000000")
  feature_hist <- feature_hist + ggplot2::ggtitle("Distribution of feature counts for all cells")
  feature_hist <- feature_hist + ggplot2::xlab("Feature counts") + ggplot2::ylab("Number of cells") + ggplot2::theme_bw()
  return(feature_hist)
}

#' plotTopGenesPerSample
#'
#' Generates a violin-beeswarm plot of top gene counts and contribution to to 
#' a cell's total expression.
#' 
#' @param object An \code{\linkS4class{EMSet}} object.
#' @return A ggplot2 object containing a violin-beeswarm plot
#' @examples
#' # Load unprocessed EMSet
#' em_set <- ascend::raw_set
#' 
#' # Plot top genes
#' topgenes <- plotTopGenesPerSample(em_set)
#' 
#' @export
#' @importFrom dplyr left_join
#' @importFrom ggbeeswarm geom_quasirandom 
#'
plotTopGenesPerSample <- function(object){
  # Silence R Check
  batch <- "Shhh"
  top_500_pct <- "Shhh"
  top_100_pct <- "Shhh"
  
  # Check if control has been added
  col_info <- as.data.frame(colInfo(object))
  col_data <- as.data.frame(SummarizedExperiment::colData(object))
  row_info <- as.data.frame(rowInfo(object))
  row_data <- as.data.frame(SummarizedExperiment::rowData(object))
  
  cell_info <- dplyr::left_join(col_info, col_data, by = "cell_barcode")
  batch_names <- unique(cell_info$batch)
  
  # What information do we want?
  # Most expressed genes
  expression_matrix <- SingleCellExperiment::counts(object)
  ordered_row_data <- row_data[order(row_data$qc_topgeneranking), ]
  top_500_counts <- expression_matrix[ordered_row_data[1:500,1], ]
  top_100_counts <- expression_matrix[ordered_row_data[1:500,1], ]
  top_500_pct <- 100 * (Matrix::colSums(top_500_counts)/Matrix::colSums(expression_matrix))
  top_100_pct <- 100 * (Matrix::colSums(top_100_counts)/Matrix::colSums(expression_matrix))
  
  # Volcano plots for control and sample
  plot_dataframe <- cell_info[ , c("cell_barcode", "batch")]
  plot_dataframe$top_500_pct <- as.vector(top_500_pct)
  plot_dataframe$top_100_pct <- as.vector(top_100_pct)
  plot_dataframe$batch <- factor(plot_dataframe$batch, levels = batch_names)
  
  gradient_ramp <- ggplot2::scale_colour_gradient(name = "% Top 100 gene expression", low="#001b7f", high="#f1d351")
  plot_title <- "% Top gene expression to total cell expression"
  
  violin_plot <- ggplot2::ggplot(plot_dataframe, ggplot2::aes(x = batch, 
                                                              y = top_500_pct, 
                                                              colour = top_100_pct))
  violin_plot <- violin_plot + ggplot2::geom_violin(size = 1, scale = "width", colour = NA) + 
    ggbeeswarm::geom_quasirandom(shape = 16, size=5, alpha=0.5, dodge.width=0.5) + ggplot2::theme_bw()
  violin_plot <- violin_plot + gradient_ramp + 
    ggplot2::scale_x_discrete(name = "Sample", limits = unique(plot_dataframe$batch)) + 
    ggplot2::xlab("Sample") + ggplot2::ylab("% Top 500 gene expression") + 
    ggplot2::ggtitle(plot_title)
  return(violin_plot)
}


#' plotLibsizePerSample
#'
#' Generates a violin-beeswarm plot of library sizes per sample.
#' 
#' @param object An \code{\linkS4class{EMSet}} object.
#' @return A ggplot2 object containing a violin-beeswarm plot
#' @examples 
#' # Load unprocessed EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot library size per sample
#' libsize_per_sample <- plotLibsizePerSample(EMSet)
#' 
#' @export
#' @importFrom dplyr left_join
#' @importFrom ggbeeswarm geom_quasirandom 
#'
plotLibsizePerSample <- function(object){
  # Silence R check
  batch <- "Shhh"
  total_counts <- "Shhh"
  feature_counts <- "Shhh"
  
  # Check if control has been added
  col_info <- as.data.frame(colInfo(object))
  col_data <- as.data.frame(SummarizedExperiment::colData(object))
  row_info <- as.data.frame(rowInfo(object))
  row_data <- as.data.frame(SummarizedExperiment::rowData(object))
  
  cell_info <- dplyr::left_join(col_info, col_data, by = "cell_barcode")
  
  # What information do we want?
  # Volcano plots for control and sample
  plot_dataframe <- cell_info[ , c("cell_barcode", "batch", "qc_libsize", "qc_nfeaturecounts")]
  colnames(plot_dataframe) <- c("cell_barcode", "batch", "total_counts", "feature_counts")
  plot_dataframe$batch <- factor(plot_dataframe$batch, levels = unique(plot_dataframe$batch)) 
  
  gradient_ramp <- ggplot2::scale_colour_gradient(low="#001b7f", high="#f1d351")
  plot_title <- sprintf("Total library size and feature counts per sample")
  
  violin_plot <- ggplot2::ggplot(plot_dataframe, ggplot2::aes(x = factor(batch), y = total_counts, colour = feature_counts))
  violin_plot <- violin_plot + ggplot2::geom_violin(size = 1, scale = "width", colour = NA) + ggbeeswarm::geom_quasirandom(shape = 16, size=5, alpha=0.5, dodge.width=0.5)
  violin_plot <- violin_plot + gradient_ramp + ggplot2::scale_x_discrete(name = "Feature Counts", limits = unique(plot_dataframe$batch)) + ggplot2::xlab("Sample") + ggplot2::ylab("Total counts")  + ggplot2::ggtitle(plot_title)
  return(violin_plot)
}


#' plotControlPctPerSample
#'
#' Generates a violin/beeswarm plot of percentage of control expression per 
#' sample. Called by \code{\link[ascend:plotGeneralQC]{plotGeneralQC}}.
#' 
#' @param object An \code{\linkS4class{EMSet}}
#' @param control Name of the control you would like to plot.
#' @return A ggplot2 object containing a violin-beeswarm plot
#' 
#' @examples
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Generate control violin
#' control_violin <- plotControlPctPerSample(EMSet, control = "Mt")
#' 
#' @export
#' @importFrom dplyr left_join
#' @importFrom ggbeeswarm geom_quasirandom 
#'
plotControlPctPerSample<- function(object, control = NULL){
  # Silence R check
  batch <- "Shhh"
  percentage <- "Shhh"
  
  # Check if control has been added
  col_info <- as.data.frame(colInfo(object))
  col_data <- as.data.frame(SummarizedExperiment::colData(object))
  row_info <- as.data.frame(rowInfo(object))
  row_data <- as.data.frame(SummarizedExperiment::rowData(object))

  if (! control %in% row_info[ , "control_group"]){
    stop("Please check that you have defined this control in your dataset.")
  }
    
  cell_info <- dplyr::left_join(col_info, col_data, by = "cell_barcode")
  
  # What information do we want?
  # Volcano plots for control and sample
  plot_dataframe <- cell_info[ , c("cell_barcode", "batch", sprintf("qc_%s_pct_counts", control), sprintf("qc_%s_ncounts", control))]
  colnames(plot_dataframe) <- c("cell_barcode", "batch", "percentage", "counts")
  plot_dataframe$batch <- factor(plot_dataframe$batch, levels = unique(plot_dataframe$batch)) 
  
  gradient_ramp <- ggplot2::scale_colour_gradient(low="#001b7f", high="#f1d351")
  plot_title <- sprintf("Percentage of reads mapped to %s genes", control)
  
  violin_plot <- ggplot2::ggplot(plot_dataframe, ggplot2::aes(x = factor(batch), y = percentage, colour = counts))
  violin_plot <- violin_plot + ggplot2::geom_violin(size = 1, scale = "width", colour = NA) + 
    ggbeeswarm::geom_quasirandom(shape = 16, size=5, alpha=0.5, dodge.width=0.5)
  violin_plot <- violin_plot + ggplot2::theme_bw() + gradient_ramp + 
    ggplot2::scale_x_discrete(limits = unique(plot_dataframe$batch)) + 
    ggplot2::xlab("Sample") + ggplot2::ylab("% Reads")  + 
    ggplot2::ggtitle(plot_title)
  return(violin_plot)
}

#' plotLibsizeBoxplot
#' 
#' Generates a barplot of each cell's library size, and arranges them in 
#' descending order.
#' 
#' @param object An \linkS4class{EMSet}.
#' @examples
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot libsize boxplot
#' libsize_barplot <- plotLibsizeBoxplot(EMSet)
#' 
#' @export
#' @importFrom dplyr left_join
#' @importFrom stats reorder
#' @return A ggplot2 glob containing the barplot.
#' 
plotLibsizeBoxplot <- function(object){
  # Silence R Check
  cell_barcode <- "Shhh"
  qc_libsize <- "Shhh"
  batch <- "Shhh"
  
  # Get metrics
  metrics_df <- as.data.frame(SummarizedExperiment::colData(object))
  
  # Get cell info
  cell_info <- as.data.frame(colInfo(object))
  
  # Need to marry the information as we have separated the stats and 
  # metadata
  
  # Combine data
  metrics_df <- dplyr::left_join(cell_info, metrics_df, by = "cell_barcode")
  metrics_df <- metrics_df[, c("cell_barcode", "batch", "qc_libsize")]
  metrics_df$batch <- factor(metrics_df$batch, levels = sort(unique(metrics_df$batch)))

  libsize_boxplot <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = batch, y = qc_libsize)) + ggplot2::geom_boxplot()
  libsize_boxplot <- libsize_boxplot + ggplot2::ggtitle("Library sizes per batch") + ggplot2::xlab("Sample")
  libsize_boxplot <- libsize_boxplot + ggplot2::ylab("Library size") + ggplot2::theme_bw()
  return(libsize_boxplot)
}

#' plotLibsizeHist
#' 
#' Generates a histogram of all cell library sizes.
#' 
#' @param object An \linkS4class{EMSet}.
#' @examples
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot libsize histograms
#' libsize_hist <- plotLibsizeHist(EMSet)
#' 
#' @export
#' @importFrom dplyr left_join
#' @return A ggplot2 glob containing the histogram.
#' 
plotLibsizeHist <- function(object){
  # Silence checks
  qc_libsize <- "Shhh"
  
  # Get metrics
  metrics_df <- as.data.frame(colData(object))
  
  # Get cell info
  cell_info <- as.data.frame(object@colInfo)
  
  # Combine data
  metrics_df <- dplyr::left_join(cell_info, metrics_df, by = "cell_barcode")
  metrics_df <- metrics_df[, c("cell_barcode", "batch", "qc_libsize")]
  libsize_hist <- ggplot2::ggplot(metrics_df, ggplot2::aes(qc_libsize))
  libsize_hist <- libsize_hist + ggplot2::geom_histogram(breaks = seq(min(metrics_df$qc_libsize), 
                                                                      max(metrics_df$qc_libsize), 
                                                                      by = 1000), fill = "#000000")
  libsize_hist <- libsize_hist + ggplot2::ggtitle("Distribution of library sizes for all cells")
  libsize_hist <- libsize_hist + ggplot2::xlab("Library size") + ggplot2::ylab("Number of cells") + ggplot2::theme_bw()
  return(libsize_hist)
}

#' plotAverageGeneCount
#' 
#' Generates a histogram of average counts for each gene.
#' 
#' @param object An \linkS4class{EMSet}.
#' @param metric Scale to plot data by - average (DEFAULT), log2 or log10.
#' 
#' @examples 
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot average gene count
#' average_genes <- plotAverageGeneCount(EMSet, metric = "average")
#' 
#' # Plot log2 average gene count
#' average_gene_2 <- plotAverageGeneCount(EMSet, metric = "log2")
#' 
#' # Plot log10 average gene count
#' average_gene_10 <- plotAverageGeneCount(EMSet, metric = "log10")
#' 
#' @export
#' @importFrom dplyr left_join
#' @return A ggplot2 glob containing the histogram.
#' 
plotAverageGeneCount <- function(object, metric = c("average", "log2", "log10")){
  # Silence R Check
  count <- "Shhh"
  
  if (is.null(metric)){
    metric <- "average"
  }
  
  # Get metrics
  gene_id_name <- colnames(rowData(object))[1]
  metrics_df <- as.data.frame(SummarizedExperiment::rowData(object))
  metrics_df <- metrics_df[ , c(gene_id_name, "qc_meancounts")]
  
  if (metric == "log2"){
    label <- Log[2]~"Average Count"
    metrics_df$count <- log2(metrics_df$qc_meancounts)
    metrics_df <- metrics_df[which(!(is.infinite(metrics_df$count))), ]
    break_width <- max(abs(metrics_df$count))*0.1
  }
  if (metric == "log10"){
    label <- Log[10]~"Average Count"
    metrics_df$count <- log10(metrics_df$qc_meancounts)
    metrics_df <- metrics_df[which(!(is.infinite(metrics_df$count))), ]
    break_width <- max(abs(metrics_df$count))*0.1
  } else{
    label <- "Average Count"
    metrics_df$count <- metrics_df$qc_meancounts
    break_width <- max(abs(metrics_df$count))*0.1
  }
  
  ac_hist <- ggplot2::ggplot(metrics_df, ggplot2::aes(count))
  ac_hist <- ac_hist + ggplot2::geom_histogram(breaks = seq(min(metrics_df$count), 
                                                            max(metrics_df$count), 
                                                            by = break_width), fill = "#000000")
  ac_hist <- ac_hist + ggplot2::ggtitle(label)
  ac_hist <- ac_hist + ggplot2::xlab(label) + ggplot2::ylab("Number of genes") + ggplot2::theme_bw()
  return(ac_hist)
}

#' plotTopGenesBoxplot
#'
#' Generates a boxplot using \link[ggplot2]{geom_boxplot} of the most expressed 
#' genes in the dataset, in a range defined by the user.
#' 
#' @param object A \code{\linkS4class{EMSet}}.
#' @param n Number of genes to be plotted.
#' @param controls Include control genes in plot (Default: TRUE).
#' 
#' @return A ggplot glob containing a box scatter plot that represents the 
#' expression of the most highly expressed genes
#' @examples
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot top gene expression
#' top_genes <- plotTopGenesBoxplot(EMSet, n = 20, controls = FALSE)
#' 
#' @importFrom tidyr gather
#' @export
plotTopGenesBoxplot <- function(object, n = 20, controls = TRUE){
  # Silence R check
  gene <- "Shhh"
  pct_expression <- "Shhh"
  cell_barcode <- "Shhh"
  
  # Prep the data to feed into the function
  row_info <- rowInfo(object)
  controls <- FALSE
  if(!controls){
    # Check if controls are excluded
    if (progressLog(object)$controls){
      all_controls <- unique(row_info$control_group[which(!is.na(row_info$control_group))])
      object <- excludeControl(object, control = all_controls)      
    }
  }
  
  # Extract required data
  plot_title <- sprintf("Top %i expressed genes", n)
  expression_matrix <- SingleCellExperiment::counts(object)
  row_info <- rowInfo(object)
  row_data <- SummarizedExperiment::rowData(object)
  col_data <- colData(object)
  
  # Sort genes by rank
  row_data <- row_data[order(row_data$qc_topgeneranking), ]
  
  # Subset n amount of genes and get related information
  plot_data <- row_data[1:n, ]
  top_gene_list <- plot_data[, 1]
  top_gene_counts <- expression_matrix[top_gene_list, ]
  total_counts_per_cell <- col_data$qc_libsize
  pct_exprs_per_cell <- top_gene_counts/total_counts_per_cell
  
  # Arrange data so it is nice
  if (is(expression_matrix, "sparseMatrix")){
    pct_exprs_per_cell <- as.matrix(Matrix::t(pct_exprs_per_cell))
  }
  
  pct_exprs_per_cell <- as.data.frame(pct_exprs_per_cell)
  pct_exprs_per_cell$cell_barcode <- rownames(pct_exprs_per_cell)
  gathered_pct_exprs_per_cell <- tidyr::gather(pct_exprs_per_cell, key = gene, value = pct_expression, -cell_barcode)
  gathered_pct_exprs_per_cell$gene <- factor(gathered_pct_exprs_per_cell$gene, levels = rev(top_gene_list))
  
  z_theme <- function() {
    # From perceptions - https://github.com/zonination/perceptions.
    
    palette <- RColorBrewer::brewer.pal("Greys", n=9)
    color.background = palette[1]
    color.grid.major = palette[5]
    color.axis.text = palette[7]
    color.axis.title = palette[7]
    color.title = palette[8]
    # Begin construction of chart
    ggplot2::theme_bw(base_size=9) +
      # Set the entire chart region to a light gray color
      ggplot2::theme(panel.background = ggplot2::element_rect(fill=color.background, color=color.background)) +
      ggplot2::theme(plot.background= ggplot2::element_rect(fill=color.background, color=color.background)) +
      ggplot2::theme(panel.border= ggplot2::element_rect(color=color.background)) +
      # Format the grid
      ggplot2::theme(panel.grid.major= ggplot2::element_line(color=color.grid.major,size=.25)) +
      ggplot2::theme(panel.grid.minor= ggplot2::element_blank()) +
      ggplot2::theme(axis.ticks= ggplot2::element_blank()) +
      # Format the legend, but hide by default
      ggplot2::theme(legend.position="none") +
      ggplot2::theme(legend.background = ggplot2::element_rect(fill=color.background)) +
      ggplot2::theme(legend.text = ggplot2::element_text(size=7,color=color.axis.title)) +
      # Set title and axis labels, and format these and tick marks
      ggplot2::theme(plot.title = ggplot2::element_text(color=color.title, size=20, vjust=1.25)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=14, color=color.axis.text)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size=14, color=color.axis.text)) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size=16, color=color.axis.title, vjust=0)) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(size=16, color=color.axis.title, vjust=1.25))
  }
  
  
  boxplot <- ggplot2::ggplot(gathered_pct_exprs_per_cell, ggplot2::aes(x = gene, y = pct_expression)) + ggplot2::geom_boxplot(ggplot2::aes(fill=gene), alpha=.5, outlier.colour = NULL)
  boxplot <- boxplot + ggplot2::coord_flip() + z_theme() + ggplot2::xlab("Gene") + ggplot2::ylab("% Expression") + ggplot2::ggtitle(plot_title)
  return(boxplot)
}

#' plotSmoothScatter
#' 
#' Generates a smooth scatter of of average counts for each gene per cell.
#' 
#' @param object An \linkS4class{EMSet}.
#' @examples
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' smooth_scatter <- plotSmoothScatter(EMSet)
#' 
#' @return A ggplot2 glob containing the histogram.
#' @importFrom graphics smoothScatter
#' @export
#' 
plotSmoothScatter <- function(object){
  metrics_df <- SummarizedExperiment::rowData(object)
  smooth_plot <- smoothScatter(log10(metrics_df$qc_meancounts), 
                               metrics_df$qc_ncells, 
                               xlab=expression(Log[10]~"Average Count"),
                               ylab="Number of expressing cells", 
                               main = expression(Log[10]~" Average gene expression across cells"))
  return(smooth_plot)
}


#' plotGeneNumber
#' 
#' Generates a histogram of the number of genes per cell.
#' 
#' @param object An \linkS4class{EMSet}
#' @return A histogram containing the distribution of number of genes per cell.
#' 
#' @examples
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot gene numbers
#' gene_number_plot <- plotGeneNumber(EMSet)
#' 
#' @importFrom SummarizedExperiment colData
#' @export
plotGeneNumber <- function(object){
  # silence R Check
  qc_ngenes <- "Shhh"
  
  metrics_df <- as.data.frame(SummarizedExperiment::colData(object))
  plot <- ggplot2::ggplot(metrics_df, ggplot2::aes(qc_ngenes))
  plot <- plot + ggplot2::geom_histogram(binwidth = 100)
  plot <- plot + ggplot2::ggtitle("Number of genes per cell") + 
    ggplot2::xlab("Number of genes") + ggplot2::ylab("Number of cells") +
    ggplot2::theme_bw()
  return(plot)
}

#' plotGeneralQC
#'
#' This function generates a series of plots that can be used to assess the
#' present quality of an EMSet.
#' 
#' The plots are as follows:
#' \itemize{
#' \item{\strong{Library Size per Cell}: A series of barplots depicting the 
#' library size of each cell in descending order.}
#' \item{\strong{Library Size}: A histogram depicting the distribution of 
#' library sizes across the dataset.}
#' \item{\strong{Average Gene Count}: A histogram depicting number of cells vs 
#' mean gene expression.}
#' \item{\strong{Average Gene Count (Log2)}: A histogram depicting number of 
#' cells vs Log2 mean gene expression.}
#' \item{\strong{Average Gene Count (Log10)}: A histogram depicting number of 
#' cells vs Log2 mean gene expression.}
#' \item{\strong{Log 10 Average Gene Count (Smooth Scatter)}: Smooth scatter 
#' plot of number of cells vs Log10 mean gene expression.}
#' \item{\strong{Top Genes Per Sample}: Violin/beehive plots depicting 
#' proportion of top genes to total cell expression per sample.}
#' \item{\strong{Top Gene Expression}: Boxplots depicting top 25 genes in terms 
#' of expression.}
#' }

#' If controls are defined, two additional plots are generated:
#' \itemize{
#' \item{\strong{Percentage Control Expression}: Histograms depicting number of 
#' cells with the contribution of control genes as a percentage of total counts.}
#' \item{\strong{Proportion of Control Expression}: Violin/beehive plots for 
#' each sample and control, depicting the proportion of controls to total 
#' expression.}
#' }
#' 
#' @param object An \code{\linkS4class{EMSet}} object.
#' @examples 
#' # Load EMSet
#' EMSet <- ascend::raw_set
#' 
#' # Plot general QC
#' general_qc_plots <- plotGeneralQC(EMSet)
#' 
#' @importFrom gridExtra grid.arrange
#' @export
#' @return A list of plot objects.
#' 
plotGeneralQC <- function(object){
  # Collect plots here!
  output_list <- list()
  
  print("Plotting library size plots...")
  # 1. Libsize barplot
  output_list[["libsize_boxplot"]] <- plotLibsizeBoxplot(object)
  
  # 2. Libsize histogram
  output_list[["libsize_histogram"]] <- plotLibsizeHist(object)
  
  print("Plotting average count plots...")
  # 3. Average count histograms
  output_list[["averagecount_histogram"]] <- plotAverageGeneCount(object, metric = "average")
  output_list[["averagecount_log2"]] <- plotAverageGeneCount(object, metric = "log2")
  output_list[["averagecount_log10"]] <- plotAverageGeneCount(object, metric = "log10")
  output_list[["averagecount_smoothScatter"]] <- plotSmoothScatter(object)
  
  # 4. Top gene expression
  print("Plotting top gene expression...")
  output_list[["topgenes_boxplot"]] <- plotTopGenesBoxplot(object, n = 25, controls = TRUE)
  output_list[["topgenes_violin"]] <- plotTopGenesPerSample(object)

    # If controls present...
  if ("controls" %in% names(progressLog(object))){
    print("Controls detected. Plotting control-specific plots...")
    
    # 1. Plot feature count histogram
    output_list[["featurecount_hist"]] <- plotFeatureHist(object)
    output_list[["ngenes_hist"]] <- plotGeneNumber(object)
    
    # 2. Proportion of controls
    control_hists <- lapply(names(progressLog(object)$set_controls), function(x) plotControlHist(object, control = x))
    names(control_hists) <- names(progressLog(object)$set_controls)
    control_violins <- lapply(names(progressLog(object)$set_controls), function(x) plotControlPctPerSample(object, control = x))
    names(control_violins) <- names(progressLog(object)$set_controls)
    output_list[["control_hists"]] <- control_hists
    output_list[["control_violins"]] <- control_violins
  }
  
  print("General QC plots complete!")
  return(output_list)
}
