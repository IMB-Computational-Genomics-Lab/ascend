library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)

devtools::load_all("~/CodeRepositories/ascend-dev")
object <- readRDS("~/Documents/NK_d30_QFNPC_updated.RDS")

# Plot manually
# Extract required values from the object
hclust.obj <- object@Clusters$Hclust
#hclust.obj$labels <- rep("", length(hclust.obj$labels))
optimal.height <- object@Clusters$OptimalTreeHeight
nclusters <- object@Clusters$NumberOfClusters
cluster.list <- object@Clusters$Clusters


# Count table
cluster.df <- as.data.frame(table(cluster.list))
dendro.obj <- as.dendrogram(hclust.obj)

# Count table
cluster.df <- as.data.frame(table(cluster.list))
clusters <- as.vector(cluster.list[order.dendrogram(dendro.obj)])

# Reorder frequencies
cluster.df <- cluster.df[match(unique(clusters), cluster.df$cluster.list), ]

group_func <- function(x, clusters = clusters){
  sizes <- sapply(x, function(y) sum(clusters == y))
  return(sizes)
}

coloured.dendro <- dendextend::color_branches(dendro.obj, clusters = clusters)
#coloured.dendro <- dendextend::rotate(coloured.dendro, names(cluster.list[order.dendrogram(dendro.obj)])
coloured.dendro %>% dendextend::branches_color(coloured.dendro, clusters = clusters)

coloured.dendro %>% 
coloured.dendro <- dendextend::set_labels(coloured.dendro, rep("", length(hclust.obj$labels)))
#k.dendro <- dendextend::colour_branches(dendro.obj, k = nclusters, groupLabels = cluster.df$Freq)


# Apply labels directly to dendrogram
# Clusters need to be the same length as the labels, hence why we can't replace them yet
label.vec <- sapply(clusters, function(x) cluster.df[which(cluster.df$cluster.list == x), "Freq"])


#coloured.dendro <- dendextend::color_branches(dendro.obj, k = nclusters, groupLabels = unique(label.vec))

#ordered.cluster <- as.vector(cluster.list[order.dendrogram(coloured.dendro)])
#cluster.df <- cluster.df[match(unique(ordered.cluster), cluster.df$cluster.list), ]
#dendro.labels <- cluster.df$Freq

# Get the list of colours in the dendrogram and count them
colour.values <- as.data.frame(table(coloured.dendro %>% dendextend::get_leaves_branches_col()))

# These sometimes don't match the order in the object, so we have to get the order from the dendrogram
colour.order <- unique(coloured.dendro %>% dendextend::get_leaves_branches_attr("col"))
colour.values <- colour.values[match(colour.order, colour.values$Var1), ]
dendro.values <- colour.values$Freq


#coloured.dendro <- dendextend::color_branches(dendro.obj, k = nclusters, groupLabels = dendro.values)

# Match to cluster labels

# Add labels to the dendro now
ordered.cluster <- as.vector(cluster.list[hclust.obj$labels[order.dendrogram(coloured.dendro)]])
cluster.df <- cluster.df[match(unique(ordered.cluster), cluster.df$cluster.list), ]
dendro.labels <- cluster.df$Freq