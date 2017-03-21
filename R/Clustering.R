## Clustering Methods
removeSmallClusters = function(g, min = 1) {
  nClusters <- clusters(g)$no
  cSizes <- clusters(g)$csize
  smallClustersId <- which(cSizes <= min)
  cellsToRemoveId <- which(clusters(g)$membership %in% smallClustersId)
  newG <- delete_vertices(g, cellsToRemoveId)
  return(newG)
}

Calculate_Cluster_Centroids <- function(X) {

  Cells_Coord <- X$Dim_Red$Cells_Principal
  Cluster_Quali <- X$cluster$Cluster_Quali
  nClusters <- X$cluster$nClusters
  Cell_Coord_Cluster <- left_join(by = "Sample", Cluster_Quali, Cells_Coord %>%
    rownames_to_column(var = "Sample")) %>% select_("-Sample")

  Coord_Centroids <- Cell_Coord_Cluster %>% group_by_("Cluster") %>% summarise_each(funs(mean))


  ClusterCoord_Axis <- Coord_Centroids %>% select(contains("Axis"))
  Cluster_Distance <- rdist(ClusterCoord_Axis)
  colnames(Cluster_Distance) <- Coord_Centroids$Cluster
  rownames(Cluster_Distance) <- Coord_Centroids$Cluster



  Y <- rdist(x1 = Coord_Centroids %>% select(contains("Axis")), x2 = X$Dim_Red$Genes_Standard) %>%
    t()
  colnames(Y) <- Coord_Centroids$Cluster
  rownames(Y) <- X$Dim_Red$Genes_Standard %>% rownames
  Gene_Cluster_Distance <- Y %>% as.data.frame() %>% rownames_to_column(var = "Genes") %>%
    as_tibble()


  Closest_Cluster <- tibble(Y %>% rownames, colnames(Y)[Y %>% apply(1, which.min)]) %>%
    set_colnames(c("Genes", "Cluster"))
  Closest_Cluster <- Closest_Cluster %>% separate(Genes, into = c("Gene", "bin"),
    sep = "-bin")

  # ____________________________________________________________________________
  # Cluster Visualisation ####

  ## ............................................................................
  ## cluster ggplot1 ####
  Graph1 <- ggplot(X$Dim_Red$Cells_Standard %>% rownames_to_column(var = "Sample") %>%
    inner_join(X$cluster$Cluster_Quali, by = "Sample"), aes(x = Axis1, y = Axis2)) +
    geom_point(aes(colour = Cluster)) + theme_bw() + guides(colour = guide_legend(title = "Cluster")) +
    ggtitle("Clustering Results in standard cell space") + theme(legend.text = element_text(colour = "black",
    size = 8))

  ## ............................................................................
  ## cluster ggplot2 ####
  Graph2 <- ggplot(X$Dim_Red$Genes_Standard, aes(x = Axis1, y = Axis2)) + geom_point(alpha = 0.5) +
    theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(data = X$cluster$Coord_Centroids,
    aes(x = Axis1, y = Axis2, colour = Cluster), size = 3, alpha = 1) + theme(legend.text = element_text(colour = "black",
    size = 8))

#   ____________________________________________________________________________
#   final                                                                   ####

### Creation of cluster object                                              ####
  X$cluster <- list(Cluster_Quali, nClusters, Coord_Centroids, Cluster_Distance,
    Gene_Cluster_Distance, Closest_Cluster, Graph1, Graph2) %>% set_names(c("Cluster_Quali",
    "nClusters", "Coord_Centroids", "Cluster_Distance", "Gene_Cluster_Distance",
    "Closest_Cluster", "Graph1", "Graph2"))
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Shiny Cluster                                                           ####
  X$cluster$Shiny <- Create_Shiny_Cluster(X)
  class(X$cluster) <- "Cluster_Object"
  return(X)



}


Cluster_knn <- function(X, mutual = TRUE, k = 3, rm = 0) {
  n <- nrow(X$Dim_Red$Cell2Cell_Distance)
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    d <- sort(X$Dim_Red$Cell2Cell_Distance[i, ])
    A[i, X$Dim_Red$Cell2Cell_Distance[i, ] <= d[k + 1]] <- 1
  }
  diag(A) <- 0
  if (mutual) {
    for (i in 1:n) {
      A[i, ] <- A[i, ] & A[, i]
      A[, i] <- A[i, ]
    }
  }
  A <- A * X$Dim_Red$Cell2Cell_Distance
  g <- graph.adjacency(A, weighted = TRUE, mode = "undirected")
  cluster_info <- g %>% clusters
  X$cluster$Cluster_Quali <- paste0("Cluster", cluster_info$membership) %>% tibble(cluster_info$membership %>%
    names) %>% set_colnames(c("Cluster", "Sample"))
  X$cluster$nClusters <- cluster_info$no
  g <- Calculate_Cluster_Centroids(X)
  return(g)
}

Cluster_max_distance <- function(X, mutual = TRUE, max.distance = 0.5, rm = 0) {
  X$Dim_Red$Cell2Cell_Distance[row(X$Dim_Red$Cell2Cell_Distance) == col(X$Dim_Red$Cell2Cell_Distance)] <- 0
  g <- graph.adjacency(X$Dim_Red$Cell2Cell_Distance, weighted = TRUE, mode = "undirected")
  g <- delete_edges(g, which(E(g)$weight > max.distance))
  g <- delete_edges(g, which(E(g)$weight == 0))
  cluster_info <- g %>% clusters
  X$cluster$Cluster_Quali <- paste0("Cluster", cluster_info$membership) %>% tibble(cluster_info$membership %>%
    names) %>% set_colnames(c("Cluster", "Sample"))
  X$cluster$nClusters <- cluster_info$no
  g <- Calculate_Cluster_Centroids(X)
  return(g)
}


Cluster_percentage <- function(X, mutual = TRUE, shortest.rank.percent = 10, rm = 0) {
  X$Dim_Red$Cell2Cell_Distance[row(X$Dim_Red$Cell2Cell_Distance) == col(X$Dim_Red$Cell2Cell_Distance)] <- 0
  g <- graph.adjacency(X$Dim_Red$Cell2Cell_Distance, weighted = TRUE, mode = "undirected")
  g <- delete_edges(g, which(E(g)$weight > sort(E(g)$weight)[shortest.rank.percent/100 *
    length(E(g)$weight)]))
  g <- delete_edges(g, which(E(g)$weight == 0))
  cluster_info <- g %>% clusters
  X$cluster$Cluster_Quali <- paste0("Cluster", cluster_info$membership) %>% tibble(cluster_info$membership %>%
    names) %>% set_colnames(c("Cluster", "Sample"))
  X$cluster$nClusters <- cluster_info$no
  g <- Calculate_Cluster_Centroids(X)
  return(g)
}

Cluster_Personalised <- function(X, Y) {
  X$cluster$Cluster_Quali <- Y
  names(X$cluster$Cluster_Quali) <- names(Y)
  X$cluster$Cluster_Quali <- tibble(names(X$cluster$Cluster_Quali), X$cluster$Cluster_Quali) %>%
    set_colnames(c("Sample", "Cluster"))
  X$cluster$nClusters <- X$cluster$Cluster_Quali$Cluster %>% unique %>% length
  X <- Calculate_Cluster_Centroids(X)
  return(X)
}

Cluster_Kmeans <- function(X, nCluster, maxIter = 10, nstart = 50) {
  Distance <- X$Dim_Red$Cell2Cell_Distance
  Cells_Coord <- X$Dim_Red$Cells_Principal
  Cluster <- Distance %>% kmeans(centers = nCluster, iter.max = maxIter, nstart = nstart)
  X$cluster$Cluster_Quali <- paste0("Cluster", Cluster$cluster)
  names(X$cluster$Cluster_Quali) <- rownames(Distance)
  X$cluster$Cluster_Quali <- tibble(names(X$cluster$Cluster_Quali), X$cluster$Cluster_Quali) %>%
    set_colnames(c("Sample", "Cluster"))
  X$cluster$nClusters <- nCluster
  X <- Calculate_Cluster_Centroids(X)
  return(X)
}

