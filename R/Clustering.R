# ____________________________________________________________________________
# 1 Cluster Statistics ####
#' Generic function to compute distances, centroids and statistics after statistics.
#'
#' @param X an MCXpress object with MCA object (after MCA step)
#' @return MCXpress object with MCA object and Clustering Object (after MCA step)
calculate_cluster_centroids <- function(X, dim) {
  ## ............................................................................
  ## A Initialisation of variables ####
  X$cluster$labels <- X$cluster$labels %>% select(Cluster, Sample)
  cells_coord <- X$MCA$cells_principal[, 1:dim]
  genes_coord <- X$MCA$genes_standard[, 1:dim]
  X$cluster$labels <- X$cluster$labels %>% select(Cluster, Sample)
  labels <- X$cluster$labels
  nClusters <- X$cluster$nClusters

  ## ............................................................................
  ## B Calculate Cluster Centroids and Distance Related #### . .
  ## . . . . . . . . . . .  . . . . . . . . . . . . . . . . . .
  ## . . . . . . ..  a Cluster Centroids Calculation ####
  cat("\n Calculating Centroids \n")
  Cell_Coord_Cluster <- inner_join(by = "Sample", labels, cells_coord %>%
    rownames_to_column(var = "Sample")) %>% select_("-Sample")
  coord_centroids <- Cell_Coord_Cluster %>%
    group_by_("Cluster") %>%
    summarise_all(funs(mean))

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  b Inter Cluster Distance ####
  cat("\n Calculating Distance between Cluster Centroids \n")
  ClusterCoord_Axis <- coord_centroids %>% select(contains("Axis"))
  cluster_distances <- rdist(ClusterCoord_Axis) %>%
    set_colnames(coord_centroids$Cluster) %>%
    set_rownames(coord_centroids$Cluster)

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  c Gene-Cluster Distance Calculation ####
  cat("\n Calculating Distance betwenn Cluster and Genes \n")
  GCD <- rdist(
    x1 = coord_centroids %>% select(contains("Axis")),
    x2 = genes_coord
  ) %>%
    t() %>%
    set_colnames(coord_centroids$Cluster) %>%
    set_rownames(genes_coord %>% rownames()) %>%
    cbind(genes_coord ^ 2 %>%
      rowSums() %>%
      sqrt() %>%
      as.data.frame() %>%
      set_colnames("Origin"))

  gene_cluster_distances <- GCD %>%
    as.data.frame() %>%
    rownames_to_column(var = "Genes") %>%
    as_tibble()

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  d Closest Cluster to a Gene ####
  closest_cluster <- tibble(GCD %>% rownames(), colnames(GCD)[GCD %>%
    apply(1, which.min)]) %>% set_colnames(c("Genes", "Cluster"))
  closest_cluster <- closest_cluster %>% separate(Genes, into = c(
    "Gene",
    "bin"
  ), sep = "-bin")

  ## ............................................................................
  ## C Visualisation Cluster ####

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  a Cell Space GGplot ####
  plot1 <- ggplot2::ggplot(cells_coord %>%
    rownames_to_column(var = "Sample") %>%
    inner_join(X$cluster$labels, by = "Sample"), ggplot2::aes(
    x = Axis1,
    y = Axis2
  )) + geom_point(ggplot2::aes(colour = Cluster)) +
    ggplot2::theme_bw() + ggplot2::guides(colour = ggplot2::guide_legend(title = "Cluster")) +
    ggplot2::ggtitle("Clustering Results in standard cell space") +
    ggplot2::theme(legend.text = ggplot2::element_text(
      colour = "black",
      size = 8
    ))

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  b Gene Space with Centroids GGplot ####
  plot2 <- ggplot2::ggplot(genes_coord, ggplot2::aes(
    x = Axis1,
    y = Axis2
  )) + ggplot2::geom_point(alpha = 0.5) + ggplot2::theme_bw() +
    ggplot2::ggtitle("Genes Low Dimensional Space") + geom_point(
      data = coord_centroids,
      aes(x = Axis1, y = Axis2, colour = Cluster), size = 3,
      alpha = 1
    ) + theme(legend.text = element_text(
      colour = "black",
      size = 8
    ))

  ## ............................................................................
  ## D Cluster Object Finalisation ####
  X$cluster <- list(
    labels, dim, nClusters, coord_centroids, cluster_distances,
    gene_cluster_distances, closest_cluster, plot1, plot2
  ) %>%
    set_names(c(
      "labels", "dim", "nClusters", "coord_centroids",
      "cluster_distances", "gene_cluster_distances", "closest_cluster",
      "plot1", "plot2"
    ))
  X$Shiny <- create_dashboard2(X)
  class(X$cluster) <- "Cluster"
  return(X)
}

# ____________________________________________________________________________
# 2 Clustering Methods ####

## ............................................................................
## A Supervised Clustering ####
#' Supervised Clustering on MCA
#'
#' @param X MCXpress object containing MCA object
#' @param Y A vector with Sample names as name and the coresponding cluster as value.
#' @param dim number of component to retain for distance calculation.
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{labels}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{gene_cluster_distances}{GSEA Parameter used}
#' \item{closest_cluster}{Indicate the closest cluster for a gene}
#' \item{coord_centroids}{GSEA Parameter used}
#' \item{plot1}{Clustering Visualisation in the MCA cell space}
#' \item{plot2}{Visualisation of the centroids in the MCA gene space}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' your_cluster <- gsub("\\d$|HFF_",replacement = "", x=colnames(MCX64553$ExpressionMatrix))
#' names(your_cluster) <- colnames(MCX64553$ExpressionMatrix)
#' MCX64553 <-  cluster_supervised(MCX64553, your_cluster,dim=2)
#' @export
cluster_supervised <- function(X, Y, dim) {
  X$cluster$labels <- Y
  names(X$cluster$labels) <- names(Y)
  X$cluster$labels <- tibble(names(X$cluster$labels), X$cluster$labels) %>%
    set_colnames(c("Sample", "Cluster"))
  X$cluster$nClusters <- X$cluster$labels$Cluster %>% n_distinct()
  X <- calculate_cluster_centroids(X, dim)
  return(X)
}

## ............................................................................
## B K-Means Clustering ####
#' K-Means Clustering on MCA
#'
#' @param X MCXpress object containing MCA object
#' @param k Number of cluster to compute
#' @param dim number of component to retain for distance calculation.
#' @param maxIter maximum number of iteration
#' @param nstart number of random sets of initial centers
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{labels}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{gene_cluster_distances}{GSEA Parameter used}
#' \item{closest_cluster}{Indicate the closest cluster for a gene}
#' \item{coord_centroids}{GSEA Parameter used}
#' \item{plot1}{Clustering Visualisation in the MCA cell space}
#' \item{plot2}{Visualisation of the centroids in the MCA gene space}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' MCX64553 <- cluster_kmeans(MCX64553, k=6, dim=2)
#' @export
cluster_kmeans <- function(X, k = 2, dim=2, maxIter = 10, nstart = 50) {
  cat("Performing Kmeans Clustering with ", k, " Cluster")

  Cluster <- X$MCA$cells_standard[, 1:dim] %>% kmeans(
    centers = k, iter.max = maxIter,
    nstart = nstart
  )
  X$cluster$labels <- paste0("Cluster", Cluster$cluster)
  names(X$cluster$labels) <- rownames(X$MCA$cells_standard)
  X$cluster$labels <- tibble(names(X$cluster$labels), X$cluster$labels) %>%
    set_colnames(c("Sample", "Cluster"))
  X$cluster$nClusters <- k
  X <- calculate_cluster_centroids(X, dim)
  return(X)
}
## ............................................................................
## C Hierarchical Clustering ####
#' Hierarchical Clustering on MCA
#'
#' @param X MCXpress object containing MCA object
#' @param method hclust method ('average', 'ward', etc..)
#' @param k integer indicating the number of cluster to obtain
#' @param dim number of component to retain for distance calculation.
#' @param h numeric indicating the heights where the tree should be cut
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{labels}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{gene_cluster_distances}{GSEA Parameter used}
#' \item{closest_cluster}{Indicate the closest cluster for a gene}
#' \item{coord_centroids}{GSEA Parameter used}
#' \item{plot1}{Clustering Visualisation in the MCA cell space}
#' \item{plot2}{Visualisation of the centroids in the MCA gene space}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' MCX64553 <- cluster_hclust(MCX64553, k=6, method="ward.D", dim=2)
#' @export
cluster_hclust <- function(X, dim=2, method = "average", k = NULL, h = NULL) {
  Distance <- X$MCA$cells_principal %>% dist()
  Cluster <- Distance %>% hclust(method = method)
  Cluster <- Cluster %>% cutree(k = k, h = h)
  X$cluster$labels <- tibble(paste0("Cluster", Cluster), (Cluster %>%
    names())) %>% set_names(c("Cluster", "Sample"))
  X$cluster$nClusters <- Cluster %>% n_distinct()
  X <- calculate_cluster_centroids(X, dim)
  return(X)
}


## ............................................................................
## D K_Medoids Clustering ####
#' K-Medoids Clusetering on the MCA Analysis
#'
#' @param X MCXpress object containing MCA object
#' @param k integer indicating the number of cluster to obtain
#' @param dim number of component to retain for distance calculation.
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{labels}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{gene_cluster_distances}{GSEA Parameter used}
#' \item{closest_cluster}{Indicate the closest cluster for a gene}
#' \item{coord_centroids}{GSEA Parameter used}
#' \item{plot1}{Clustering Visualisation in the MCA cell space}
#' \item{plot2}{Visualisation of the centroids in the MCA gene space}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' MCX64553 <- cluster_k_medoids(MCX64553, k=6,dim=2)
#' @export
cluster_k_medoids <- function(X, k = 2, dim=2) {
  Distance <-
    Cluster <- X$MCA$cells_standard[, 1:dim] %>%
    pam(k = k, cluster.only = TRUE)
  X$cluster$labels <- tibble(paste0("Cluster", Cluster), (Cluster %>%
    names())) %>% set_names(c("Cluster", "Sample"))
  X$cluster$nClusters <- Cluster %>% n_distinct()
  X <- calculate_cluster_centroids(X, dim)
  return(X)
}

#' Model Clustering on MCA Coordinates
#'
#' @param X MCXpress object containing MCA object
#' @param k integer indicating the number of cluster to obtain
#' @param dim number of component to retain for distance calculation.
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{labels}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{gene_cluster_distances}{GSEA Parameter used}
#' \item{closest_cluster}{Indicate the closest cluster for a gene}
#' \item{coord_centroids}{GSEA Parameter used}
#' \item{plot1}{Clustering Visualisation in the MCA cell space}
#' \item{plot2}{Visualisation of the centroids in the MCA gene space}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' MCX64553 <- cluster_mclust(MCX64553, k=6, dim=2)
#' @export
cluster_mclust <- function(X, k, dim=2) {
  Cluster <- X$MCA$cells_standard[, 1:dim] %>%
    mclust::Mclust(G = k) %>%
    use_series(classification) %>%
    set_names(nm = X$MCA$cells_standard %>% rownames())
  X$cluster$labels <- tibble(paste0("Cluster", Cluster), (Cluster %>%
    names())) %>% set_names(c("Cluster", "Sample"))
  X$cluster$nClusters <- Cluster %>% n_distinct()
  X <- calculate_cluster_centroids(X, dim)
  return(X)
}



## ............................................................................
## B K-Means Clustering ####
#' K-Means Clustering on MCA
#'
#' @param X MCXpress object containing MCA object
#' @param dim number of component to retain for distance calculation.
#' @param max_distance Distance threshold for edges conservation
#' @param ratio Percentages of the closest cells to be linked to
#' @param k number of nearest neighbour
#' @param mutual nearest neighbour mutual or not
#' @param methods clustering methods for graph
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{labels}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{gene_cluster_distances}{GSEA Parameter used}
#' \item{closest_cluster}{Indicate the closest cluster for a gene}
#' \item{coord_centroids}{GSEA Parameter used}
#' \item{plot1}{Clustering Visualisation in the MCA cell space}
#' \item{plot2}{Visualisation of the centroids in the MCA gene space}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' MCX64553 <- cluster_igraph(MCX64553,dim=3, k=3, mutual=TRUE,methods="fast_greedy")
#' @export
cluster_igraph <- function(X, dim=2, max_distance= Inf, ratio=1, k=NULL, mutual=TRUE, methods="louvain") {
  cat("Performing Graph Clustering")
  Distance <- X$MCA$cells_standard[, 1:dim] %>% dist() %>% as.matrix()
  if (k %>% is.null() %>% not()) {
    Distance <- Distance %>% knnalgorithm(k = k, mutual = mutual)
  }
  g <- igraph::graph.adjacency(Distance, weighted = TRUE, mode = "undirected")
  g <- igraph::delete.edges(g, which(igraph::E(g)$weight > max_distance))
  g <- igraph::delete.edges(g, which(igraph::E(g)$weight > sort(igraph::E(g)$weight)[ratio * length(igraph::E(g)$weight)]))
  g <- igraph::delete.edges(g, which(igraph::E(g)$weight == 0))
  if (methods == "louvain") {
    clusgraph <- g %>% igraph::cluster_fast_greedy()
  }
  if (methods == "fast_greedy") {
    clusgraph <- g %>% igraph::cluster_fast_greedy()
  }
  if (methods == "label_prop") {
    clusgraph <- g %>% igraph::cluster_label_prop()
  }
  if (methods == "walktrap") {
    clusgraph <- g %>% igraph::cluster_walktrap()
  }
  X$cluster$labels <- clusgraph$membership %>% paste0("Cluster", .) %>% set_names(X$ExpressionMatrix %>% colnames())
  names(X$cluster$labels) <- rownames(X$MCA$cells_standard)
  X$cluster$labels <- tibble::tibble(names(X$cluster$labels), X$cluster$labels) %>%
    set_colnames(c("Sample", "Cluster"))
  X$cluster$nClusters <- X$cluster$labels$Cluster %>% n_distinct()
  X <- calculate_cluster_centroids(X, dim)
  return(X)
}


knnalgorithm <- function(distance, mutual=TRUE, k=3) {
  n <- nrow(distance)
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    d <- sort(distance[i, ])
    A[i, distance[i, ] <= d[k + 1]] <- 1
  }
  diag(A) <- 0
  if (mutual) {
    for (i in 1:n) {
      A[i, ] <- A[i, ] & A[, i]
      A[, i] <- A[i, ]
    }
  }
  A <- A * distance
  return(A)
}
