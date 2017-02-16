##Clustering Methods
removeSmallClusters = function(g, min = 1) {
  nClusters <- clusters(g)$no
  cSizes <- clusters(g)$csize
  smallClustersId <- which(cSizes <= min)
  cellsToRemoveId <-
    which(clusters(g)$membership %in% smallClustersId)
  newG <- delete_vertices(g, cellsToRemoveId)
  return(newG)
}


Calculate_Cluster_Centroids <- function(cluster, X) {
  cluster_info  <-  cluster %>%  clusters
  Cells_Coord   <-  X$Dim_Red$Cells_Principal
  X$cluster$Cluster_Quali   <-  paste0('Cluster', cluster_info$membership)
  names(X$cluster$Cluster_Quali)  <- cluster_info$membership %>%  names
  X$cluster$Cluster_Quali <-
    tibble(names(X$cluster$Cluster_Quali), X$cluster$Cluster_Quali) %>%  set_colnames(c("Sample", "Cluster"))
  Cell_Coord_Cluster <-
    left_join(by = 'Sample',
              X$cluster$Cluster_Quali,
              Cells_Coord %>%  rownames_to_column(var = 'Sample')) %>%  select_("-Sample")
  Coord_Centroids <-
    Cell_Coord_Cluster %>%  group_by_("Cluster") %>%  summarise_each(funs(mean))
  X$cluster$nClusters <- cluster_info$no
  X$cluster$Coord_Centroids <- Coord_Centroids
  X$cluster$Distance <-
    Calculate_Cluster_Cluster_Distance(X$cluster$Coord_Centroids)
  X <- Calculate_Distance_Centroids_to_Gene(X)
  X$cluster$Graph1 <-
    ggplot(X$Dim_Red$Cells_Standard, aes_string(x ="Axis1", y ="Axis2")) + geom_point(aes(colour =
                                                                               X$cluster$Cluster_Quali %>% as.vector())) + theme_bw() + guides(colour =
                                                                                                                                                 guide_legend(title = "Cluster")) + ggtitle("Cells Low Dimensional Space") +
    theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Graph2 <-
    ggplot(X$Dim_Red$Genes_Standard, aes_string(x ="Axis1", y ="Axis2")) + geom_point(alpha =
                                                                           0.5) + theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(
                                                                             data = X$cluster$Coord_Centroids,
                                                                             aes_string(x = "Axis1", y = "Axis2", colour = "Cluster"),
                                                                             size = 3,
                                                                             alpha = 1
                                                                           ) + theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Shiny <- Create_Shiny_Cluster(X)
  class(X$cluster) <- "Cluster_Object"
  return(X)
}


Cluster_knn <- function(X,
                        mutual = TRUE,
                        k = 3,
                        rm = 0) {
  n <- nrow(X$Dim_Red$Cell2Cell_Distance)
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    d <- sort(X$Dim_Red$Cell2Cell_Distance[i,])
    A[i, X$Dim_Red$Cell2Cell_Distance[i,] <= d[k + 1]] <- 1
  }
  diag(A) <- 0
  if (mutual) {
    for (i in 1:n) {
      A[i,] <- A[i,] & A[, i]
      A[, i] <- A[i,]
    }
  }
  A <- A * X$Dim_Red$Cell2Cell_Distance
  g <- graph.adjacency(A, weighted = TRUE, mode = "undirected")
  g <- Calculate_Cluster_Centroids(g, X)
  return(g)
}

Cluster_max_distance <-
  function(X,
           mutual = TRUE,
           max.distance = 0.5,
           rm = 0) {
    X$Dim_Red$Cell2Cell_Distance[row(X$Dim_Red$Cell2Cell_Distance) == col(X$Dim_Red$Cell2Cell_Distance)] <-
      0
    g <-
      graph.adjacency(X$Dim_Red$Cell2Cell_Distance ,
                      weighted = TRUE,
                      mode = "undirected")
    g <- delete_edges(g, which(E(g)$weight > max.distance))
    g <- delete_edges(g, which(E(g)$weight == 0))
    g <- Calculate_Cluster_Centroids(g, X)
    return(g)
  }


Cluster_percentage <-
  function(X,
           mutual = TRUE,
           shortest.rank.percent = 10,
           rm = 0) {
    X$Dim_Red$Cell2Cell_Distance[row(X$Dim_Red$Cell2Cell_Distance) == col(X$Dim_Red$Cell2Cell_Distance)] <-
      0
    g <-
      graph.adjacency(X$Dim_Red$Cell2Cell_Distance ,
                      weighted = TRUE,
                      mode = "undirected")
    g <-
      delete_edges(g, which(E(g)$weight > sort(E(g)$weight)[shortest.rank.percent /
                                                              100 * length(E(g)$weight)]))
    g <- delete_edges(g, which(E(g)$weight == 0))
    g <- Calculate_Cluster_Centroids(g, X)
    return(g)
  }

Cluster_Personalised <- function(X, Y) {
  Distance <- X$Dim_Red$Cell2Cell_Distance
  Cells_Coord <- X$Dim_Red$Cells_Principal
  X$cluster$Cluster_Quali <- paste0('Cluster_', Y)
  names(X$cluster$Cluster_Quali) <- rownames(Distance)
  X$cluster$Cluster_Quali <-
    tibble(names(X$cluster$Cluster_Quali), X$cluster$Cluster_Quali) %>%  set_colnames(c("Sample", "Cluster"))
  Cell_Coord_Cluster <-
    left_join(by = 'Sample',
              X$cluster$Cluster_Quali,
              Cells_Coord %>%  rownames_to_column(var = 'Sample')) %>%  select(-Sample)
  Coord_Centroids <-
    Cell_Coord_Cluster %>%  group_by(Cluster) %>%  summarise_each(funs(mean))
  X$cluster$nClusters <-
    X$cluster$Cluster_Quali$Cluster %>% unique %>% length
  X$cluster$Coord_Centroids <- Coord_Centroids
  X$cluster$Distance <-
    Calculate_Cluster_Cluster_Distance(X$cluster$Coord_Centroids)
  X <- Calculate_Distance_Centroids_to_Gene(X)
  X$cluster$Graph1 <-
    ggplot(X$Dim_Red$Cells_Standard, aes(x = Axis1, y = Axis2)) + geom_point(aes(colour =
                                                                               X$cluster$Cluster_Quali %>% as.vector())) + theme_bw() + guides(colour =
                                                                                                                                                 guide_legend(title = "Cluster")) + ggtitle("Cells Low Dimensional Space") +
    theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Graph2 <-
    ggplot(X$Dim_Red$Genes_Standard, aes(x = Axis1, y = Axis2)) + geom_point(alpha =
                                                                           0.5) + theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(
                                                                             data = X$cluster$Coord_Centroids,
                                                                             aes(x = Axis1, y = Axis2, colour = Cluster),
                                                                             size = 3,
                                                                             alpha = 1
                                                                           ) + theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Distance_Var <-
    X$cluster$Gene_Cluster_Distance %>% gather(-Genes, key = "Cluster", value = "Distance") %>% group_by(Genes) %>%  summarise(Dis_Var =
                                                                                                                                 var(Distance)) %>% mutate(Rank_Var = -Dis_Var %>% rank) %>% arrange(Rank_Var)
  X$cluster$Shiny <- Create_Shiny_Cluster(X)
  return(X)
}

Cluster_Kmeans <- function(X,
                           nCluster,
                           maxIter = 10,
                           nstart = 50)
{
  Distance <- X$Dim_Red$Cell2Cell_Distance
  Cells_Coord <- X$Dim_Red$Cells_Principal
  Cluster <-
    Distance %>% kmeans(centers = nCluster ,
                        iter.max = maxIter,
                        nstart = nstart)
  X$cluster$Cluster_Quali <- paste0('Cluster', Cluster$cluster)
  names(X$cluster$Cluster_Quali) <- rownames(Distance)
  X$cluster$Cluster_Quali <-
    tibble(names(X$cluster$Cluster_Quali), X$cluster$Cluster_Quali) %>%  set_colnames(c("Sample", "Cluster"))
  Cell_Coord_Cluster <-
    left_join(by = 'Sample',
              X$cluster$Cluster_Quali,
              Cells_Coord %>%  rownames_to_column(var = 'Sample')) %>%  select(-Sample)
  Coord_Centroids <-
    Cell_Coord_Cluster %>%  group_by(Cluster) %>%  summarise_each(funs(mean))
  X$cluster$nClusters <- nCluster
  X$cluster$Coord_Centroids <- Coord_Centroids
  X$cluster$TEST<- Cluster$cluster
  X$cluster$Distance <-
    Calculate_Cluster_Cluster_Distance(X$cluster$Coord_Centroids)
  X <- Calculate_Distance_Centroids_to_Gene(X)
  X$cluster$Calinski_Harabasz<-cluster.stats(d = (Distance*1000), clustering = Cluster$cluster)$ch
  X$cluster$Graph1 <-
    ggplot(X$Dim_Red$Cells_Standard, aes(x = Axis1, y = Axis2)) + geom_point(aes(colour =
                                                                               X$cluster$Cluster_Quali %>% as.vector())) + theme_bw() + guides(colour =
                                                                                                                                                 guide_legend(title = "Cluster")) + ggtitle("Cells Low Dimensional Space") +
    theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Graph2 <-
    ggplot(X$Dim_Red$Genes_Standard, aes(x = Axis1, y = Axis2)) + geom_point(alpha =
                                                                           0.5) + theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(
                                                                             data = X$cluster$Coord_Centroids,
                                                                             aes(x = Axis1, y = Axis2, colour = Cluster),
                                                                             size = 3,
                                                                             alpha = 1
                                                                           ) + theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Shiny <- Create_Shiny_Cluster(X)
  X$cluster$Distance_Var <-
    X$cluster$Gene_Cluster_Distance %>% gather(-Genes, key = "Cluster", value = "Distance") %>% group_by(Genes) %>%  summarise(Dis_Var =
                                                                                                                                 var(Distance)) %>% mutate(Rank_Var = -Dis_Var %>% rank) %>% arrange(Rank_Var)
  class(X$cluster) <- "Cluster_Object"
  return(X)
}

#Gene To Centroids distance calculation
Find_Closest_Cluster <- function(X) {
  Closest_Cluster <- colnames(X)[X %>% apply(1, which.min)]
  Gene <- X %>% rownames
  DF <- tibble(Gene, Closest_Cluster)
  DF <-
    DF %>%  separate(Gene, into = c('Gene', 'bin'), sep = '-bin') %>%  filter(bin ==
                                                                                1) %>%  select(-bin)
  return(DF)
}


Calculate_Distance_Centroids_to_Gene <- function(X) {
  A <- X$cluster$Coord_Centroids %>% select(contains('Axis'))
  B <- X$Dim_Red$Genes_Standard
  Y <- rdist(x1 = A, x2 = B) %>% t()
  colnames(Y) <- X$cluster$Coord_Centroids$Cluster
  rownames(Y) <- B %>% rownames
  X$cluster$Gene_Cluster_Distance <-
    Y  %>% as.data.frame() %>%  rownames_to_column(var = 'Genes') %>%  as_tibble()
  X$cluster$Closest_Cluster <- Find_Closest_Cluster(Y)
  return(X)
}

Calculate_Cluster_Cluster_Distance <- function(ClusterCoord) {
  ClusterCoord_Axis <- ClusterCoord %>%  select(contains('Axis'))
  Distance <- rdist(ClusterCoord_Axis)
  colnames(Distance) <- ClusterCoord$Cluster
  rownames(Distance) <- ClusterCoord$Cluster
  return(Distance)
}




