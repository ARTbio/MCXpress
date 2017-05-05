#   ____________________________________________________________________________
#   1 Cluster Statistics                                                    ####
#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
Calculate_Cluster_Centroids <- function(X) {
##  ............................................................................
##  A Initialisation of variables                                           ####
    cells_coord <- X$Dim_Red$Cells_Principal
    genes_coord <- X$Dim_Red$Genes_Standard
    Cluster_Quali<-X$cluster$Cluster_Quali
    nClusters<- X$cluster$nClusters

##  ............................................................................
##  B Calculate Cluster Centroids and Distance Related                      ####
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Cluster Centroids Calculation                                         ####
    cat("\n Calculating Centroids \n")
    Cell_Coord_Cluster <- inner_join(by = "Sample",
        Cluster_Quali, cells_coord %>% rownames_to_column(var = "Sample")) %>%
        select_("-Sample")
    Coord_Centroids <- Cell_Coord_Cluster %>% group_by_("Cluster") %>% summarise_each(funs(mean))

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### b Inter Cluster Distance                                                ####
    cat("\n Calculating Distance between Cluster Centroids \n")
    ClusterCoord_Axis <- Coord_Centroids %>% select(contains("Axis"))
    Cluster_Distance <- rdist(ClusterCoord_Axis) %>%  set_colnames(Coord_Centroids$Cluster) %>% set_rownames(Coord_Centroids$Cluster)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### c Gene-Cluster Distance Calculation                                     ####
    cat("\n Calculating Distance betwenn Cluster and Genes \n")
    GCD <- rdist(x1 = Coord_Centroids %>%
                                     select(contains("Axis")),
                                   x2 = genes_coord
                                   ) %>% t() %>%
      set_colnames(Coord_Centroids$Cluster) %>%
      set_rownames(genes_coord %>% rownames) %>%  cbind(genes_coord^2 %>% rowSums() %>% sqrt() %>% as.data.frame() %>% set_colnames("Origin"))

    Gene_Cluster_Distance<-GCD %>%
      as.data.frame() %>%
      rownames_to_column(var = "Genes") %>%
      as_tibble()

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### d Closest Cluster to a Gene                                             ####
    closest_Cluster <- tibble(GCD %>% rownames, colnames(GCD)[GCD %>% apply(1, which.min)]) %>% set_colnames(c("Genes", "Cluster"))
    closest_Cluster <- closest_Cluster %>% separate(Genes, into = c("Gene", "bin"),sep = "-bin")

##  ............................................................................
##  C Visualisation Cluster                                                 ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Cell Space GGplot                                                     ####
    graph1 <- ggplot(cells_coord %>%  rownames_to_column(var="Sample") %>%  inner_join(X$cluster$Cluster_Quali, by="Sample"),
        aes(x = Axis1, y = Axis2)) + geom_point(aes(colour = Cluster)) + theme_bw() + guides(colour = guide_legend(title = "Cluster")) +
        ggtitle("Clustering Results in standard cell space") + theme(legend.text = element_text(colour = "black",
        size = 8))

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### b Gene Space with Centroids GGplot                                      ####
   graph2 <- ggplot(genes_coord,
        aes(x = Axis1, y = Axis2)) + geom_point(alpha = 0.5) +
        theme_bw() + ggtitle("Genes Low Dimensional Space") +
        geom_point(data = Coord_Centroids,
            aes(x = Axis1, y = Axis2, colour = Cluster),
            size = 3, alpha = 1) + theme(legend.text = element_text(colour = "black",
        size = 8))

##  ............................................................................
##  D Cluster Object Finalisation                                           ####
    X$cluster<- list(Cluster_Quali, nClusters, Coord_Centroids, Cluster_Distance, Gene_Cluster_Distance, closest_Cluster, graph1, graph2) %>%
    set_names(c("Cluster_Quali","nClusters","Coord_Centroids", "Cluster_Distance", "Gene_Cluster_Distance", "closest_Cluster", "graph1","graph2"))
    X$cluster$Shiny <- Create_Shiny_Cluster(X)
    class(X$cluster) <- "Cluster_Object"
    return(X)
    }

#   ____________________________________________________________________________
#   2 Clustering Methods                                                    ####
##  ............................................................................
##  A k-nearest neighbors                                                   ####
#' k-nearest neighbors clustering
#'
#' @param X MCXpress object containing MCXmca object
#' @param mutual logical specifying if the connection between neighbors must be mutual
#' @param k maximum order of neighbors
#' @param rm minimum cluster size, all cluster size below will be designated as Cluster_0
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_knn <- function(X, mutual = TRUE, k = 3, rm = 0) {
    n <- nrow(X$Dim_Red$Cell2Cell_Distance)
    A <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        d <- sort(X$Dim_Red$Cell2Cell_Distance[i, ])
        A[i, X$Dim_Red$Cell2Cell_Distance[i, ] <= d[k +
            1]] <- 1
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
        names)%>% set_colnames(c("Cluster","Sample"))
    X$cluster$nClusters <- cluster_info$no
    g <- Calculate_Cluster_Centroids(X)
    return(g)
}
##  ............................................................................
##  B Maximum Distance                                                      ####
#' Title
#'
#' @param X MCXpress object containing a MCXmca object
#' @param mutual logical specifying if the connection between neighbors must be mutual
#' @param max.distance
#' @param rm
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_max_distance <- function(X, mutual = TRUE,
    max.distance = 0.5, rm = 0) {
    X$Dim_Red$Cell2Cell_Distance[row(X$Dim_Red$Cell2Cell_Distance) ==
        col(X$Dim_Red$Cell2Cell_Distance)] <- 0
    g <- graph.adjacency(X$Dim_Red$Cell2Cell_Distance,
        weighted = TRUE, mode = "undirected")
    g <- delete_edges(g, which(E(g)$weight > max.distance))
    g <- delete_edges(g, which(E(g)$weight == 0))
    cluster_info <- g %>% clusters
        X$cluster$Cluster_Quali <- paste0("Cluster", cluster_info$membership) %>% tibble(cluster_info$membership %>%
        names)%>% set_colnames(c("Cluster","Sample"))
    X$cluster$nClusters <- cluster_info$no
    g <- Calculate_Cluster_Centroids(X)
    return(g)
}
##  ............................................................................
##  C Percentage                                                            ####

#' Title
#'
#' @param X MCXpress object containing a MCXmca object
#' @param mutual
#' @param shortest.rank.percent
#' @param rm
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_percentage <- function(X, mutual = TRUE, shortest.rank.percent = 10,
    rm = 0) {
    X$Dim_Red$Cell2Cell_Distance[row(X$Dim_Red$Cell2Cell_Distance) ==
        col(X$Dim_Red$Cell2Cell_Distance)] <- 0
    g <- graph.adjacency(X$Dim_Red$Cell2Cell_Distance,
        weighted = TRUE, mode = "undirected")
    g <- delete_edges(g, which(E(g)$weight > sort(E(g)$weight)[shortest.rank.percent/100 *
        length(E(g)$weight)]))
    g <- delete_edges(g, which(E(g)$weight == 0))
    cluster_info <- g %>% clusters
        X$cluster$Cluster_Quali <- paste0("Cluster", cluster_info$membership) %>% tibble(cluster_info$membership %>%
        names)%>% set_colnames(c("Cluster","Sample"))
    X$cluster$nClusters <- cluster_info$no
    g <- Calculate_Cluster_Centroids(X)
    return(g)
}

##  ............................................................................
##  D Supervised Clustering                                                 ####
#' Supervised Clustering
#'
#' @param X MCXpress object containing Dim_Red object
#' @param Y A vector with Sample names as name and the coresponding cluster as value.
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_supervised<- function(X, Y){

    X$cluster$Cluster_Quali        <- Y

    names(X$cluster$Cluster_Quali) <- names(Y)
    X$cluster$Cluster_Quali        <- tibble(names(X$cluster$Cluster_Quali),
        X$cluster$Cluster_Quali) %>% set_colnames(c("Sample","Cluster"))
    X$cluster$nClusters            <- X$cluster$Cluster_Quali$Cluster %>%  unique %>%  length
    X                              <- Calculate_Cluster_Centroids(X)
    return(X)
}

##  ............................................................................
##  E K-Means Clustering                                                    ####
#' Title
#'
#' @param X MCXpress object containing a MCXmca object
#' @param nCluster
#' @param maxIter
#' @param nstart
#'
#' @returnMCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_kmeans <- function(X, k=2, maxIter = 10, nstart = 50){
  cat("Performing Kmeans Clustering with ", k, " Cluster")
    Distance <- X$Dim_Red$Cell2Cell_Distance
    Cluster <- Distance %>% kmeans(centers = k,
        iter.max = maxIter, nstart = nstart)
    X$cluster$Cluster_Quali <- paste0("Cluster", Cluster$cluster)
    names(X$cluster$Cluster_Quali) <- rownames(Distance)
    X$cluster$Cluster_Quali <- tibble(names(X$cluster$Cluster_Quali),
        X$cluster$Cluster_Quali) %>% set_colnames(c("Sample",
        "Cluster"))
    X$cluster$nClusters <- k
    X <- Calculate_Cluster_Centroids(X)
    return(X)
}
##  ............................................................................
##  F Hierarchical Clustering                                               ####
#' Title
#'
#' @param X MCXpress object containing a MCXmca object
#' @param method hclust method ("average", "ward", etc..)
#' @param k integer indicating the number of cluster to obtain
#' @param h numeric indicating the heights where the tree should be cut
#' @param members
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_hclust <- function(X, method="average", k=NULL, h=NULL, members=NULL) {
    Distance <- X$Dim_Red$Cell2Cell_Distance %>%  as.dist
    Cluster <- Distance %>% hclust(method=method, members=members)
    Cluster<-Cluster %>% cutree(k = k, h = h)
    X$cluster$Cluster_Quali <- tibble(paste0("Cluster", Cluster),(Cluster %>% names)) %>% set_names(c("Cluster", "Sample"))
    X$cluster$nClusters <- Cluster %>% unique %>%  length
    X <- Calculate_Cluster_Centroids(X)
    return(X)
}


##  ............................................................................
##  G Hierarchical Clustering                                               ####
#' Title
#'
#' @param X MCXpress object containing a MCXmca object
#' @param method hclust method ("average", "ward", etc..)
#' @param k integer indicating the number of cluster to obtain
#' @param h numeric indicating the heights where the tree should be cut
#' @param members
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' @export
cluster_k_medoids <- function(X, k = 2) {
    Distance <- X$Dim_Red$Cell2Cell_Distance
    Cluster <- Distance %>% as.dist %>%  pam(k=k, cluster.only = T)
    X$cluster$Cluster_Quali <- tibble(paste0("Cluster", Cluster),(Cluster %>% names)) %>% set_names(c("Cluster", "Sample"))
    X$cluster$nClusters <- Cluster %>% unique %>%  length
    X <- Calculate_Cluster_Centroids(X)
    return(X)
}



##  ............................................................................
#   ____________________________________________________________________________
#   3 Misc Cluster                                                          ####
##  ............................................................................
##  A Remove Small Clusters                                                 ####
#' Title
#'
#' @param g
#' @param min
#'
#' @return
#' @export
#'
#' @examples
removeSmallClusters = function(g, min = 1) {
    nClusters       <- clusters(g)$no
    cSizes          <- clusters(g)$csize
    smallClustersId <- which(cSizes <= min)
    cellsToRemoveId <- which(clusters(g)$membership %in% smallClustersId)
    newG            <- delete_vertices(g, cellsToRemoveId)
    return(newG)
}
##  ............................................................................
