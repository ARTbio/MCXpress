#   ____________________________________________________________________________
#   1 Cluster Statistics                                                    ####
#' Generic function to compute distances, centroids and statistics after statistics.
#'
#' @param X an MCXpress object with MCA object (after MCA step)
#' @return MCXpress object with MCA object and Clustering Object (after MCA step)
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
    X$Shiny <- Create_Dashboard2(X)
    class(X$cluster) <- "Cluster_Object"
    return(X)
    }

#   ____________________________________________________________________________
#   2 Clustering Methods                                                    ####

##  ............................................................................
##  A Supervised Clustering                                                 ####
#' Supervised Clustering on MCA
#'
#' @param X MCXpress object containing Dim_Red object
#' @param Y A vector with Sample names as name and the coresponding cluster as value.
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{Cluster_Quali}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{Gene_Cluster_Distance}{GSEA Parameter used}
#' \item{closest_Cluster}{Indicate the closest cluster for a gene}
#' \item{Coord_Centroids}{GSEA Parameter used}
#' \item{graph1}{Clustering Visualisation in the MCA cell space}
#' \item{graph2}{Visualisation of the centroids in the MCA gene space}
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
##  B K-Means Clustering                                                    ####
#' K-Means Clustering on MCA
#'
#' @param X MCXpress object containing a MCXmca object
#' @param k Number of cluster to compute
#' @param maxIter maximum number of iteration
#' @param nstart number of random sets of initial centers
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{Cluster_Quali}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{Gene_Cluster_Distance}{GSEA Parameter used}
#' \item{closest_Cluster}{Indicate the closest cluster for a gene}
#' \item{Coord_Centroids}{GSEA Parameter used}
#' \item{graph1}{Clustering Visualisation in the MCA cell space}
#' \item{graph2}{Visualisation of the centroids in the MCA gene space}
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
##  C Hierarchical Clustering                                               ####
#' Hierarchical Clustering on MCA
#'
#' @param X MCXpress object containing a MCXmca object
#' @param method hclust method ("average", "ward", etc..)
#' @param k integer indicating the number of cluster to obtain
#' @param h numeric indicating the heights where the tree should be cut
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{Cluster_Quali}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{Gene_Cluster_Distance}{GSEA Parameter used}
#' \item{closest_Cluster}{Indicate the closest cluster for a gene}
#' \item{Coord_Centroids}{GSEA Parameter used}
#' \item{graph1}{Clustering Visualisation in the MCA cell space}
#' \item{graph2}{Visualisation of the centroids in the MCA gene space}
#' @export
cluster_hclust <- function(X, method="average", k=NULL, h=NULL) {
    Distance <- X$Dim_Red$Cell2Cell_Distance %>%  as.dist
    Cluster <- Distance %>% hclust(method=method)
    Cluster<-Cluster %>% cutree(k = k, h = h)
    X$cluster$Cluster_Quali <- tibble(paste0("Cluster", Cluster),(Cluster %>% names)) %>% set_names(c("Cluster", "Sample"))
    X$cluster$nClusters <- Cluster %>% unique %>%  length
    X <- Calculate_Cluster_Centroids(X)
    return(X)
}


##  ............................................................................
## D K_Medoids Clustering                                               ####
#' K-Medoids Clusetering on the MCA Analysis
#'
#' @param X MCXpress object containing a MCXmca object
#' @param k integer indicating the number of cluster to obtain
#'
#' @return MCXpress object containing a MCXmca and MCXcluster object
#' \item{Cluster_Quali}{Clustering results}
#' \item{nClusters}{Number of Cluster}
#' \item{Gene_Cluster_Distance}{GSEA Parameter used}
#' \item{closest_Cluster}{Indicate the closest cluster for a gene}
#' \item{Coord_Centroids}{GSEA Parameter used}
#' \item{graph1}{Clustering Visualisation in the MCA cell space}
#' \item{graph2}{Visualisation of the centroids in the MCA gene space}
#' @export
cluster_k_medoids <- function(X, k = 2) {
    Distance <- X$Dim_Red$Cell2Cell_Distance
    Cluster <- Distance %>% as.dist %>%  pam(k=k, cluster.only = T)
    X$cluster$Cluster_Quali <- tibble(paste0("Cluster", Cluster),(Cluster %>% names)) %>% set_names(c("Cluster", "Sample"))
    X$cluster$nClusters <- Cluster %>% unique %>%  length
    X <- Calculate_Cluster_Centroids(X)
    return(X)
}

