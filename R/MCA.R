##MCA
# wants <- c("survival")
# has   <- wants %in% rownames(installed.packages())
# if(any(!has)) install.packages(wants[!has])
####Library
library(diptest)
library(abind)
library(fields)
library(igraph)
library(ade4)
library(reshape2)
library(fpc)

select <- dplyr::select
transmute <- dplyr::transmute



Initialise_MCXpress <- function(X, min_reads = 1) {
  if (X %>% is.matrix){
    if((X %>% rownames %>% is.null)|(X %>%  colnames %>% is.null)){errormessage<-"Gene name should be the rownames of the matrix and Sample name the column name"
    stop(errormessage)} else{
  MCXpress <- list()
  X <- as.matrix(unique(X[apply(X, 1, var) > 0, ]))
  X <- X[str_length(X %>% rownames) > 0, ]
  X <- X[rowSums(X) > min_reads, ]
  MCXpress$ExpressionMatrix<-X
  class(MCXpress)<- "MCXpress_object"
  return(MCXpress)}}else{errormessage<-"The input is not a matrix"
  stop(errormessage)}
}

###Gene Selection
Select_Most_Variable_Gene <- function(X, ngenes) {
  if(X %>%  class %>%  equals("MCXpress_object")){
  exp_matrix<-X$ExpressionMatrix
  means <- exp_matrix %>%  rowMeans
  vars  <- apply(exp_matrix, 1, var)
  cv2   <- vars / means ^ 2
  minMeanForFit <- unname(quantile(means[which(cv2 > .3)], 0.95))
  useForFit <- means >= minMeanForFit
  fit <-
    glmgam.fit(cbind(a0 = 1, a1tilde = 1 / means[useForFit]), cv2[useForFit])
  a0  <- unname(fit$coefficients["a0"])
  a1  <- unname(fit$coefficients["a1tilde"])
  fit$coefficients
  xg  <-
    exp(seq(min(log(means[means > 0])), max(log(means)), length.out = 1000))
  vfit  <- a1 / xg + a0
  df  <- ncol(exp_matrix) - 1
  afit  <- a1 / means + a0
  varFitRatio <- vars / (afit * means ^ 2)
  varorder  <- order(varFitRatio, decreasing = T)
  oed <- exp_matrix[varorder, ]
  X$ExpressionMatrix <- exp_matrix[varorder[1:ngenes],]
  return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
  stop(errormessage)}
}


Select_Highest_Mean <- function(X, ngenes) {
  if(X %>%  class %>%  equals("MCXpress_object")){
  exp_matrix<- X$ExpressionMatrix
  High_Means <-
    (exp_matrix %>%  rowMeans %>% sort(decreasing = TRUE) %>%  names %>% na.omit)
  High_Means <-  High_Means[str_length(High_Means) > 0]
  High_Means <-  High_Means[1:ngenes]
  X$ExpressionMatrix<-exp_matrix[High_Means, ]
  return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
  stop(errormessage)}
}

Select_Diptest <- function(X, pval = 0.1) {
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<- X$ExpressionMatrix
  dipPval <-
    as.numeric(lapply(1:dim(exp_matrix)[1], function(x)
      dip.test(as.numeric(exp_matrix[x, ]))$p.value))
  dipPval <- matrix(c(dipPval, seq(1:length(dipPval))), ncol = 2)
  dipPvalOrderedIndex <- dipPval %>% order
  dipPvalSigIndex <-
    dipPvalOrderedIndex[dipPval[dipPvalOrderedIndex] < pval]
  X$ExpressionMatrix<- exp_matrix[dipPvalSigIndex, ]
  return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
  stop(errormessage)}
}

#Gene Discretization Methods

Discretisation_Bsplines <- function(X, nbins = 2) {
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<-X$ExpressionMatrix
  nc = ncol(exp_matrix)
  ng = nrow(exp_matrix)
  discreteMatrix = array(0, dim = c(ng, nc, nbins))
  for (counterGene in seq(ng)) {
    j = bs(
      scale(exp_matrix[counterGene, ]),
      df = nbins,
      degree = 2,
      intercept = T
    )
    for (counterCondition in seq(nc)) {
      if (nbins == 2) {
        newJ <- cbind(rowSums(j[, 1:2]), j[, 3:dim(j)[2]])
        discreteMatrix[counterGene, counterCondition, ] = newJ[counterCondition, ]
      }
      else {
        discreteMatrix[counterGene, counterCondition, ] = j[counterCondition, ]
      }
    }
  }
  rownames(discreteMatrix) <- rownames(exp_matrix)
  colnames(discreteMatrix) <- colnames(exp_matrix)
  # Note: this is a fuzzy coded expression matrix
  Disjunctive_Matrix <- discreteMatrix %>% Create_Disjunctive_Matrix
  X$Disjunctive_Matrix<-Disjunctive_Matrix
  return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
  stop(errormessage)}
}

Discretisation_Range_01 <- function(X) {
  exp_matrix<-X$ExpressionMatrix
  Scale01 <-
    (apply(
      X = exp_matrix,
      1,
      FUN = function(x) {
        (x - min(x)) / (max(x) - min(x))
      }
    )) #Scale from 0 to 1
  Scale01bin1 <-
    t(Scale01) #Put in two different bins (bin1=probability of being expressed bin2=not expressed)
  Scale01bin2 <- t(1 - Scale01)
  discreteMatrix <-
    (abind(Scale01bin1, Scale01bin2, along = 3)) #create 3 dim array dim1=genes dim2=cells dim3=bins
  rownames(discreteMatrix) <-
    rownames(exp_matrix) #assign the coresponding row and colnames
  colnames(discreteMatrix) <- colnames(exp_matrix)
  Disjunctive_Matrix <- discreteMatrix %>% Create_Disjunctive_Matrix
  X$Disjunctive_Matrix<-Disjunctive_Matrix
  return(X)
}
#Dimension reduction methods

Create_Disjunctive_Matrix <- function(ExpressionMatrix) {
  cat('Creating Disjunctive Matrix...',  "\n")
  MCA_Matrix <- NULL
  for (i in 1:dim(ExpressionMatrix)[3]) {
    Bin_MCA <- ExpressionMatrix[, , i]
    rownames(Bin_MCA) <-
      paste (rownames(Bin_MCA), paste("-bin", i, sep = ""), sep = "")
    MCA_Matrix <- rbind(MCA_Matrix, Bin_MCA)
  }
  MCA_Matrix <- MCA_Matrix[order(rownames(MCA_Matrix)), ] %>% t()
  return(MCA_Matrix)
}


Dimension_reduction_MCA <- function(X, Dim = 5) {
  cat('Beginning MCA...\n')
  dis_matrix<- X$Disjunctive_Matrix
  MCA <- dis_matrix %>% as.data.frame(row.names = dis_matrix %>% rownames)
  MCA <- prep.fuzzy.var(MCA, rep(2, (MCA %>% ncol) / 2))
  MCA_results <-
    dudi.fca(MCA, nf = Dim, scannf = FALSE)
  Component <- paste0("PC", 1:Dim)
  X$Dim_Red$Cells_Standard  <-
    MCA_results$li %>%  set_colnames(Component) %>%  as.data.frame(row.names = MCA_results$li %>% rownames)
  X$Dim_Red$Cells_Principal <-
    MCA_results$l1 %>%  set_colnames(Component) %>%  as.data.frame(row.names = MCA_results$l1 %>% rownames)
  X$Dim_Red$Genes_Standard  <-
    MCA_results$co %>%  set_colnames(Component) %>%  as.data.frame(row.names = MCA_results$co %>% rownames)
  X$Dim_Red$Genes_Principal <-
    MCA_results$c1 %>%  set_colnames(Component) %>%  as.data.frame(row.names = MCA_results$c1 %>% rownames)
  #Distance Calculation
  cat('Calculating Cell to Cell Distance...\n')
  X$Dim_Red$Cell2Cell_Distance <-
    rdist(X$Dim_Red$Cells_Standard[, 1:Dim])
  X$Dim_Red$Cell2Cell_Distance %<>%  set_colnames(X$Dim_Red$Cells_Principal %>%  rownames)
  X$Dim_Red$Cell2Cell_Distance %<>%  set_rownames(X$Dim_Red$Cells_Principal %>%  rownames)
  X$Dim_Red$Cell2Gene_Distance <- rdist(X$Dim_Red$Cells_Principal, X$Dim_Red$Genes_Standard) %>%  set_colnames(X$Dim_Red$Genes_Standard %>% rownames) %>%  set_rownames(X$Dim_Red$Cells_Principal %>% rownames)

  #End Distance Calculation
  X$Dim_Red$Eigen_Value     <-
    MCA_results$eig %>% set_names(Component)
  X$Dim_Red$GraphEigen <-
    ((MCA_results$eig * 100) / (MCA_results$eig %>% sum)) %>%  as.data.frame() %>%  mutate(PC =
                                                                                             (paste0("PC", c(
                                                                                               1:length(MCA_results$eig)
                                                                                             )))) %>% mutate(Cumul = ((MCA_results$eig * 100) / (MCA_results$eig %>% sum)) %>% cumsum)
  colnames(X$Dim_Red$GraphEigen)[1] <- "EigenValue"
  X$Dim_Red$GraphEigen <- X$Dim_Red$GraphEigen[c(2, 1, 3)]
  X$Dim_Red$GraphEigen <-
    X$Dim_Red$GraphEigen[1:(X$Dim_Red$Cells_Standard %>%  nrow), ] %>% na.omit
  X$Dim_Red$Cumul <-
    X$Dim_Red$GraphEigen %>%  gather(EigenValue, Cumul, key = "Type", value = "Value")
  X$Dim_Red$GraphEigen <-
    X$Dim_Red$GraphEigen %>%  ggplot(aes(x = PC, y = Cumul, text = PC)) + geom_bar(stat =
                                                                                     "identity") + scale_x_discrete(limits = X$Dim_Red$GraphEigen$PC)
  X$Dim_Red$Methods <- "MCA"

  #Calculate Wilcoxon for Cell2CellDistance
  cat('Calculating Distance Wilcoxon...\n')
  Wilcox <- vector(length = (X$Dim_Red$Cells_Standard %>% ncol) - 1)
  for (i in 1:((X$Dim_Red$Cells_Standard %>% ncol) - 1)) {
    Wil <-
      wilcox.test(rdist((X$Dim_Red$Cells_Standard)[, 1:i]) %>% as.vector ,
                  rdist((X$Dim_Red$Cells_Standard)[, 1:(i + 1)]) %>% as.vector,
                  alternative = "less")
    Wilcox[i] <- Wil$p.value
    names(Wilcox)[i] <- paste0("PC", i, "-PC", i + 1)
  }
  Wilcoxon_df<-tibble(Wilcox %>%  names, Wilcox) %>%  set_colnames(c("Comp","p_value"))
  X$Dim_Red$Wilcoxon <- Wilcox
  #End Calculate Wilcoxon for Cell2CellDistance

  #Calculate Correlation PC and Genes
  X$Dim_Red$PC_Gene_Cor <-
    cor(dis_matrix, X$Dim_Red$Cells_Principal) %>% as.data.frame(row.names = X$Dim_Red$PC_Gene_Cor %>% rownames) %>% rownames_to_column(var =
                                                                                                                                        "Genes") %>%  as_tibble() %>%  gather(-Genes, key = "Component", value= "Cor")
  #End Calculate Correlation PC and Genes

  X$Dim_Red$Graph <-
    X$Dim_Red$Cells_Principal %>% ggplot(aes(
      x = PC1,
      y = PC2,
      text = row.names.data.frame(X$Dim_Red$Cells_Principal)
    )) + geom_point() + theme_light()
  X$Dim_Red$Shiny <- X %>% Create_Shiny_Dim_Red
  class(X$Dim_Red) <- "Dim_Red_Object"
  class(X) <- "MCXpress_object"
  cat('MCA is finished \n')
  return(X)
}


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
  names(X$cluster$Cluster_Quali)  <-  cluster_info %$% names(membership)
  X$cluster$Cluster_Quali <-
    tibble(names(X$cluster$Cluster_Quali), X$cluster$Cluster_Quali) %>%  set_colnames(c("Sample", "Cluster"))
  Cell_Coord_Cluster <-
    left_join(by = 'Sample',
              X$cluster$Cluster_Quali,
              Cells_Coord %>%  rownames_to_column(var = 'Sample')) %>%  select(-Sample)
  Coord_Centroids <-
    Cell_Coord_Cluster %>%  group_by(Cluster) %>%  summarise_each(funs(mean))
  X$cluster$nClusters <- cluster_info$no
  X$cluster$Coord_Centroids <- Coord_Centroids
  X$cluster$Distance <-
    Calculate_Cluster_Cluster_Distance(X$cluster$Coord_Centroids)
  X <- Calculate_Distance_Centroids_to_Gene(X)
  X$cluster$Graph1 <-
    ggplot(X$Dim_Red$Cells_Standard, aes(x = PC1, y = PC2)) + geom_point(aes(colour =
                                                                               X$cluster$Cluster_Quali %>% as.vector())) + theme_bw() + guides(colour =
                                                                                                                                                 guide_legend(title = "Cluster")) + ggtitle("Cells Low Dimensional Space") +
    theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Graph2 <-
    ggplot(X$Dim_Red$Genes_Standard, aes(x = PC1, y = PC2)) + geom_point(alpha =
                                                                           0.5) + theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(
                                                                             data = X$cluster$Coord_Centroids,
                                                                             aes(x = PC1, y = PC2, colour = Cluster),
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
    ggplot(X$Dim_Red$Cells_Standard, aes(x = PC1, y = PC2)) + geom_point(aes(colour =
                                                                               X$cluster$Cluster_Quali %>% as.vector())) + theme_bw() + guides(colour =
                                                                                                                                                 guide_legend(title = "Cluster")) + ggtitle("Cells Low Dimensional Space") +
    theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Graph2 <-
    ggplot(X$Dim_Red$Genes_Standard, aes(x = PC1, y = PC2)) + geom_point(alpha =
                                                                           0.5) + theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(
                                                                             data = X$cluster$Coord_Centroids,
                                                                             aes(x = PC1, y = PC2, colour = Cluster),
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
    ggplot(X$Dim_Red$Cells_Standard, aes(x = PC1, y = PC2)) + geom_point(aes(colour =
                                                                               X$cluster$Cluster_Quali %>% as.vector())) + theme_bw() + guides(colour =
                                                                                                                                      guide_legend(title = "Cluster")) + ggtitle("Cells Low Dimensional Space") +
    theme(legend.text = element_text(colour = "black", size = 8))
  X$cluster$Graph2 <-
    ggplot(X$Dim_Red$Genes_Standard, aes(x = PC1, y = PC2)) + geom_point(alpha =
                                                                           0.5) + theme_bw() + ggtitle("Genes Low Dimensional Space") + geom_point(
                                                                             data = X$cluster$Coord_Centroids,
                                                                             aes(x = PC1, y = PC2, colour = Cluster),
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
  A <- X$cluster$Coord_Centroids %>% select(contains('PC'))
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
  ClusterCoord_PC <- ClusterCoord %>%  select(contains('PC'))
  Distance <- rdist(ClusterCoord_PC)
  colnames(Distance) <- ClusterCoord$Cluster
  rownames(Distance) <- ClusterCoord$Cluster
  return(Distance)
}

Functional_Analysis_GSEA <-
  function(X,
           GMTfile,
           nperm = 1000,
           maxSize = 1000,
           minSize = 0,
           nproc = 1) {
        PCResults<- tibble()
    All_PCRanking<- tibble()
    for(i in 1:(min(X$Dim_Red$Eigen_Value%>% length,10)) ){
      cat(paste0("Processing PC", i, '...\n'))
      PCRanking<- X$Dim_Red$PC_Gene_Cor %>%
        filter(Component==paste0("PC",i)) %>%
        separate(col= Genes, into = c("Genes","bin"),sep = "-bin", convert = TRUE) %>%
        filter(bin==1) %>%
        mutate(AbsCor = Cor %>%  abs) %>%
        arrange(Cor) %>%
        mutate(Ranking=rank(Cor))
      GSEA_PCRank<-PCRanking$Ranking %>% as.matrix %>%  as.vector() %>%   set_names(value = PCRanking$Genes)
      GSEAResults <-
        fgsea(
          pathways = GMTfile,
          stats = GSEA_PCRank,
          nperm = nperm,
          maxSize = maxSize,
          minSize = minSize,
          nproc = nproc,
          BPPARAM = SerialParam(),
          gseaParam = 1
        )
      GSEAResults %<>%  as_tibble()%>% mutate(PC = rep(paste0("PC",i), GSEAResults %>%  nrow()))
      All_PCRanking <- bind_rows(All_PCRanking, PCRanking)
      PCResults <- bind_rows(PCResults, GSEAResults)
    }
    X$Functionnal_Analysis$GSEA_Results_PC<- PCResults[, c(ncol(PCResults), 1:(ncol(PCResults) - 1))]
    X$Functionnal_Analysis$GSEA_Results_PC[c("ES", "NES", "padj","pval")]<- X$Functionnal_Analysis$GSEA_Results_PC[c("ES", "NES", "padj","pval")] %>% format(scientific = TRUE, digits = 2)
    X$Functionnal_Analysis$RankingPC<- All_PCRanking %>%  select(-bin)
    cat("Beginning Gene Set Enrichment Analysis For Cluster...\n")
    Results <- tibble()
    All_Ranking <- tibble()
    Gene_All_Clusters_Distance <-
      X$cluster$Gene_Cluster_Distance %>% separate(col = Genes,
                                                   sep = '-bin',
                                                   into = c('Genes', 'bin')) %>% filter(bin == 1) %>% dplyr::select(-bin) %>%  gather(contains("Cluster"), key = Cluster, value =
                                                                                                                                        Distance)
    ID <- Gene_All_Clusters_Distance$Cluster %>% unique()
    for (i in 1:X$cluster$nClusters) {
      cat(paste0("Processing Cluster", i, '...\n'))
      Ranking <- NULL
      Gene_Cluster_Distance <-
        Gene_All_Clusters_Distance %>%  filter(Cluster == ID[i]) %>%  arrange(desc(Distance))
      Ranking_Table <-
        Gene_Cluster_Distance %>% select(Cluster, Genes, Distance) %>% mutate(Ranking =
                                                                                rank(Distance, ties.method = "random"))
      Ranking<-Ranking_Table$Ranking %>% rev
      names(Ranking) <- Ranking_Table$Genes
      GSEAResults <-
        fgsea(
          GMTfile,
          Ranking,
          nperm = nperm,
          maxSize = maxSize,
          minSize = minSize,
          nproc = nproc,
          BPPARAM = SerialParam(),
          gseaParam = 1
        )
      GSEAResults <-
        GSEAResults %>%  as_tibble()  %>% mutate(Cluster = rep(ID[i], GSEAResults %>%  nrow()))
      All_Ranking <- bind_rows(All_Ranking, Ranking_Table)
      Results <- bind_rows(Results, GSEAResults)
    }

    X$Functionnal_Analysis$GSEA_Results <-
      Results[, c(ncol(Results), 1:(ncol(Results) - 1))]
    X$Functionnal_Analysis$GSEA_Results[c("ES", "NES", "padj","pval")]<- X$Functionnal_Analysis$GSEA_Results[c("ES", "NES", "padj","pval")] %>% format(scientific = TRUE, digits = 2)
    X$Functionnal_Analysis$Ranking  <-  All_Ranking %>% mutate(Ranking=Ranking %>% rev)
    X$Functionnal_Analysis$Distance_Var <-
      X$cluster$Gene_Cluster_Distance %>% gather(-Genes, key = "Cluster", value = "Distance") %>% group_by(Genes) %>%  summarise(Dis_Var =
                                                                                                                                   var(Distance)) %>% mutate(Rank_Var = -Dis_Var %>% rank)

    X$Functionnal_Analysis$Shiny  <-Create_Shiny_Functionnal_Analysis(X)
    A<-X$Functionnal_Analysis$Ranking %>%  select(-Distance)
    B<-X$Functionnal_Analysis$RankingPC %>%  select(-Cor, -AbsCor) %>%  extract(,c(2,1,3))
    colnames(A)[1]<-"Group"
    colnames(B)[1]<-"Group"
    Grouped<-bind_rows(A,B)
    X$Functionnal_Analysis$Grouped<-Grouped
    X$Functionnal_Analysis$Info <-
      paste0("number of permutation:",
             nperm,
             "maxSize:",
             maxSize,
             "minSize:",
             minSize)
    cat(paste0('DONE\n'))
    class(X$Functionnal_Analysis) <- "FA_Object"
    return(X)
  }


Reactome_Category_Generator <- function(Info_File, Hierarchy_File) {
  colnames(Info_File)[1] = 'value'
  MAS <-  Hierarchy_File[grep("HSA", Hierarchy_File$X1), ]
  A   <-  !(Hierarchy_File$X1 %in% Hierarchy_File$X2)
  B   <-  Hierarchy_File$X1[A]
  C   <-  B[grep("HSA", B)] %>%  unique %>%  as_tibble
  D   <-  left_join(C, Info_File, by = "value") %>%  select(1:2)
  X   <-  MAS$X2 %>%  unique
  W   <-  vector(length = X %>%  length)
  for (i in 1:(X %>% length)) {
    Y <- X[i]
    Flag <- Y %in% MAS$X2
    while (Flag == TRUE) {
      Y <- MAS$X1[grep(Y, MAS$X2)]
      Flag <- Y %in% MAS$X2
    }
    W[i] <- Y
  }
  Q <- W %>%  as_tibble()
  colnames(Q) <- 'value'
  Q <- Q %>%  add_column(X)
  XC <- left_join(Q, D) %>% select(-value) %>% arrange(X2)
  colnames(XC)[1] <- 'value'
  Standard <- bind_rows(XC, D)
  Final <-
    left_join(Standard, Info_File, by = 'value') %>%  select(-value)
  colnames(Final) <- c('Category', 'pathway', 'species')
  return(Final)
}

#Object Print Methods

print.Dim_Red_Object <- function(obj, ...) {
  cat(obj$Methods, 'Dimension Reduction Results', "\n", "\n")
  NAME <-
    c(
      "",
      "$Cells_Standard",
      "$Cells_Principal",
      "$Genes_Standard",
      "$Genes_Principal",
      "$Graph",
      "$Eigen_Value",
      "$Shiny"
    )
  DESCRIPTION <-
    c(
      "",
      "Cells raw coordinate",
      "Cells normalised coordinate",
      "Genes raw coordinate",
      "Genes normalised coordinate",
      "Plot of Cell space",
      "eigenvalue",
      "Interactive Plot"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}

print.Cluster_Object <- function(obj, ...) {
  cat('Clustering Results', "\n", "\n")
  NAME <-
    c(
      "",
      "$nClusters",
      "$Cluster_Quali",
      "$Closest_Cluster",
      "$Coord_Centroids",
      "$Graph1",
      "$Graph2",
      "$Shiny"
    )
  DESCRIPTION <-
    c(
      "",
      "Number of Cluster",
      "A vector indicating the cluster for each samples",
      "Table indicating for each genes the cluster that is expressing the most",
      "Coordinates of the Cluster Centroids",
      "Cluster Plot in Cell Space",
      "Centroids of cluster plotted in Genespace",
      "Interactive Plot"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}

print.FA_Object <- function(obj, ...) {
  cat('Functionnal Analysis Results', "\n", "\n")
  NAME <- c("", "$Ranking", "$GSEA_Results", "$Shiny")
  DESCRIPTION <-
    c(
      "",
      "Table with Gene Ranking for each cluster",
      "fgsea package Gene Set Enrichment Analysis Results for each cluster",
      "Interactive Plot"
    )
  tibble(NAME, DESCRIPTION) %>%  print.data.frame(row.names = F, right = F)
}
##A small tweak on fgsea plotEnrichment function so it is compatible with plotly
plotlyEnrichment<-function(pathway,stats, gseaParam=1)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  name<-statsAdj[as.vector(na.omit(match(pathway, names(statsAdj))))] %>% sort %>% names %>% rev
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL

  data_markers <- data.frame(Genes=rep(name, each=100),x=rep(pathway, 100) %>% sort, y=rep(seq(from=diff/2, to=-diff/2, length.out = 100),pathway %>% length))
  plot_ly(data = toPlot) %>%  add_lines( x=~x, y=~y, mode="lines", line = list(color = 'rgb(154, 240, 24)', width = 2), name= "Enrichment") %>%
    add_lines(x=~x ,y = ~min(bottoms), line = list(color = 'rgb(110, 193, 248)', width = 2, dash = 'dash'), name= "Lower limit", hoverinfo="text", text=~round(min(bottoms), digits = 4)) %>%
    add_lines(x=~x ,y = ~max(tops), line = list(color = 'rgb(250, 150, 10)', width = 2, dash = 'dash'), name="Upper limit",  hoverinfo="text", text=~round(max(tops), digits = 4)) %>%
    add_markers(data=data_markers , x=~x ,y =~y, marker=list(color = 'rgb(0, 10, 10)', size=1), name="Genes", hoverinfo="text", text=~paste0(Genes,'</br>',"Rank:", x), showlegend=FALSE)
}

Reactome_GridPlot <-
  function(X, Info_File, Hierarchy_File, GMTfile, p = 0.25) {
    GMT <- GMTfile %>%  attributes() %>%  as_tibble
    colnames(GMT) <- 'pathway'
    Category <- Reactome_Category_Generator(Info_File, Hierarchy_File)
    GMT_Category <- inner_join(GMT, Category, by = 'pathway')
    Reactome_Table <-
      inner_join(X$Functionnal_Analysis$GSEA_Results, GMT_Category, by = 'pathway')
    Reactome_Table_filtered <-
      Reactome_Table %>%  filter(!is.na(padj)) %>%  filter(padj < p) %>%  mutate(Padj =
                                                                                   format(padj, scientific = TRUE, digits = 2))
    Reactome_Grid <- NULL
    Reactome_Grid <-
      ggplot(Reactome_Table_filtered,
             aes(
               Category,
               ES,
               colour = NES,
               text = pathway,
               text2 = Padj
             )) + geom_jitter(alpha = 0.7, size = (round(log(
               Reactome_Table_filtered$size
             ))) / 3) + facet_grid(Cluster ~ .,
                                   scales = "fixed",
                                   shrink = FALSE,
                                   space = "free_x") + theme_light() + theme(panel.grid.major.y = element_line(colour =
                                                                                                                 "gray")) + theme(axis.text.x = element_text(
                                                                                                                   face = "bold",
                                                                                                                   angle = 60,
                                                                                                                   hjust = 1,
                                                                                                                   size = 8
                                                                                                                 )) + scale_colour_gradient2(
                                                                                                                   low = "blue",
                                                                                                                   mid = "black",
                                                                                                                   high = "red",
                                                                                                                   midpoint = 0
                                                                                                                 )
    Reactome_Grid <-
      Reactome_Grid + geom_vline(xintercept = seq(1.5, length(
        unique(Reactome_Table_filtered$Category)
      ), 1), color = "gray")
    Reactome_Grid <- Reactome_Grid
    return(Reactome_Grid)
  }
