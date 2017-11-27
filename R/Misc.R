#' Initialisation of the MCXpress Object
#'
#' Create a MCXpressObject from an Expression Matrix
#'
#' @param X DESCRIPTION.
#' @param min_reads DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' @export
#'
Initialise_MCXpress <- function(X)
{
  if (X %>% is.matrix)
  {
    if ((X %>% rownames %>% is.null) | (X %>% colnames %>%
                                        is.null))
    {
      errormessage <-
        "Gene name should be the rownames of the matrix and Sample name the column name"
      stop(errormessage)
    } else
    {
      MCXpress <- list()
      X <- as.matrix((X[apply(X, 1, var) > 0,]))
      X <- X[str_length(X %>% rownames) > 0,]

      X <- X %>% subset(X %>% rownames %>% duplicated %>% not)
      colnames(X) <- gsub(pattern = " ",
                          replacement = "_",
                          x = colnames(X))
      colnames(X) <- gsub(pattern = "-",
                          replacement = "_",
                          x = colnames(X))
      colnames(X) <- gsub(pattern = "\\.",
                          replacement = "_",
                          x = colnames(X))
      MCXpress$ExpressionMatrix <- X
      class(MCXpress) <- "MCXpress"
      return(MCXpress)
    }
  } else
  {
    errormessage <- "The input is not a matrix"
    stop(errormessage)
  }
}



Df_To_Mat <- function(x) {
  Dframe <-
    x[, -1] %>% as.matrix %>% set_rownames(x[, 1] %>% as.matrix %>%   as.vector())
  return(Dframe)
}

Mat_To_Df <- function(x, namecol = "Sample") {
  Mat <-
    x %>%  as.data.frame() %>% tibble::rownames_to_column(var = namecol) %>%  as_tibble
  return(Mat)
}

Vec_To_Df <- function(x, namecol = c("X1", "X2")) {
  Dframe <-
    x %>%  as.data.frame() %>% tibble::rownames_to_column() %>%  as_tibble
  return(Dframe)
}

Num_Axis_Broken_Stick <- function(X) {
  vari <- X$MCA$explained_eigen_variance$Explained_Variance
  vari <- vari / 100
  broken <- sum(1 / 1:(vari %>% length)) / (vari %>% length)
  (vari > broken) %>%  which %>%  length
}

Num_Axis_Eigen_Distance <- function(X) {
  eig <- X$MCA$eigen_value
  A <- (eig[1] - eig[length(eig)]) / (1 - eig %>% length)
  B <- eig[1] - A
  Affine <- function(x) {
    (A * x + B - eig) %>% which.max()
  }
  return(Affine(1:(eig %>% length)))
}


#' Finding Top Expressed Genes for Each Cluster
#'
#' Search for each cluster the Top n closest genes in the MCA components and visualise it on a heatmap.
#'
#' @param X An MCXpress Object
#' @param n Integer specifying the number of genes to select for each cluster, default set at 5.
#' @param plotly logical, if FALSE gives a static heatmap, if TRUE gives an interactive version with heatmaply.
#'
#' @return
#' A Heatmap Visualisation
#' @export
#'
#' @examples
#'MCX64553 <- Initialise_MCXpress(GSE64553)
#'MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#'MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#'MCX64553 <- MCA(MCX64553, Dim = 5)
#'MCX64553 <- cluster_kmeans(MCX64553, k=6)
#'MCX64553 %>% Heatmap_Cluster_SC(n = 5, plotly = F)
Heatmap_Cluster_SC <- function(x, n = 5, plotly = F) {
  if (X$cluster %>% is.null) {
    stop("Clustering must be done before using this function")
  }
  OrderingCluster <-
    X$cluster$cluster_distances %>%  as.dist %>%  hclust(method = "ward.D") %>%  use_series(order)
  (X$cluster$cluster_distances %>% rownames)[OrderingCluster]
  X$cluster$labels$Cluster <-
    factor(X$cluster$labels$Cluster,
           levels = (X$cluster$cluster_distances %>%  rownames)[OrderingCluster])
  Cell <- X$cluster$labels %>%  dplyr::arrange(Cluster)
  TableCell <- c(0, Cell$Cluster %>% table %>% cumsum)
  Lab <-
    ((TableCell[-1] + TableCell[-(TableCell %>% length)]) / 2) %>% add(0.5) %>%  round(digits = 0)
  Lab2 <- rep(NA, times = Cell$Sample %>% length)
  Lab2[Lab] <- Lab %>%  names
  CellOrder <- Cell$Sample
  Gathering <- X$cluster$gene_cluster_distances %>%
    separate(Genes, sep = "-bin", into = c("Genes", "bin")) %>%
    dplyr::filter(bin == 1) %>%
    select(-bin) %>%
    gather(value = "Distance" , key = "Cluster" , -Genes) %>%
    dplyr::arrange(Cluster, Distance) %>%
    dplyr::filter(Cluster != "Origin") %>%
    group_by(Cluster) %>%
    top_n(n,-Distance)
  Gathering$Cluster <-
    factor(
      Gathering$Cluster,
      levels = X$cluster$labels %>%  dplyr::arrange(Cluster) %>%  use_series(Cluster) %>%  unique
    )
  GenesOrder <-  Gathering %>%  dplyr::arrange(Cluster) %$% Genes
  row_side_colors <-
    Cell %>%  select(Cluster)     %>%  dplyr::arrange(Cluster)
  col_side_colors <-
    Gathering %>%  select(Cluster) %>% dplyr::arrange(Cluster)
  HeatMatrix <- X$ExpressionMatrix[GenesOrder, CellOrder] %>%  t
  colsep <- col_side_colors %>%  table %>% cumsum()
  rowsep <- row_side_colors %>%  table %>% cumsum()
  if (plotly)
  {
    heatmaply::heatmaply(
      x = HeatMatrix,
      dendrogram = "none",
      Rowv = F,
      Colv = F,
      color = heat.colors(100) %>%  rev,
      col_side_colors  = col_side_colors,
      row_side_colors  = row_side_colors,
      row_side_palette = Spectral,
      col_side_palette = Spectral,
      heatmap_layers = c(
        lapply((rowsep %>% max) - rowsep, function(i)
          geom_hline(yintercept = i + 0.5, color = "black")),
        lapply((colsep %>% max) - colsep, function(i)
          geom_vline(xintercept = i + 0.5, color = "black"))
      ),
      na.rm = F,
      showticklabels = c(TRUE, FALSE)
    ) %>%  layout(showlegend = FALSE)
  } else
  {
    heatmap.2(
      HeatMatrix,
      Rowv = "none",
      Colv = "none",
      dendrogram = "none",
      trace = "none",
      colsep = c(0, (col_side_colors %>% table %>% cumsum())),
      rowsep = c(0, (row_side_colors %>% table() %>%  cumsum)),
      sepcolor = "black",
      col = heat.colors(100) %>%  rev,
      RowSideColors = rep(
        rainbow(row_side_colors %>% table() %>%  length),
        times = row_side_colors %>% table()
      ),
      ColSideColors = rep(
        rainbow(col_side_colors %>% table() %>%  length),
        times = col_side_colors %>% table()
      ),
      srtCol = 45,
      cexRow = 2,
      labRow = Lab2,
      cexCol = 1,
      keysize = 1,
      margins = c(10, 10),
      density.info = "none",
      na.color = "gray"
    )
  }
}

#' Title
#'
#' @param X MCXpress Object after Clustering
#' @param n Number of Genes to Display per Cluster
#' @param plotly Logical indicating if interactive html visualisation should be used
#' @return
#'
#' @export
#'
#' @examples
#'MCX64553 <- Initialise_MCXpress(GSE64553)
#'MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#'MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#'MCX64553 <- MCA(MCX64553, Dim = 5)
#'MCX64553 <- cluster_kmeans(MCX64553, k=6)
#'MCX64553 %>% Heatmap_Cluster(n = 5, plotly = F)
Heatmap_Cluster <- function(X, n = 5, plotly = T) {
  if (X$cluster %>% is.null) {
    stop("Clustering must be done before using this function")
  }
  Gathering <- X$cluster$gene_cluster_distances %>%
    separate(Genes, sep = "-bin", into = c("Genes", "bin")) %>%
    dplyr::filter(bin == 1) %>%
    select(-bin) %>%
    gather(value = "Distance" , key = "Cluster" , -Genes) %>%
    dplyr::arrange(Cluster, Distance) %>%
    dplyr::filter(Cluster != "Origin") %>%
    group_by(Cluster) %>%
    top_n(n,-Distance) %>% ungroup
  GenesOrder <-  Gathering %>%  dplyr::arrange(Cluster) %$% Genes
  HeatMatrix <- X$ExpressionMatrix[GenesOrder,] %>%  t
  HeatDF <-
    HeatMatrix %>%  t %>%  Mat_To_Df("Genes") %>% gather("Cells", "Value",-Genes) %>%  inner_join(X$cluster$labels %>% dplyr::rename(Cells = Sample))
  HeatMat <-
    HeatDF %>%  group_by(Cluster, Genes) %>%  summarise(MEAN = mean(Value)) %>%  spread(value = MEAN, key = Cluster)
  HeatMat <- HeatMat %>% Df_To_Mat() %>% t
  HeatMat <- HeatMat[, GenesOrder]
  row_side_colors <-
    HeatMat %>%  rownames %>% factor %>%  as.data.frame()
  col_side_colors <-
    Gathering %>%  select(Cluster) %>% dplyr::arrange(Cluster)
  colsep <- col_side_colors %>%  table %>% cumsum()
  rowsep <- row_side_colors %>%  table %>% cumsum()
  if (plotly)
  {
    heatmaply::heatmaply(
      x = HeatMat,
      seriate = "none",
      dendrogram = "none",
      Rowv = T,
      Colv = T,
      color =heat.colors(100) %>%  rev,
      col_side_colors  = col_side_colors,
      row_side_colors  = row_side_colors,
      row_side_palette = heatmaply::Spectral,
      col_side_palette = heatmaply::Spectral,
      heatmap_layers = c(
        lapply((rowsep %>% max) - rowsep, function(i)
          geom_hline(yintercept = i + 0.5, color = "black")),
        lapply((colsep %>% max) - colsep, function(i)
          geom_vline(xintercept = i + 0.5, color = "black"))
      ),
      na.rm = F,
      showticklabels = c(TRUE, TRUE),
      plot_method = "ggplot",
    ) %>%  layout(showlegend = FALSE)
  } else
  {
    heatmap.2(
      HeatMat,
      Rowv = "none",
      Colv = "none",
      dendrogram = "none",
      trace = "none",
      colsep = c(0, (col_side_colors %>% table %>% cumsum)),
      rowsep = c(0, (row_side_colors %>% table %>%  cumsum)),
      sepcolor = "black",
      col = heat.colors(100) %>%  rev,
      RowSideColors =
        heatmaply::Spectral(HeatMat %>% nrow),
      ColSideColors = rep(
        heatmaply::Spectral(HeatMat %>%  nrow),
        each = n
      ),
      srtCol = 45,
      cexRow = 1,
      labRow = rowlab,
      cexCol = 1,
      keysize = 1,
      margins = c(6, 6),
      density.info = "none",
      na.color = "gray"
    )
  }
}





#' Title
#'
#' @param X
#' @param pval
#' @param es
#' @param nes
#' @param color
#' @param title
#' @param rmna
#' @param metrics
#' @param plotly
#'
#' @return
#' @export
#'
#' @examples
GSEA_Heatmap_Cluster <-
  function(X, pval = 0.05, es = 0 , nes = -20, color = cm.colors(100), title = "", rmna = T,nPath=5, metrics = "NES", plotly = F, margin=c(0,0,0,0), cexCol=1, cexRow=1) {
    DF <- lapply(c("padj", "ES", "NES"), function(val) {
      df <-
        X$GSEA$GSEA_Results[!(names(X$GSEA$GSEA_Results) == "Origin")] %>%
        purrr::map_df(
          .f = function(x) {
            x %>% extract2(val)
          }
        ) %>%
        mutate(Pathway = X$GSEA$Pathways) %>%
        dplyr::select(Pathway, everything())
      df2 <-
        df %>% gather_(key = "Cells", value = val, colnames(df)[-1])
    })
    DF <- DF[[1]] %>%  inner_join(DF[[2]], c("Pathway", "Cells")) %>%  inner_join(DF[[3]], c("Pathway", "Cells"))
    if (metrics == "ES") {
      DF$ES[(DF$padj > pval) | (DF$ES < es)] <- NA
      MAT <-
        DF %>%  dplyr::select(Pathway, Cells, ES) %>%  spread(Cells, ES)
    }
    if (metrics == "NES") {
      DF$NES[(DF$padj > pval) | (DF$NES < nes)] <- NA
      MAT <-
        DF %>%  dplyr::select(Pathway, Cells, NES) %>%  spread(Cells, NES)
    }
    MATPadj <-
      DF %>%  dplyr::select(Pathway, Cells, padj) %>%  spread(Cells, padj)
    InMat <-
      MAT  %>% dplyr::select(-Pathway) %>%  as.matrix.data.frame %>%  set_rownames(MAT$Pathway)
    if (rmna) {
      InMat <- InMat[(InMat %>% is.na %>% rowSums) != (InMat %>% ncol),]
    }
    InMat <- GSEA_select_variable(InMat, nPath = nPath)
    RowOrder <-
      InMat %>% is.na %>%  ifelse(0, InMat) %>%  dist %>%  hclust(method = "ward.D") %>% use_series(order)
    InMat <- InMat[RowOrder,]
    InMat <- InMat %>% t
    if(plotly){heatmaply::heatmaply(
      InMat,
      seriate = "none",
      row_dend_left = F,
      Rowv = F,
      Colv = F,
      grid_color = "black",
      dendrogram = "none",
      draw_cellnote = F,
      cexCol = cexCol,
      cexRow = cexRow,
      margins = margin,
      na.value = "black",
      colors = color,
      na.rm = FALSE,
      row_side_colors = InMat %>%  rownames %>%  factor %>%  data.frame(),
      row_side_palette = rainbow
    ) %>% layout(showlegend = FALSE)}
    else{
      heatmap.2(
        InMat,
        Rowv = "none",
        Colv = "none",
        dendrogram = "none",
        trace = "none",
        colsep = 0:(InMat %>%  ncol),
        rowsep = 0:(InMat %>%  nrow),
        sepcolor = "black",
        col = color,
        RowSideColors = rainbow(InMat %>%  nrow),
        srtCol = 45,
        cexRow = cexRow,
        cexCol = cexCol,
        keysize = 1,
        key.xlab = metrics,
        key.title = "Enrichment",
        margins = margin,
        density.info = "none",
        na.color = "gray",
        main = title
      )
    }

    }

#' Heatmap of GSEA at Single Cell Level
#'
#' @param x An MCXpress object
#' @param pval Numeric indicating the threshold pvalue
#' @param es  Numeric indicating enrichment score threshold
#' @param nes Numeric normalised enrichment score threshold
#' @param color color palette to use for the heatmap
#' @param title Plot Title
#' @param rmna Logical indicating if any empty column should be removed.
#' @param metrics NES or ES
#' @param plotly Logical Indicating if the static version or the interactive version should be plotted.
#'
#' @return
#' @export
#'
#' @examples
GSEA_Heatmap_SC <-
  function(X,
           pval = 0.05,
           es = 0 ,
           nes = -20,
           color = cm.colors(100),
           title = "",
           rmna = T,
           metrics = "NES",
           plotly = T,
           cexRow=1,
           cexCol=1,
           nPath=10,
           jaccard=T){
    cluster <- X$cluster$labels %>%  dplyr::rename(Cells=Sample)
    DF      <- create_gsea_matrix(X, pval=pval, es=es, nes=nes, metrics = metrics, SC=T)
    if(metrics=="NES"){MAT<- DF %>%  dplyr::select(Pathway, Cells, NES) %>%  spread(Cells, NES)}
    else{MAT <- DF %>%  dplyr::select(Pathway, Cells, ES) %>%  spread(Cells, ES)}
    DFORDER <- DF %>%  inner_join(cluster, by="Cells")%>% arrange(Cluster)
    newOrder <- DFORDER  %>%
      split(DFORDER$Cluster) %>%
      map(function(x) {
        clus_gsea_matrix <- x %>%
          select(Pathway, Cells, ES) %>%
          spread(Cells, ES)
        clus_gsea_matrix <-
          clus_gsea_matrix %>%  select(-Pathway) %>% as.matrix %>% apply(2, as.numeric) %>%  set_rownames(clus_gsea_matrix$Pathway)
        clus_gsea_matrix[clus_gsea_matrix %>% is.na()] <- 0
        clus_gsea_matrix[, clus_gsea_matrix %>% t %>%  dist("euclidean") %>%  hclust(method ="ward.D") %>% use_series(order)] %>% colnames
      }) %>% Reduce(c,.)

    InMat <-
      MAT  %>% dplyr::select(-Pathway) %>%  as.matrix.data.frame %>%  set_rownames(MAT$Pathway)
    InMat <- GSEA_remove_na(InMat, rmna=rmna)
    InMat <- GSEA_select_variable(InMat, nPath = nPath)
    if(jaccard){
    col_order <-   Jaccard(X$SC_GSEA$GMTfile[InMat %>%  rownames])%>% hclust("ward.D2") %>% use_series(order)
    }
    else
    {col_order <-
      InMat %>% is.na %>%  ifelse(0, InMat) %>%  dist %>%  hclust(method = "ward.D") %>% use_series(order)}
    Clust2 <- cluster %$% set_names(Cluster, Cells)
    Tab <- Clust2 %>%  table
    Tab <- Tab[Tab %>% names %>% gtools::mixedsort()]
    SepROW <- c(0, Tab %>% cumsum())
    LabROW1 <-
      ((SepROW[-(SepROW %>% length)] + Tab %>% cumsum()) / 2) %>% round(0) %>% magrittr::set_names(value = Tab %>% names)
    LabROW2 <- rep(NA,  X$ExpressionMatrix %>% ncol)
    LabROW2[LabROW1] <- LabROW1 %>% names
    row_side_color <-
      rep(Tab %>%  names, times = Tab) %>%  factor %>%  data.frame()

    InMat <- InMat %>% t
    InMat <- InMat[newOrder, col_order]
    labrowPlotly <- LabROW2 %>% subset(.,LabROW2 %>% is.na %>% not) %>% set_names((InMat %>% rownames)[LabROW2 %>% is.na %>% not])
    if (plotly) {
      heatmaply::heatmaply(
        InMat,
        row_dend_left = F,
        Rowv = F,
        Colv = F,
        dendrogram = "none",
        draw_cellnote = F,
        cexCol = cexCol,
        cexRow = cexRow,
        margins = c(200, 200, 50, 200),
        na.value = "black",
        colors = color,
        na.rm = FALSE,
        heatmap_layers = c(lapply(max(SepROW) - SepROW, function(i)
          geom_hline(yintercept = i + 0.5, color = "black")),
          scale_y_discrete(labels=LabROW2[LabROW2 %>% is.na %>% not()],breaks=(InMat %>% rownames())[LabROW2 %>% is.na %>% not()])),
        row_side_colors = row_side_color,
        row_side_palette = rainbow,
        showticklabels = c(T, T)
      ) %>%
        layout(showlegend = FALSE)
    }
    else{
      heatmap.2(
        InMat,
        Rowv = "none",
        Colv = "none",
        dendrogram = "none",
        trace = "none",
        colsep = 0:(InMat %>%  ncol),
        rowsep = SepROW,
        sepcolor = "black",
        col = color,
        labRow = LabROW2,
        RowSideColors = rep(rainbow(Tab %>%  length), times = Tab),
        srtCol = 45,
        cexRow = cexRow,
        cexCol = cexCol,
        keysize = 1,
        key.xlab = metrics,
        key.title = "Enrichment",
        margins = c(10, 10),
        density.info = "none",
        na.color = "gray",
        main = title
      )
    }
  }

create_gsea_matrix <- function(X, pval=0.05, es=-1, nes=0, metrics="NES", SC=T){
  list_gsea_df <- lapply(c("padj", "ES", "NES"), function(val) {
    if(SC){
      gsea_results <- X$SC_GSEA$GSEA_Results
      path <- X$SC_GSEA$Pathways
    }
    else{
      gsea_results <- X$GSEA$GSEA_Results
      path <- X$GSEA$Pathways
    }
    df <-
      gsea_results[!(names(gsea_results) == "Origin")] %>%
      purrr::map_df(
        .f = function(x) {
          x %>% extract2(val)
        }
      ) %>%
      mutate(Pathway = path) %>%
      dplyr::select(Pathway, everything())
    df <-
      df %>% gather_(key = "Cells", value = val, colnames(df)[-1])
  })
  list_gsea_df <- list_gsea_df[[1]] %>%  inner_join(list_gsea_df[[2]], c("Pathway", "Cells")) %>%  inner_join(list_gsea_df[[3]], c("Pathway", "Cells"))
  if (metrics == "ES") {
    list_gsea_df$ES[(list_gsea_df$padj > pval) | (list_gsea_df$ES < es)] <- NA
    }
  if (metrics == "NES") {
    list_gsea_df$NES[(list_gsea_df$padj > pval) | (list_gsea_df$NES < nes)] <- NA
    }
  return(list_gsea_df)
}


GSEA_Heatmap_SC2 <-
  function(X,
           pval = 0.05,
           es = 0 ,
           nes = -20,
           color = cm.colors(100),
           title = "",
           rmna = T,
           metrics = "NES",
           plotly = F,
           cexRow=1,
           cexCol=1,
           nPath=20){
    cluster <- X$cluster$labels %>%  dplyr::rename(Cells=Sample)
    DF      <- create_gsea_matrix(X, pval=pval, es=es, nes=nes, metrics = metrics, SC=T)
    if(metrics=="NES"){MAT<- DF %>%  dplyr::select(Pathway, Cells, NES) %>%  spread(Cells, NES)}
    else{MAT<- DF %>%  dplyr::select(Pathway, Cells, ES) %>%  spread(Cells, ES)}
    InMat <-
      MAT  %>% dplyr::select(-Pathway) %>%  as.matrix.data.frame %>%  set_rownames(MAT$Pathway)
    if (rmna) {
      InMat <- InMat[(InMat %>% is.na %>% rowSums) != (InMat %>% ncol),]
    }
    InMat <- InMat %>% t
    InMat[InMat %>% is.na] <- 0
    InMat <- InMat[,InMat %>%  apply(2, var) %>%  sort(T) %>%  head(nPath) %>%  names]
    A <- cluster %>%  arrange(Cells) %>%  use_series(Cluster) %>% unique %>%  length %>%  rainbow %>% set_names(cluster %>%  arrange(Cells) %>%  use_series(Cluster) %>% unique)
    B <- A[cluster %>%  arrange(Cells) %>%  use_series(Cluster)]
    if (plotly) {
      heatmaply::heatmaply(
        InMat,
        row_dend_left = T,
        Rowv = T,
        Colv = T,
        hclust_method="ward.D2",
        dendrogram = "both",
        draw_cellnote = F,
        margins = c(200, 200, 50, 200),
        na.value = "black",
        row_side_colors =cluster %>%  arrange(Cells) %>%  select(Cluster),
        colors = color,
        na.rm = F,
        showticklabels = c(T, F)
      ) %>%
        plotly::layout(showlegend = FALSE)
    }
    else{
      heatmap.2(
        InMat,
        Rowv = T,
        Colv = T,
        dendrogram = "both",
        hclustfun = function(x) hclust(x, method = "ward.D2"),
        trace = "none",
        labRow = "",
        col = color,
        srtCol = 45,
        RowSideColors = B,
        keysize = 1,
        key.xlab = metrics,
        key.title = "Enrichment",
        margins = c(10, 10),
        density.info = "none",
        na.color = "gray",
        main = title
      )
    }
  }

Jaccard <- function(GeneSet){sapply(X = GeneSet, FUN = function(Hall){GeneSet %>%  sapply(FUN= function(Hall2, Hall){
  1-((intersect(Hall,Hall2) %>% length)/(union(Hall,Hall2) %>% length))},Hall=Hall)}) %>%  as.dist
}

GSEA_select_variable <- function(x, nPath){
  y <- x
  y[y %>% is.na] <- 0
  z <- y %>%  apply(1, var) %>%  sort(T) %>%  head(nPath) %>%  names
  x <- x[z,]
}

GSEA_remove_na <- function(x,rmna=T){if(rmna){
  x <- x[(x %>% is.na %>% rowSums) != (x %>% ncol),]
} else{x}}

