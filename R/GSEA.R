#' Gene Set Enrichment Analysis of Clusters and Axis
#'
#' Performs a gene set enrichment analysis with the fgsea bioconductor package
#' on the cluster and the axis using the ranking of the distance between the clusters centroids and the genes.
#'
#' @param X A MCXPress object containing a cluster and Dimred object
#' @param GMTfile .gmt file containing a list of pathway with their represnting
#'   gene.
#' @param nperm the number of permutation.
#' @param minSize threshold for the minimum number of genes that must be in the pathway.
#' @param maxSize threshold for the maximum number of genes that must be in the pathway.
#' @param nproc number of processor to use for the enrichment algorithm.
#' @param nbin an integer indicating which bin to use for the enrichment analysis.
#' @param naxis number of axis to perform the enrichment analysis.
#' @param gseaParam an integer to weight the enrichment analysis.
#' @return Return a MCXpress object containing a GSEA object
#' \item{GSEA_Results}{List containing dataframes of the GSEA analysis for each Cluster}
#' \item{Ranking}{List containing vectors of gene ranking according to their distance in MCA for each Cluster}
#' \item{GSEA_Results_Axis}{List containing dataframes of the GSEA analysis for each Axis}
#' \item{Ranking_Axis}{List containing vectors of gene ranking according to their MCA cells coordinate and genes expression correlation for each Axis}
#' \item{GMTfile}{GMTfile used for GSEA}
#' \item{Pathway}{Name of the genesets contained in the gmtfile}
#' \item{gseaParam}{GSEA Parameter used}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553, Dim = 5)
#' MCX64553 <- cluster_kmeans(MCX64553, k=6)
#' MCX64553 <- GSEA(X = MCX64553, GMTfile = reactome_gmtfile, nperm = 10000)
#' @export
# ____________________________________________________________________________
# Gene set enrichment analysis ####
GSEA <- function(X, GMTfile, nperm = 1000, minSize = 15, maxSize = 500,
    nproc = 4, nbin = 1, naxis = 2, gseaParam = 0)
    {

    ## ............................................................................
    ## A Axis Correlation Gene Set Enrichment Analysis ####


    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ### . . . . . . . ..  a Filtering of bin ####
    df <- X$MCA$Axis_Gene_Cor[, 1:(naxis + 1)] %>% separate(col = Genes,
        into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
        filter(bin == 1) %>% select(-bin)

    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ### . . . . . . . ..  b Correlation Ranking Calculation ####
    cat("Calculating ranking of genes correlation for each axis \n")
    axis_rank <- df %>% select(-Genes) %>% purrr::map(function(x)
    {x %>% abs %>% set_names(df$Genes) %>% rank(ties.method = "random") %>%
            sort
    })

    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ### . . . . . . . ..  c Fast Gene Set Enrichment Analysis on
    ### each Axis ####

    cat("Beginning enrichment analysis for axis \n")
    pb <- txtProgressBar(width = 50, style=3, char = "+")
    axis_gsea <- axis_rank %>% purrr::map2(.y = seq(from=0, to=1, length.out = (axis_rank %>%  length)), .f = function(x, y)
        {
        setTxtProgressBar(pb, y)
        gsea <- fgsea(pathways = GMTfile, stats = x, nperm = nperm,
            maxSize = maxSize, minSize = minSize, nproc = nproc,
            BPPARAM = SerialParam(), gseaParam = gseaParam) %>%
            as_tibble
        ins <- gsea %>% magrittr::extract(, 2:5) %>% round(digits = 5)
        val <- gsea %>% magrittr::inset(, 2:5, value = ins)
        return(val)
    })
    close(pb)

    ## ............................................................................
    ## B Gene Set Enrichment Analysis on Cluster ####

    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ### . . . . . . . ..  a Filter Bin ####
    df2 <- X$cluster$gene_cluster_distances %>% tidyr::separate(col = Genes,
        into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
        dplyr::filter(bin %in% nbin)

    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ### . . . . . . . ..  b Calculate Rank for Cluster ####

    cat("\nCalculating ranking of genes for each clusters \n")
    cluster_rank <- df2 %>% dplyr::select(-Genes,-bin) %>% lapply(FUN = function(x)
    {
        1-(x %>% set_names(df2$Genes) %>% (function(x){2*(x-min(x))/(max(x)-min(x))}) %>% sort)
    })

    ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    ### . . . . . . . ..  c fgsea analysis for cluster ####

    cat("Beginning enrichment analysis for clusters\n\n")
    pb <- txtProgressBar(width = 50, style=3, char = "+")
    cluster_gsea <- mapply(x=cluster_rank, y = seq(from=0, to=1, length.out = (cluster_rank %>%  length)), FUN= function(x,y)
        {
        setTxtProgressBar(pb, y)
        gsea <- fgsea(pathways = GMTfile, stats = x, nperm = nperm,
            maxSize = maxSize, minSize = minSize, nproc = nproc,
          BPPARAM = SerialParam(), gseaParam = gseaParam) %>%
            as_tibble()
        ins <- gsea %>% magrittr::extract(, 2:5) %>% round(digits = 5)
        val <- gsea %>% magrittr::inset(, 2:5, value = ins)
        return(val)
    }, SIMPLIFY = F)
    close(pb)

    ## ............................................................................
    ## C GSEA finalisation ####

    cat(paste0("Creating Enrichment Analysis Object\n"))
    X$GSEA$GSEA_Results_Axis <- axis_gsea
    X$GSEA$RankingAxis <- axis_rank
    X$GSEA$GSEA_Results <- cluster_gsea
    X$GSEA$Ranking <- cluster_rank
    X$GSEA$GMTfile <- GMTfile
    X$GSEA$Pathways <- axis_gsea$Axis1$pathway
    X$GSEA$AllRanking <- axis_rank %>% append(cluster_rank)
    X$GSEA$gseaParam <- gseaParam
    X$Shiny <- create_dashboard3(X)
    cat(paste0("Enrichment Analysis Completed\n"))
    class(X$GSEA) <- "GSEA"
    return(X)
}



parallel_fgsea <- function(x, a, b, c, d, e, f){
  gsea <- tibble::as_tibble(fgsea::fgsea(stats = x, pathways = a, nperm = b,
                                         minSize =  c, maxSize = d, nproc = e,
                                         BPPARAM = BiocParallel::SerialParam(), gseaParam =f))
  ins <-  round(magrittr::extract(gsea, 2:5),digits = 5)
  val <-  magrittr::inset(gsea, 2:5, value = ins)
  return(val)
}


#'GSEA parallel computing version
#'
#'Compute GSEA on cluster, number of axis retained is the same as the one used for the clustering step. 
#'
#' @param X A MCXPress object containing a Dimred object
#' @param GMTfile named list of pathway with their representing gene. 
#' @param nperm the number of permutation.
#' @param minSize threshold for the minimum number of genes that must be in the pathway.
#' @param maxSize threshold for the maximum number of genes that must be in the pathway.
#' @param nproc number of processor to use for the enrichment algorithm.
#' @param nbin an integer indicating which bin to use for the enrichment analysis.
#' @param naxis number of axis to perform the enrichment analysis.
#' @param gseaParam an integer to weight the enrichment analysis.
#' @return Return a MCXpress object containing a GSEA object
#' \item{GSEA_Results}{List containing dataframes of the GSEA analysis for each Cluster}
#' \item{Ranking}{List containing vectors of gene ranking according to their distance in MCA for each Cluster}
#' \item{GSEA_Results_Axis}{List containing dataframes of the GSEA analysis for each Axis}
#' \item{Ranking_Axis}{List containing vectors of gene ranking according to their MCA cells coordinate and genes expression correlation for each Axis}
#' \item{GMTfile}{GMTfile used for GSEA}
#' \item{Pathway}{Name of the genesets contained in the gmtfile}
#' \item{gseaParam}{GSEA Parameter used}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553, Dim = 5)
#' MCX64553 <- cluster_kmeans(MCX64553, k=6)
#' MCX64553 <- GSEA(X = MCX64553, GMTfile = reactome_gmtfile, nperm = 10000)
#' @export
GSEAparall <- function(X, GMTfile, nperm = 1000, minSize = 15, maxSize = 500,
                 nproc = 4, nbin = 1, naxis = 2, gseaParam = 0)
{
  ## ............................................................................
  ## A Axis Correlation Gene Set Enrichment Analysis ####


  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  a Filtering of bin ####
  df <- X$MCA$Axis_Gene_Cor[, 1:(naxis + 1)] %>% separate(col = Genes,
                                                          into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
    filter(bin == 1) %>% select(-bin)

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  b Correlation Ranking Calculation ####
  cat("Calculating ranking of genes correlation for each axis \n")
  axis_rank <- df %>% select(-Genes) %>% purrr::map(function(x)
  {x %>% abs %>% set_names(df$Genes) %>% rank(ties.method = "random") %>%
      sort
  })

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  c Fast Gene Set Enrichment Analysis on
  ### each Axis ####

  cat("Beginning enrichment analysis for axis \n")
  pb <- txtProgressBar(width = 50, style=3, char = "+")
  axis_gsea <- axis_rank %>% purrr::map2(.y = seq(from=0, to=1, length.out = (axis_rank %>%  length)), .f = function(x, y)
  {
    setTxtProgressBar(pb, y)
    gsea <- fgsea(pathways = GMTfile, stats = x, nperm = nperm,
                  maxSize = maxSize, minSize = minSize, nproc = nproc,
                  BPPARAM = SerialParam(), gseaParam = gseaParam) %>%
      as_tibble
    ins <- gsea %>% magrittr::extract(, 2:5) %>% round(digits = 5)
    val <- gsea %>% magrittr::inset(, 2:5, value = ins)
    return(val)
  })
  close(pb)

  ## ............................................................................
  ## B Gene Set Enrichment Analysis on Cluster ####

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  a Filter Bin ####
  df2 <- X$cluster$gene_cluster_distances %>% tidyr::separate(col = Genes,
                                                              into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
    dplyr::filter(bin %in% nbin)

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  b Calculate Rank for Cluster ####
  Scaling2 <- function(x){2*(x-min(x))/(max(x)-min(x))}
  cat("\nCalculating ranking of genes for each clusters \n")
  # cluster_rank <- df2 %>% dplyr::select(-Genes,-bin) %>% bplapply(FUN = function(x,y){
  #   1-(sort(Scaling2(set_names(x,y))))
  # },BPPARAM = SnowParam(workers = nproc, tasks=nproc, progressbar = T), y=df2$Genes)
  cat("\nCalculating ranking of genes for each clusters \n")
  cluster_rank <- df2 %>% dplyr::select(-Genes,-bin) %>% lapply(FUN = function(x)
  {
    1-(x %>% set_names(df2$Genes) %>% (function(x){2*(x-min(x))/(max(x)-min(x))}) %>% sort)
  })
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  c fgsea analysis for cluster ####

  cat("Beginning enrichment analysis for clusters\n\n")
  cluster_gsea <-
    cluster_rank %>% BiocParallel::bplapply(
      FUN = parallel_fgsea,
      BPPARAM = SnowParam(workers = nproc, tasks=nproc, progressbar = T),
      a = GMTfile,
      b = nperm,
      c = minSize,
      d = maxSize,
      e = nproc,
      f = gseaParam
    )


  ## ............................................................................
  ## C GSEA finalisation ####

  cat(paste0("Creating Enrichment Analysis Object\n"))
  X$GSEA$GSEA_Results_Axis <- axis_gsea
  X$GSEA$RankingAxis <- axis_rank
  X$GSEA$GSEA_Results <- cluster_gsea
  X$GSEA$Ranking <- cluster_rank
  X$GSEA$GMTfile <- GMTfile
  X$GSEA$Pathways <- axis_gsea$Axis1$pathway
  X$GSEA$AllRanking <- axis_rank %>% append(cluster_rank)
  X$GSEA$gseaParam <- gseaParam
  X$Shiny <- create_dashboard3(X)
  cat(paste0("Enrichment Analysis Completed\n"))
  class(X$GSEA) <- "GSEA"
  return(X)
}



#' GSEA at single cell level
#'
#' Compute the enrichment score of every cell to a given pathway set using the gene ranking derived from MCA.
#'
#' @param X A MCXPress object containing a Dimred object
#' @param GMTfile named list of pathway with their representing gene. 
#' @param nperm the number of permutation.
#' @param minSize threshold for the minimum number of genes that must be in the pathway.
#' @param maxSize threshold for the maximum number of genes that must be in the pathway.
#' @param nproc number of processor to use for the enrichment algorithm.
#' @param nbin an integer indicating which bin to use for the enrichment analysis.
#' @param naxis number of axis to perform the enrichment analysis.
#' @param gseaParam an integer to weight the enrichment analysis.
#' @param dim Number of Axis to retain for the calculation gene ranking for each cells
#' @return Return a MCXpress object containing a GSEA object
#' \item{GSEA_Results}{List containing dataframes of the GSEA analysis for each Cluster}
#' \item{Ranking}{List containing vectors of gene ranking according to their distance in MCA for each Cluster}
#' \item{GSEA_Results_Axis}{List containing dataframes of the GSEA analysis for each Axis}
#' \item{Ranking_Axis}{List containing vectors of gene ranking according to their MCA cells coordinate and genes expression correlation for each Axis}
#' \item{GMTfile}{GMTfile used for GSEA}
#' \item{Pathway}{Name of the genesets contained in the gmtfile}
#' \item{gseaParam}{GSEA Parameter used}
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553)
#' MCX64553 <- SC_GSEAparall(X = MCX64553, GMTfile = reactome_gmtfile, nperm = 10000, dim=3)
#' @export
SC_GSEAparall <- function(X, GMTfile, nperm = 1000, minSize = 15, maxSize = 500,
                       nproc = 4, nbin = 1, naxis = 2, gseaParam = 0, dim=2)
{
  ## ............................................................................
  ## A Axis Correlation Gene Set Enrichment Analysis ####

  ## ............................................................................
  ## B Gene Set Enrichment Analysis on Cluster ####

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  a Filter Bin ####
  cells_coord <-  X$MCA$cells_principal[,1:dim]
  genes_coord <-  X$MCA$genes_standard[,1:dim]
  df2 <-cell_gene_distances <- fields::rdist(x1 = cells_coord %>% select(contains("Axis")) %>%  as.matrix %>%  set_rownames(cells_coord %>%  rownames),
                                             x2 = genes_coord %>%  as.matrix) %>% t() %>% set_colnames(cells_coord %>%  rownames) %>%
    set_rownames(genes_coord %>% rownames) %>% as.data.frame() %>% rownames_to_column(var = "Genes") %>%
    as_tibble() %>% tidyr::separate(col = Genes, into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
    dplyr::filter(bin %in% nbin)

  if (nbin %>%  length %>%  equals(1) %>%  not){
    df2 <- df2 %>% group_by(Genes) %>% summarise_at(-1,.funs = min) %>% ungroup %>%  mutate(bin="mix")
  }

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  b Calculate Rank for Cluster ####
  Scaling2 <- function(x){2*(x-min(x))/(max(x)-min(x))}
  cat("\nCalculating ranking of genes for each clusters \n")
  # cluster_rank <- df2 %>% dplyr::select(-Genes,-bin) %>% bplapply(FUN = function(x,y){
  #   1-(sort(Scaling2(set_names(x,y))))
  # },BPPARAM = SnowParam(workers = nproc, tasks=nproc, progressbar = T), y=df2$Genes)
  cat("\nCalculating ranking of genes for each clusters \n")
  cluster_rank <- df2 %>% dplyr::select(-Genes,-bin) %>% lapply(FUN = function(x)
  {
    1-(x %>% set_names(df2$Genes) %>% (function(x){2*(x-min(x))/(max(x)-min(x))}) %>% sort)
  })
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ### . . . . . . . ..  c fgsea analysis for cluster ####

  cat("Beginning enrichment analysis for cells\n\n")
  cluster_gsea <-
    cluster_rank %>% BiocParallel::bplapply(
      FUN = parallel_fgsea,
      BPPARAM = BiocParallel::SnowParam(workers = nproc, tasks=nproc, progressbar = T),
      a = GMTfile,
      b = nperm,
      c = minSize,
      d = maxSize,
      e = nproc,
      f = gseaParam
    )


  ## ............................................................................
  ## C GSEA finalisation ####

  cat(paste0("Creating SIigle Cell Enrichment Analysis Object\n"))
  X$SC_GSEA$GSEA_Results <- cluster_gsea
  X$SC_GSEA$Ranking <- cluster_rank
  X$SC_GSEA$Pathways <- X$SC_GSEA$GSEA_Results %>%  extract2(1) %>%  use_series(pathway)
  X$SC_GSEA$gseaParam <- gseaParam
  X$SC_GSEA$GMTfile <- GMTfile
  cat(paste0("Single Cell Enrichment Analysis Completed\n"))
  class(X$SC_GSEA) <- "GSEA"
  return(X)
}







#' Interactive plot enrichment score of a pathway
#'
#' Plotly version of the plotEnrichment function of the fgsea package
#'
#' @param pathway A single geneset
#' @param stats Ranking of the genes
#' @param gseaParam GSEA parameter value
#' @return enrichment plot with plotly
plotlyEnrichment <- function(pathway, stats, gseaParam = 0)
{
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    name <- rnk[(rnk %>% names) %in% pathway] %>% sort %>% names
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
    plot_ly(data = toPlot) %>% add_lines(x = ~x, y = ~y, mode = "lines",
        line = list(color = "rgb(150, 240, 30)", width = 2),
        name = "Enrichment") %>% add_lines(x = ~x, y = ~min(bottoms),
        line = list(color = "rgb(110, 190, 250)", width = 2,
            dash = "dash"), name = "Lower limit", hoverinfo = "text",
        text = ~round(min(bottoms), digits = 4)) %>% add_lines(x = ~x,
        y = ~max(tops), line = list(color = "rgb(250, 150, 10)",
            width = 2, dash = "dash"), name = "Upper limit",
        hoverinfo = "text", text = ~round(max(tops), digits = 4)) %>%
        add_segments(x = ~pathway, xend = ~pathway, y = ~diff/2,
            yend = ~-diff/2, line = list(color = "rgb(0, 10, 10)",
                width = 1), text = ~paste0(name), showlegend = FALSE, hoverinfo = "text")
}

