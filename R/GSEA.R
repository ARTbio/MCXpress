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
#   ____________________________________________________________________________
#   Gene set enrichment analysis                                            ####
GSEA <- function(X, GMTfile, nperm = 1000,
  minSize = 15, maxSize = 500, nproc = 1, nbin = 1, naxis = 2,
  gseaParam = 0) {

##  ............................................................................
##  A Axis Correlation Gene Set Enrichment Analysis                         ####


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Filtering of bin                                                      ####
  df <- X$MCA$Axis_Gene_Cor[, 1:(naxis + 1)] %>% separate(col = Genes,
    into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
    filter(bin == 1) %>% select(-bin)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### b Correlation Ranking Calculation                                       ####
  cat("Calculating ranking of genes correlation for each axis \n")
  axis_rank <- df %>% select(-Genes) %>% purrr::map(function(x) {
    x %>% abs %>% set_names(df$Genes) %>% rank(ties.method = "first") %>%
      sort
  })

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### c Fast Gene Set Enrichment Analysis on each Axis                        ####

  cat("Beginning enrichment analysis for axis \n")
  axis_gsea <- axis_rank %>% purrr::map2(.y = axis_rank %>%
    names, .f = function(x, y) {
    cat(paste("processing:", y, "\n"))
    gsea <- fgsea(pathways = GMTfile, stats = x, nperm = nperm,
      maxSize = maxSize, minSize = minSize, nproc = nproc,
      BPPARAM = SerialParam(), gseaParam = gseaParam) %>%
      as_tibble
    ins <- gsea %>% magrittr::extract(, 2:5) %>% round(digits = 5)
    val <- gsea %>% magrittr::inset(, 2:5, value = ins)
    return(val)
  })

##  ............................................................................
##  B Gene Set Enrichment Analysis on Cluster                               ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Filter Bin                                                            ####
 df2 <- X$cluster$gene_cluster_distances %>% tidyr::separate(col = Genes,
    into = c("Genes", "bin"), sep = "-bin", convert = TRUE) %>%
    dplyr::filter(bin %in% nbin) %>% dplyr::group_by(Genes) %>% dplyr::summarise_at(.cols = -(1:2),
    .funs = min)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### b Calculate Rank for Cluster                                              ####

  cat("\nCalculating ranking of genes for each clusters \n")
  cluster_rank <- df2 %>% select(-Genes) %>% purrr::map(function(x) {
    x %>% abs %>% multiply_by(-1) %>% set_names(df2$Genes) %>%
      rank(ties.method = "first") %>% sort
  })

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### c fgsea analysis for cluster                                            ####

  cat("Beginning enrichment analysis for clusters\n")
  cluster_gsea <- cluster_rank %>% purrr::map2(.y = cluster_rank %>%
    names, .f = function(x, y) {
    cat(paste("processing:", y, "\n"))
    gsea <- fgsea(pathways = GMTfile, stats = x, nperm = nperm,
      maxSize = maxSize, minSize = minSize, nproc = nproc,
      BPPARAM = SerialParam(), gseaParam = gseaParam) %>%
      as_tibble()
    ins <- gsea %>% magrittr::extract(, 2:5) %>% round(digits = 5)
    val <- gsea %>% magrittr::inset(, 2:5, value = ins)
    return(val)
  })

##  ............................................................................
##  C GSEA finalisation                                                       ####

  X$GSEA$GSEA_Results_Axis <- axis_gsea
  X$GSEA$RankingAxis <- axis_rank
  X$GSEA$GSEA_Results <- cluster_gsea
  X$GSEA$Ranking <- cluster_rank
  X$GSEA$GMTfile <- GMTfile
  X$GSEA$Pathways <- axis_gsea$Axis1$pathway
  X$GSEA$AllRanking <- X$GSEA$RankingAxis %>%
    append(X$GSEA$Ranking)
  X$GSEA$gseaParam <- gseaParam
  X$Shiny <- create_dashboard3(X)
  cat(paste0("Enrichment Analysis Completed\n"))
  class(X$GSEA) <- "GSEA_Object"
  return(X)
}

#' Interactive plot enrichment score of a pathway
#'
#' Plotly version of the plotEnrichment function of the fgsea package
#'
#' @param pathway A single geneset
#' @param stats Ranking of the genes
#' @param gseaParam GSEA parameter value
#'
#' @return RETURN_DESCRIPTION
plotlyEnrichment <- function(pathway, stats, gseaParam = 0) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  name <- rnk[(rnk %>% names) %in% pathway] %>% sort %>%  names
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
    line = list(color = "rgb(154, 240, 24)", width = 2),
    name = "Enrichment") %>% add_lines(x = ~x, y = ~min(bottoms),
    line = list(color = "rgb(110, 193, 248)", width = 2,
      dash = "dash"), name = "Lower limit", hoverinfo = "text",
    text = ~round(min(bottoms), digits = 4)) %>% add_lines(x = ~x,
    y = ~max(tops), line = list(color = "rgb(250, 150, 10)",
      width = 2, dash = "dash"), name = "Upper limit",
    hoverinfo = "text", text = ~round(max(tops), digits = 4)) %>% add_segments(x =~pathway, xend=~pathway ,y =~diff/2 ,yend =~-diff/2, line=list(color = "rgb(0, 10, 10)",width=1), text=~paste0(name,
      "</br>", "Rank:", pathway), showlegend = FALSE, hoverinfo="text")
}

