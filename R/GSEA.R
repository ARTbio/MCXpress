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
#' @param bin an integer indicating which bin to use for the enrichment analysis.
#' @param naxis number of axis to perform the enrichment analysis.
#' @param gseaParam an integer to weight the enrichment analysis.
#' @return Return a MCXpress object containing a
#' @examples
#   ____________________________________________________________________________
#   Gene set enrichment analysis                                            ####
GSEA <- function(X, GMTfile, nperm = 1000,
  minSize = 15, maxSize = 500, nproc = 1, nbin = 1, naxis = 2,
  gseaParam = 0) {

##  ............................................................................
##  A Axis Correlation Gene Set Enrichment Analysis                         ####


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Filtering of bin                                                      ####
  df <- X$Dim_Red$Axis_Gene_Cor[, 1:(naxis + 1)] %>% separate(col = Genes,
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
 df2 <- X$cluster$Gene_Cluster_Distance %>% tidyr::separate(col = Genes,
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

  X$Functionnal_Analysis$GSEA_Results_Axis <- axis_gsea
  X$Functionnal_Analysis$RankingAxis <- axis_rank
  X$Functionnal_Analysis$GSEA_Results <- cluster_gsea
  X$Functionnal_Analysis$Ranking <- cluster_rank
  X$Functionnal_Analysis$GMTfile <- GMTfile
  X$Functionnal_Analysis$Pathways <- axis_gsea$Axis1$pathway
  X$Functionnal_Analysis$AllRanking <- X$Functionnal_Analysis$RankingAxis %>%
    append(X$Functionnal_Analysis$Ranking)
  X$Functionnal_Analysis$gseaParam <- gseaParam
  X$Functionnal_Analysis$Shiny <- Create_Shiny_Functionnal_Analysis(X)
  cat(paste0("Enrichment Analysis Completed\n"))
  class(X$Functionnal_Analysis) <- "FA.object"
  return(X)
}

