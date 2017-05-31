#' @import magrittr
#' @import plotly
#' @import ggplot2
#' @import shiny
#' @import shinydashboard
#' @import shinythemes
#' @import BiocParallel
#' @import fgsea
#' @import statmod
#' @import diptest
#' @importFrom splines bs
#' @importFrom cluster pam
#' @importFrom  fields rdist
#' @import tibble
#' @import dplyr
#' @importFrom tidyr spread gather separate gather_ spread_ separate_
#' @importFrom purrr map map_df
#' @importFrom stringr str_length
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @importFrom stats cor kmeans na.omit quantile setNames var wilcox.test as.dist cutree hclust
#' @import lazyeval
NULL

#'Expression Matrix from GEO ID: GSE64553
#'
#'The dataset contains only the value for HFF samples
#'and one outlier has been removed (HFF_PD58_woRotenone2).
#'
#'@format A matrix with 47966 genes (row) and 35 individuals (col)
#'@source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64553}
"GSE64553"

#'gmtfile from reactome pathway
#'
#'Contains all the geneset of the reactome pathway for human.
#'@format A list of 1892 reactome pathway geneset with hg19 symbol.
#'@source \url{http://www.reactome.org/pages/download-data/}
"reactome_gmtfile"

