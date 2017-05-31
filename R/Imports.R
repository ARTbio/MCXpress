#' @import magrittr
#' @import plotly
#' @import shiny
#' @import fgsea
#' @import shinythemes
#' @import shinydashboard
#' @import ggplot2 except=c(last_plot)
#' @import diptest
#' @import dplyr
#' @import tibble
#' @importFrom tidyr spread gather separate gather_ spread_ separate_
#' @import BiocParallel
#' @importFrom purrr map dmap
#' @importFrom stringr str_length
#' @import statmod
#' @importFrom cluster pam
#' @importFrom splines bs
#' @importFrom  fields rdist
#' @importFrom stats cor kmeans na.omit quantile setNames var wilcox.test as.dist cutree hclust
#' @import lazyeval
#' @importFrom utils head setTxtProgressBar txtProgressBar
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

