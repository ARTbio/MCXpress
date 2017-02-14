#' @import magrittr
#' @import tidyverse
#' @import plotly except=last_plot
#' @import shiny
#' @import fgsea
#' @import shinythemes
#' @import ggplot2
#' @import diptest
#' @import dplyr
#' @importFrom tidyr spread gather separate
#' @import tibble
#' @importFrom DT datatable
#' @import BiocParallel
#' @import stringr
#' @import statmod
#' @import splines
#' @importFrom  fields rdist
#' @importFrom igraph graph.adjacency delete_edges clusters delete_vertices
#' @importFrom stats cor kmeans na.omit quantile setNames var wilcox.test
#' @import lazyeval
#' @importFrom ade4 as.dudi prep.fuzzy.var dudi.fca
#' @importFrom utils head
NULL

#'Expression Matrix from GEO ID: GSE64553
#'
#'The dataset contains only the value for HFF samples
#'and one outlier has been removed.
#'
#'@format A matrix with 47966 genes (row) and 35 individuals (col)
#'@source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64553}
"GSE64553"
