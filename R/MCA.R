##MCA
# wants <- c("survival")
# has   <- wants %in% rownames(installed.packages())
# if(any(!has)) install.packages(wants[!has])
####Library
library(tidyverse)
library(statmod)
library(diptest)
library(splines)
library(abind)
library(Rtsne)
library(plotly)
library(shiny)
library(fields)
library(igraph)
library(fgsea)
library(BiocParallel)
library(stringr)
library(magrittr)
library(ade4)
library(reshape2)
library(DT)
library(fpc)
library(shinythemes)

select <- dplyr::select
transmute <- dplyr::transmute

#' Performs Multiple Corespondence Analysis on a object of class MCXpress containing a Disjunctive Matrix.
#'
#' @param X A MCXpress class object containing a Disjunctive Matrix.
#' @param Dim number of Axis to be kept for the variables and individual coordinates and used for the distance calculation.
#' @return
#' @seealso
Dimension_reduction_MCA <- function(X, Dim = 5) {
  cat('Beginning MCA...\n')
  dis_matrix<- X$Disjunctive_Matrix
  MCA <- dis_matrix %>% as.data.frame(row.names = dis_matrix %>% rownames)
  MCA <- prep.fuzzy.var(MCA, rep(2, (MCA %>% ncol) / 2))
  MCA_results <-dudi.fca(MCA, nf = Dim, scannf = FALSE)
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
  X$Dim_Red$Cell2Cell_Distance <-X$Dim_Red$Cells_Standard[,1:Dim] %>% rdist
  X$Dim_Red$Cell2Cell_Distance %<>%  set_colnames(X$Dim_Red$Cells_Principal %>%  rownames)
  X$Dim_Red$Cell2Cell_Distance %<>%  set_rownames(X$Dim_Red$Cells_Principal %>%  rownames)
  X$Dim_Red$Cell2Gene_Distance <- rdist(X$Dim_Red$Cells_Principal, X$Dim_Red$Genes_Standard) %>%  set_colnames(X$Dim_Red$Genes_Standard %>% rownames) %>%  set_rownames(X$Dim_Red$Cells_Principal %>% rownames)
  #End Distance Calculation
  X$Dim_Red$Eigen_Value     <-
    MCA_results$eig %>% set_names(Component)
  EigCum <-
    ((MCA_results$eig * 100) / (MCA_results$eig %>% sum)) %>%  as.data.frame() %>%  mutate(PC =
                                                                                             (paste0("PC", c(
                                                                                               1:length(MCA_results$eig)
                                                                                             )))) %>% mutate(Cumul = ((MCA_results$eig * 100) / (MCA_results$eig %>% sum)) %>% cumsum)
  colnames(EigCum)[1] <- "EigenValue"
  EigCum <- EigCum[c(2, 1, 3)]
  EigCum <-EigCum[1:(X$Dim_Red$Cells_Standard %>%  nrow),] %>% na.omit
  X$Dim_Red$Cumul <- EigCum %>%  gather(EigenValue, Cumul, key = "Type", value = "Value")
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


