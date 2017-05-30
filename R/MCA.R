
#   ____________________________________________________________________________
#   MCA Dimension Reduction                                                 ####
#' Dimensional reduction using multiple corespondence analysis.
#'
#' Apply fuzzy multiple correspondance analysis on the disjunctive matrix to reduce the number of dimension
#' and keep the axis with the most informative variance
#'
#' @param X A MCXpress object containing a Disjunctive Matrix produced by one of the Discretisation method.
#' @param Dim Number of Dimension to keep.
#' @return A MCXpress object containing a MCXmca object.
#' \itemize{
#'   \item{"parameter 1"}{Stuff}
#'   \item{"parameter 2"}{Stuff}
#'  }
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- Discretisation_Range_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553, Dim = 5)
#' @export
MCA <- function(X, Dim = (X$ExpressionMatrix %>%  ncol)-1){
    cat("Peforming MCA...\n")
##  ............................................................................
##  A MCA Algorithm                                                           ####
    pb<-txtProgressBar()
    Acol <- colSums(X$Disjunctive_Matrix)
    Y <- sweep(X$Disjunctive_Matrix,2, sqrt(Acol) ,"/")
    setTxtProgressBar(pb, 0.1)
    Arow <- rowSums(X$Disjunctive_Matrix)
    Z <- Y/sqrt(Arow)

    if((X$Disjunctive_Matrix %>%  ncol) %>% is_weakly_greater_than(X$Disjunctive_Matrix %>%  nrow)){
      setTxtProgressBar(pb, 0.2)
      Ztcross <- tcrossprod(Z)
      setTxtProgressBar(pb, 0.3)
      Eigen <- eigen(Ztcross)
      Eig <- Eigen$values[-1]
      Eig<-subset(Eig,Eig >0)
      V <- Eigen$vectors[,2:(Eig %>% length %>% add(1))]
      setTxtProgressBar(pb, 0.4)
    } else{
      setTxtProgressBar(pb, 0.2)
      Zcross <- crossprod(Z)
      setTxtProgressBar(pb, 0.3)
      Eigen  <- eigen(Zcross)
      Eig <- Eigen$values[-1]
      Eig<-subset(Eig,Eig >0)
      U   <- Eigen$vectors[,2:(Eig %>% length %>% add(1))]
      if(Dim>(Eig %>% length)){Dim <- Eig %>% length}
      Dt  <- (Eig %>%  sqrt %>% raise_to_power(-1) %>%   diag)
      V   <- (Z %*% U) %*% Dt
      setTxtProgressBar(pb, 0.4)
    }
##  ............................................................................
##  B Coordinates Calculation                                               ####

    D <- Eig %>%  diag
    Dc <- (1/(sqrt(Acol/sum(X$Disjunctive_Matrix))))
    Component <- paste0("Axis", 1:Dim)
    setTxtProgressBar(pb, 0.6)

    X$MCA$Cells_Principal <- (sqrt(nrow(X$Disjunctive_Matrix)) *
        V)[, 1:Dim] %>% set_colnames(Component) %>%
        set_rownames(X$Disjunctive_Matrix %>% rownames) %>% data.frame

    setTxtProgressBar(pb, 0.7)

    X$MCA$Cells_Standard <- (sqrt(nrow(X$Disjunctive_Matrix)) *
        (V %*% sqrt(D)))[, 1:Dim] %>% set_colnames(Component) %>%
        set_rownames(X$Disjunctive_Matrix %>% rownames) %>% data.frame

    setTxtProgressBar(pb, 0.8)
    X$MCA$Genes_Standard <- ((t(Z) %*% V) * Dc)[,
        1:Dim] %>% set_colnames(Component) %>% set_rownames(X$Disjunctive_Matrix %>%
        colnames) %>% data.frame
    setTxtProgressBar(pb, 0.9)
    X$MCA$Genes_Principal <- sweep(X$MCA$Genes_Standard,
        2, sqrt(Eig)[1:(Dim)], "/") %>%  set_colnames(Component) %>% set_rownames(X$Disjunctive_Matrix %>%
        colnames) %>% data.frame
    setTxtProgressBar(pb, 1)
    close(pb)
##  ............................................................................
##  C Distance Calculation                                                  ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Cell-Cell Distance                                                    ####
    cat("Calculating Cell to Cell Distance...\n")
    X$MCA$Cell2Cell_Distance <- X$MCA$Cells_Standard[,
        1:Dim] %>% rdist %>% set_colnames(X$MCA$Cells_Principal %>%
        rownames) %>% set_rownames(X$MCA$Cells_Principal %>%
        rownames)
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### b Cell-Gene Distance                                                    ####
    # X$MCA$Cell2Gene_Distance <- rdist(X$MCA$Cells_Principal,
    #     X$MCA$Genes_Standard) %>% set_colnames(X$MCA$Genes_Standard %>%
    #     rownames) %>% set_rownames(X$MCA$Cells_Principal %>%
    #     rownames)

    X$MCA$Eigen_Value <- Eig[1:Dim] %>% set_names(Component)
    Percentage_Variance <- (Eig %>% prop.table)*100
    Percentage_Variance_Cum <- Percentage_Variance %>%  cumsum()
    Axes <- paste0("Axis",1:length(Eig))
    X$MCA$Explained_Eigen_Variance <- tibble(factor(Axes, levels = Axes), Percentage_Variance,Percentage_Variance_Cum) %>% set_colnames(c("Axis","Explained_Variance","Cumulative"))
    X$MCA$Methods <- "MCA"


    # End Calculate Wilcoxon for Cell2CellDistance

    # Calculate Correlation Axis and Genes

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Axis and Gene Correlation                                               ####
    X$MCA$Axis_Gene_Cor <- cor(X$Disjunctive_Matrix,
        X$MCA$Cells_Principal, method="spearman") %>% data.frame() %>%
        rownames_to_column(var = "Genes") %>% as_tibble()
    X$MCA$Axis_Gene_Cor[,-1] <- X$MCA$Axis_Gene_Cor[,-1] %>%  dmap(.f = round , digits=3)
    # End Calculate Correlation Axis and Genes
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Cell Principal Space                                                    ####

    X$MCA$Plot <- X$MCA$Cells_Principal %>%
        ggplot2::ggplot(ggplot2::aes(x = Axis1, y = Axis2)) + ggplot2::geom_point() +
        ggplot2::theme_light()
    X$Shiny <- X %>% Create_Dashboard1()
    class(X$MCA) <- "MCA_Object"
    class(X) <- "MCXpress_object"
    cat("MCA is finished \n")
    return(X)
}
