
#   ____________________________________________________________________________
#   MCA Dimension Reduction                                                 ####
#' Dimensional reduction using multiple corespondence analysis
#'
#' @param X A MCXpress object containing a Disjunctive Matrix produced by one of the Discretisation method.
#' @param Dim Number of Dimension to use for the cell to cell distance calculation.
#'
#' @return A MCXpress objetc containing a MCXmca object.
#' @export
#'
#' @examples
reduce_dimension_mca <- function(X, Dim = (X$ExpressionMatrix %>%  ncol)-1){
    cat("Beginning MCA...\n")
##  ............................................................................
##  A MCA Algorithm                                                           ####
    Acol <- colSums(X$Disjunctive_Matrix)
    Arow <- rowSums(X$Disjunctive_Matrix)
    Y    <- apply(X$Disjunctive_Matrix, 1, FUN = function(x){
            x/sqrt(Acol)
        })
        Z <- apply(Y, 1, FUN = function(x){
            x/sqrt(Arow)
        })
    Eigen <- eigen(tcrossprod(Z))

##  ............................................................................
##  B Coordinates Calculation                                               ####
    V <- Eigen$vectors[,-1]
    Eig <- Eigen$values[-1]
    D <- Eig %>%  diag
    Dc <- (1/(sqrt(Acol/sum(X$Disjunctive_Matrix))))
    Component <- paste0("Axis", 1:Dim)
    X$Dim_Red$Cells_Principal <- (sqrt(nrow(X$Disjunctive_Matrix)) *
        V)[, 1:Dim] %>% set_colnames(Component) %>%
        set_rownames(X$Disjunctive_Matrix %>% rownames) %>% data.frame
    X$Dim_Red$Cells_Standard <- (sqrt(nrow(X$Disjunctive_Matrix)) *
        (V %*% sqrt(D)))[, 1:Dim] %>% set_colnames(Component) %>%
        set_rownames(X$Disjunctive_Matrix %>% rownames) %>% data.frame
    X$Dim_Red$Genes_Standard <- ((t(Z) %*% V) * Dc)[,
        1:Dim] %>% set_colnames(Component) %>% set_rownames(X$Disjunctive_Matrix %>%
        colnames) %>% data.frame
    X$Dim_Red$Genes_Principal <- sweep(X$Dim_Red$Genes_Standard,
        2, sqrt(Eig)[2:(Dim + 1)], "/") %>%  set_colnames(Component) %>% set_rownames(X$Disjunctive_Matrix %>%
        colnames) %>% data.frame

##  ............................................................................
##  C Distance Calculation                                                  ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### a Cell-Cell Distance                                                    ####
    cat("Calculating Cell to Cell Distance...\n")
    X$Dim_Red$Cell2Cell_Distance <- X$Dim_Red$Cells_Standard[,
        1:Dim] %>% rdist %>% set_colnames(X$Dim_Red$Cells_Principal %>%
        rownames) %>% set_rownames(X$Dim_Red$Cells_Principal %>%
        rownames)
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### b Cell-Gene Distance                                                    ####
    X$Dim_Red$Cell2Gene_Distance <- rdist(X$Dim_Red$Cells_Principal,
        X$Dim_Red$Genes_Standard) %>% set_colnames(X$Dim_Red$Genes_Standard %>%
        rownames) %>% set_rownames(X$Dim_Red$Cells_Principal %>%
        rownames)

    X$Dim_Red$Eigen_Value <- Eig[1:Dim] %>% set_names(Component)
    EigCum <- ((Eig * 100)/(Eig %>% sum)) %>% as.data.frame() %>%
        mutate(Axis = (paste0("Axis", c(1:length(Eig))))) %>%
        mutate(Cumul = ((Eig * 100)/(Eig %>% sum)) %>%
            cumsum)
    colnames(EigCum)[1] <- "EigenValue"
    EigCum <- EigCum[c(2, 1, 3)]
    EigCum <- EigCum[1:(X$Dim_Red$Cells_Standard %>%
        nrow), ] %>% na.omit
    X$Dim_Red$Cumul <- EigCum %>% gather(EigenValue,
        Cumul, key = "Type", value = "Value")
    X$Dim_Red$Methods <- "MCA"

    # Calculate Wilcoxon for Cell2CellDistance
    # cat("Calculating Distance Wilcoxon...\n")
    # Wilcox <- vector(length = (X$Dim_Red$Cells_Standard %>%
    #     ncol) - 1)
    # for (i in 1:((X$Dim_Red$Cells_Standard %>% ncol) -
    #     1)) {
    #     dis1 <- X$Dim_Red$Cells_Standard[, 1:i] %>%
    #         rdist %>% as.matrix %>% as.vector
    #     dis2 <- X$Dim_Red$Cells_Standard[, 1:(i + 1)] %>%
    #         rdist %>% as.matrix %>% as.vector
    #     Wilcox[i] <- wilcox.test(dis1, dis2, alternative = "less",
    #         paired = TRUE) %>% extract2("p.value")
    #     names(Wilcox)[i] <- paste0("Axis", i, "-Axis",
    #         i + 1)
    # }
    # X$Dim_Red$Wilcoxon <- tibble(Wilcox %>% names,
    #     Wilcox) %>% set_colnames(c("Axis", "p_value"))
    # End Calculate Wilcoxon for Cell2CellDistance

    # Calculate Correlation Axis and Genes

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Axis and Gene Correlation                                               ####

    X$Dim_Red$Axis_Gene_Cor <- cor(X$Disjunctive_Matrix,
        X$Dim_Red$Cells_Principal) %>% data.frame() %>%
        rownames_to_column(var = "Genes") %>% as_tibble()
    X$Dim_Red$Axis_Gene_Cor[,-1] <- X$Dim_Red$Axis_Gene_Cor[,-1] %>%  dmap(.f = round , digits=3)
    # End Calculate Correlation Axis and Genes


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Cell Principal Space                                                    ####

    X$Dim_Red$Plot <- X$Dim_Red$Cells_Principal %>%
        ggplot(aes(x = Axis1, y = Axis2)) + geom_point() +
        theme_light()
    X$Dim_Red$Shiny <- X %>% Create_Shiny_Dim_Red
    class(X$Dim_Red) <- "Dim_Red_Object"
    class(X) <- "MCXpress_object"
    cat("MCA is finished \n")
    return(X)
}
