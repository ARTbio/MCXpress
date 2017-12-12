# ____________________________________________________________________________
# MCA Dimension Reduction ####
#' Dimensional reduction using multiple corespondence analysis.
#'
#' Apply fuzzy multiple correspondance analysis on the disjunctive matrix to reduce the number of dimension
#' and keep the axis with the most informative variance
#'
#' @param X A MCXpress object containing a Disjunctive Matrix produced by one of the Discretisation method.
#' @param Dim Number of Dimension to keep.
#' @return A MCXpress object containing a MCXmca object.
#' \itemize{
#'   \item{'cells_principal'}{Stuff}
#'   \item{'parameter 2'}{Stuff}
#'  }
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' MCX64553 <- discretisation_01(MCX64553, scaled=FALSE)
#' MCX64553 <- MCA(MCX64553, Dim = 5)
#' @export
MCA <- function(X, Dim = (X$ExpressionMatrix %>% dim() %>% min()) - 1) {
  cat("Peforming MCA...\n")
  ## ............................................................................  A
  ## MCA Algorithm ####
  pb <- txtProgressBar(width = 50, style = 3, char = "+")
  Acol <- colSums(X$disjunctive_matrix)
  Y <- sweep(X$disjunctive_matrix, 2, sqrt(Acol), "/")
  setTxtProgressBar(pb, 0.1)
  Arow <- rowSums(X$disjunctive_matrix)
  Z <- Y / sqrt(Arow)

  if ((X$disjunctive_matrix %>% ncol()) %>% is_weakly_greater_than(X$disjunctive_matrix %>%
    nrow())) {
    setTxtProgressBar(pb, 0.2)
    Ztcross <- tcrossprod(Z)
    setTxtProgressBar(pb, 0.3)
    Eigen <- eigen(Ztcross)
    Eig <- Eigen$values[-1]
    Eig <- subset(Eig, Eig > 0)
    V <- Eigen$vectors[, 2:(Eig %>% length() %>% add(1))]
    setTxtProgressBar(pb, 0.4)
  } else {
    setTxtProgressBar(pb, 0.2)
    Zcross <- crossprod(Z)
    setTxtProgressBar(pb, 0.3)
    Eigen <- eigen(Zcross)
    Eig <- Eigen$values[-1]
    Eig <- subset(Eig, Eig > 0)
    U <- Eigen$vectors[, 2:(Eig %>% length() %>% add(1))]
    if (Dim > (Eig %>% length())) {
      Dim <- Eig %>% length()
    }
    Dt <- (Eig %>% sqrt() %>% raise_to_power(-1) %>% diag())
    V <- (Z %*% U) %*% Dt
    setTxtProgressBar(pb, 0.4)
  }
  ## ............................................................................  B
  ## Coordinates Calculation ####

  D <- Eig %>% diag()
  Dc <- (1 / (sqrt(Acol / sum(X$disjunctive_matrix))))
  Component <- paste0("Axis", 1:Dim)
  setTxtProgressBar(pb, 0.6)

  X$MCA$cells_principal <- (sqrt(nrow(X$disjunctive_matrix)) * V)[, 1:Dim] %>%
    set_colnames(Component) %>%
    set_rownames(X$disjunctive_matrix %>% rownames()) %>%
    data.frame()

  setTxtProgressBar(pb, 0.7)

  X$MCA$cells_standard <- (sqrt(nrow(X$disjunctive_matrix)) * (V %*% sqrt(D)))[
    ,
    1:Dim
  ] %>% set_colnames(Component) %>% set_rownames(X$disjunctive_matrix %>%
    rownames()) %>% data.frame()

  setTxtProgressBar(pb, 0.8)
  X$MCA$genes_standard <- ((t(Z) %*% V) * Dc)[, 1:Dim] %>%
    set_colnames(Component) %>%
    set_rownames(X$disjunctive_matrix %>% colnames()) %>%
    data.frame()
  setTxtProgressBar(pb, 0.9)
  X$MCA$genes_principal <- sweep(X$MCA$genes_standard, 2, sqrt(Eig)[1:(Dim)], "/") %>%
    set_colnames(Component) %>%
    set_rownames(X$disjunctive_matrix %>% colnames()) %>%
    data.frame()
  setTxtProgressBar(pb, 1)
  close(pb)
  ## ............................................................................  C
  cat("Eigen Value Table Generation \n \n")
  X$MCA$eigen_value <- Eig[1:Dim] %>% set_names(Component)
  Percentage_Variance <- (Eig %>% prop.table()) * 100
  Percentage_Variance_Cum <- Percentage_Variance %>% cumsum()
  Axes <- paste0("Axis", 1:length(Eig))
  X$MCA$explained_eigen_variance <- tibble::tibble(
    factor(Axes, levels = Axes), Percentage_Variance,
    Percentage_Variance_Cum
  ) %>% set_colnames(c(
    "Axis", "Explained_Variance",
    "Cumulative"
  ))
  X$MCA$Methods <- "MCA"

  # Calculate Correlation Axis and Genes
  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  . . . . . . . ..
  ### Axis and Gene Correlation ####
  cat("Calculating Spearman Correlation s\n \n")
  X$MCA$Axis_Gene_Cor <- cor(X$disjunctive_matrix, X$MCA$cells_principal[, 1:min(X$MCA$cells_principal %>% ncol(), 5)], method = "spearman") %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "Genes") %>%
    tibble::as_tibble()
  X$MCA$Axis_Gene_Cor[, -1] <- X$MCA$Axis_Gene_Cor[, -1] %>% sapply(
    FUN = round,
    digits = 3
  )
  # End Calculate Correlation Axis and Genes . . . . . . . . .  . . . . . . . . . .
  # . . . . . . . . . . . . . . . . . . ..  Cell Principal Space ####

  X$MCA$plot <- X$MCA$cells_principal %>% ggplot2::ggplot(ggplot2::aes(
    x = Axis1,
    y = Axis2
  )) + ggplot2::geom_point() + ggplot2::theme_light()
  X$Shiny <- X %>% create_dashboard1()
  class(X$MCA) <- "MCA"
  class(X) <- "MCXpress"
  cat("MCA is finished \n")
  return(X)
}
