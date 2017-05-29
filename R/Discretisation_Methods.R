Discretisation_Bsplines <- function(X) {
  nbins <-  2
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<-X$ExpressionMatrix
    nc = ncol(exp_matrix)
    ng = nrow(exp_matrix)
    discreteMatrix = array(0, dim = c(ng, nc, nbins))
    for (counterGene in seq(ng)) {
      j = bs(
        scale(exp_matrix[counterGene, ]),
        df = nbins,
        degree = 2,
        intercept = T
      )
      for (counterCondition in seq(nc)) {
        if (nbins == 2) {
          newJ <- cbind(rowSums(j[, 1:2]), j[, 3:dim(j)[2]])
          discreteMatrix[counterGene, counterCondition, ] = newJ[counterCondition, ]
        }
        else {
          discreteMatrix[counterGene, counterCondition, ] = j[counterCondition, ]
        }
      }
    }
    rownames(discreteMatrix) <- rownames(exp_matrix)
    colnames(discreteMatrix) <- colnames(exp_matrix)
    discreteMatrix<-discreteMatrix[,,rev(1:dim(discreteMatrix)[3])]
    # Note: this is a fuzzy coded expression matrix
    Disjunctive_Matrix <- discreteMatrix %>% Create_Disjunctive_Matrix
    X$Disjunctive_Matrix<-Disjunctive_Matrix
    return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)}
}

Discretisation_Bsplines2 <- function(X) {
  nbins <-  4
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<-X$ExpressionMatrix
    nc = ncol(exp_matrix)
    ng = nrow(exp_matrix)
    discreteMatrix = array(0, dim = c(ng, nc, nbins))
    TEST <- apply(X = exp_matrix,FUN = function(x){
      bs(
        (x),
        df = 3,
        degree = 2,
        intercept = T
      )} , MARGIN = 1)

    for (gene_count in 1:ng) {
      j = bs(
        (exp_matrix[gene_count,]),
        df = nbins,
        degree = 2,
        intercept = T
      )
      for (counterCondition in (1:nc)) {
          newJ <- cbind(rowSums(j[, 1:2]), j[, 3:dim(j)[2]])
          discreteMatrix[gene_count, counterCondition,] = newJ[counterCondition,]
      }
    }
    rownames(discreteMatrix) <- rownames(exp_matrix)
    colnames(discreteMatrix) <- colnames(exp_matrix)
    discreteMatrix<-discreteMatrix[,,rev(1:dim(discreteMatrix)[3])]
    # Note: this is a fuzzy coded expression matrix
    Disjunctive_Matrix <- discreteMatrix %>% Create_Disjunctive_Matrix
    X$Disjunctive_Matrix<-Disjunctive_Matrix
    return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)}
}





#' Scale the expression matrix from 0 to 1
#'
#' @param X MCXpress object containing an ExpressionMatrix
#' @param scaled if FALSE scaling is done using as maximum value the highest value in the entire matrix, if TRUE the scaling is done gene wise using the maximum expression of each gene
#'
#' @return Disjunctive Matrix with 2 bins for each gene.
#' @export
#'
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- Discretisation_Range_01(MCX64553)
#' MCX64553$Disjunctive_Matrix
Discretisation_Range_01 <- function(X, scaled=TRUE){
  if(scaled==TRUE){
  Scale01 <-X$ExpressionMatrix %>% apply(MARGIN = 1, FUN = function(x) {(x - min(x)) / (max(x) - min(x))}) %>% t #Scale from 0 to 1
  Scale01bin1 <-  Scale01 %>%  set_rownames(Scale01 %>%  rownames %>% paste0("-bin1"))
  Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%  rownames %>% paste0("-bin2"))
  }
  else{
  maximum     <- X$ExpressionMatrix %>% max
  minimum     <- X$ExpressionMatrix %>% min
  Scale01     <- (X$ExpressionMatrix - minimum)/(maximum-minimum)
  Scale01bin1 <-  Scale01 %>%  set_rownames(Scale01 %>%  rownames %>% paste0("-bin1"))
  Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%  rownames %>% paste0("-bin2"))
  }
  Disjunctive_Matrix  <-Scale01bin1 %>%  rbind(Scale01bin2)
  Disjunctive_Matrix  <-Disjunctive_Matrix[Disjunctive_Matrix %>% rownames %>%  sort,] %>%  t
  X$Disjunctive_Matrix<-Disjunctive_Matrix
  return(X)
}

ParDiscretisation_Range_01 <- function(X, scaled=TRUE, nproc=NULL){
  clproc=makeCluster(nproc)
  on.exit(stopCluster(clproc))
  if(scaled==TRUE){
  Scale01 <-X$ExpressionMatrix %>% parApply(cl=clproc ,MARGIN = 1, FUN = function(x) {(x - min(x)) / (max(x) - min(x))}) %>% t #Scale from 0 to 1
  Scale01bin1 <-  Scale01 %>%  set_rownames(Scale01 %>%  rownames %>% paste0("-bin1"))
  Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%  rownames %>% paste0("-bin2"))
  }
  else{
  maximum     <- X$ExpressionMatrix %>% max
  minimum     <- X$ExpressionMatrix %>% min
  Scale01     <- (X$ExpressionMatrix - minimum)/(maximum-minimum)
  Scale01bin1 <-  Scale01 %>%  set_rownames(Scale01 %>%  rownames %>% paste0("-bin1"))
  Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%  rownames %>% paste0("-bin2"))
  }
  Disjunctive_Matrix  <-Scale01bin1 %>%  rbind(Scale01bin2)
  Disjunctive_Matrix  <-Disjunctive_Matrix[Disjunctive_Matrix %>% rownames %>%  sort,] %>%  t
  X$Disjunctive_Matrix<-Disjunctive_Matrix

  return(X)
}

Create_Disjunctive_Matrix <- function(ExpressionMatrix) {
  cat('Creating Disjunctive Matrix...',  "\n")
  MCA_Matrix <- NULL
  for (i in 1:dim(ExpressionMatrix)[3]) {
    Bin_MCA <- ExpressionMatrix[, , i]
    rownames(Bin_MCA) <-
      paste (rownames(Bin_MCA), paste("-bin", i, sep = ""), sep = "")
    MCA_Matrix <- rbind(MCA_Matrix, Bin_MCA)
  }
  MCA_Matrix <- MCA_Matrix[order(rownames(MCA_Matrix)), ] %>% t()
  return(MCA_Matrix)
}
