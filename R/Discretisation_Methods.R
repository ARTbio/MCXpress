Discretisation_Bsplines <- function(X, nbins = 2) {
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

Discretisation_Range_01 <- function(X) {
  Scale01 <-X$ExpressionMatrix %>% apply(MARGIN = 1, FUN = function(x) {(x - min(x)) / (max(x) - min(x))}) %>% t #Scale from 0 to 1
  Scale01bin1 <-  Scale01 %>%  set_rownames(Scale01 %>%  rownames %>% paste0("-bin1"))
  Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%  rownames %>% paste0("-bin2"))
  Disjunctive_Matrix <-Scale01bin1 %>%  rbind(Scale01bin2)
  Disjunctive_Matrix<-Disjunctive_Matrix[Disjunctive_Matrix %>% rownames %>%  sort,] %>%  t
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
