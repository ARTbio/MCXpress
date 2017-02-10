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
    # Note: this is a fuzzy coded expression matrix
    Disjunctive_Matrix <- discreteMatrix %>% Create_Disjunctive_Matrix
    X$Disjunctive_Matrix<-Disjunctive_Matrix
    return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)}
}

Discretisation_Range_01 <- function(X) {
  exp_matrix<-X$ExpressionMatrix
  Scale01 <-
    (apply(
      X = exp_matrix,
      1,
      FUN = function(x) {
        (x - min(x)) / (max(x) - min(x))
      }
    )) #Scale from 0 to 1
  Scale01bin1 <-
    t(Scale01) #Put in two different bins (bin1=probability of being expressed bin2=not expressed)
  Scale01bin2 <- t(1 - Scale01)
  discreteMatrix <-
    (abind(Scale01bin1, Scale01bin2, along = 3)) #create 3 dim array dim1=genes dim2=cells dim3=bins
  rownames(discreteMatrix) <-
    rownames(exp_matrix) #assign the coresponding row and colnames
  colnames(discreteMatrix) <- colnames(exp_matrix)
  Disjunctive_Matrix <- discreteMatrix %>% Create_Disjunctive_Matrix
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
