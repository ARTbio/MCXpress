



#'Select most variable genes
#'
#'@param X A MCXpress class object
#'@param ngenes Number of genes to keep
#'@return A MCXpress object with an updated Expression_Matrix
#'@export
#'@seealso
select_most_variable_genes <- function(X, ngenes) {
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<-X$ExpressionMatrix
    means <- exp_matrix %>%  rowMeans
    vars  <- apply(exp_matrix, 1, var)
    cv2   <- vars / means ^ 2
    minMeanForFit <- unname(quantile(means[which(cv2 > .3)], 0.95))
    useForFit <- means >= minMeanForFit
    fit <-
      glmgam.fit(cbind(a0 = 1, a1tilde = 1 / means[useForFit]), cv2[useForFit])
    a0  <- unname(fit$coefficients["a0"])
    a1  <- unname(fit$coefficients["a1tilde"])
    fit$coefficients
    xg  <-
      exp(seq(min(log(means[means > 0])), max(log(means)), length.out = 1000))
    vfit  <- a1 / xg + a0
    df  <- ncol(exp_matrix) - 1
    afit  <- a1 / means + a0
    varFitRatio <- vars / (afit * means ^ 2)
    varorder  <- order(varFitRatio, decreasing = T)
    oed <- exp_matrix[varorder, ]
    X$ExpressionMatrix <- exp_matrix[varorder[1:ngenes],]
    return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)}
}


Select_Highest_Mean <- function(X, ngenes) {
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<- X$ExpressionMatrix
    High_Means <-
      (exp_matrix %>%  rowMeans %>% sort(decreasing = TRUE) %>%  names %>% na.omit)
    High_Means <-  High_Means[str_length(High_Means) > 0]
    High_Means <-  High_Means[1:ngenes]
    X$ExpressionMatrix<-exp_matrix[High_Means, ]
    return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)}
}

Select_Diptest <- function(X, pval = 0.1) {
  if(X %>%  class %>%  equals("MCXpress_object")){
    exp_matrix<- X$ExpressionMatrix
    dipPval <-
      as.numeric(lapply(1:dim(exp_matrix)[1], function(x)
        dip.test(as.numeric(exp_matrix[x, ]))$p.value))
    dipPval <- matrix(c(dipPval, seq(1:length(dipPval))), ncol = 2)
    dipPvalOrderedIndex <- dipPval %>% order
    dipPvalSigIndex <-
      dipPvalOrderedIndex[dipPval[dipPvalOrderedIndex] < pval]
    X$ExpressionMatrix<- exp_matrix[dipPvalSigIndex, ]
    return(X)}else{errormessage<-"The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)}
}
