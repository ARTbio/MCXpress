

#'Select most variable genes
#'
#'@param X A MCXpress class object
#'@param ngenes Number of genes to keep, if not specified will take an  optimal number.
#'@return A MCXpress object with an updated Expression_Matrix
#'@examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#'@export
select_most_variable_genes <- function(X, ngenes = NULL)
{
    if (X %>% class %>% equals("MCXpress_object"))
    {
        exp_matrix <- X$ExpressionMatrix
        means <- exp_matrix %>% rowMeans
        vars <- exp_matrix %>% apply(1, var)
        cv2 <- vars/means^2
        minMeanForFit <- unname(quantile(means[which(cv2 > 0.3)],0.95))
        useForFit <- means >= minMeanForFit
        fit <- glmgam.fit(cbind(a0 = 1, a1tilde = 1/means[useForFit]),
                          cv2[useForFit])
        a0 <- unname(fit$coefficients["a0"])
        a1 <- unname(fit$coefficients["a1tilde"])
        fit$coefficients
        xg <- exp(seq(min(log(means[means > 0])), max(log(means)),
                      length.out = 1000))
        vfit <- a1/xg + a0
        df <- ncol(exp_matrix) - 1
        afit <- a1/means + a0
        varFitRatio <- vars/(afit * means^2)
        varorder <- order(varFitRatio, decreasing = TRUE)
        oed <- exp_matrix[varorder, ]
        if (ngenes %>% is.null())
        {
          X$ExpressionMatrix <- exp_matrix[varorder[1:((varFitRatio >
                                                          1) %>% which %>% length)], ]
        } else
        {
          X$ExpressionMatrix <- exp_matrix[varorder[1:ngenes],
                                           ]
        }
        return(X)
    } else
    {
      errormessage <- "The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
      stop(errormessage)
    }
}

#' Hartigan's diptest for gene selection
#'
#'
#'
#' @param X MCXpress Object with an Expression Matrix
#' @param percentage numeric between 0 and 1 indicating the percentages of cells that must express the genes.
#' @param threshold gene expression lower threshold.
#' @return
#' A filtered Expression Matrix in the MCXpressObject.
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' @export
select_diptest <- function(X, pval = 0.05)
{
  if (X %>% class %>% equals("MCXpress_object"))
  {
    exp_mat <- X$ExpressionMatrix
    Dip_Test <- function(x){
      x %>%
        as.numeric() %>%
        diptest::dip.test() %>%
        use_series(p.value)
    }
    dipPval <-exp_mat  %>% apply(MARGIN = 1,FUN = Dip_Test)
    X$ExpressionMatrix <- exp_mat[(dipPval < pval) %>% which,]
    return(X)
    } else
  {
    errormessage <- "The input is not an object of class MCXpress, apply the function Initialise_MCXpress on your expression matrix first."
    stop(errormessage)
  }
}

#' Filter genes that are only expressed in a few cells
#'
#'
#'
#' @param X MCXpress Object with an Expression Matrix
#' @param percentage numeric between 0 and 1 indicating the percentages of cells that must express the genes.
#' @param threshold gene expression lower threshold.
#' @return
#' A filtered Expression Matrix in the MCXpressObject.
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
#' @export
filter_outlier <- function(X, percentage = 0.05, threshold = 1)
{
    quant <- function(x)
    {
        (x %>% quantile(prob = 1 - percentage)) > threshold
    }
    noout <- X$ExpressionMatrix %>% apply(MARGIN = 1, FUN = quant) %>%
        which
    X$ExpressionMatrix <- X$ExpressionMatrix[noout, ]
    return(X)
}
