
#' Initialisation of the MCXpress Object
#'
#' Create a MCXpressObject from an Expression Matrix
#'
#' @param X DESCRIPTION.
#' @param min_reads DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' @export
#'
Initialise_MCXpress <- function(X)
{
  if (X %>% is.matrix)
  {
    if ((X %>% rownames %>% is.null) | (X %>% colnames %>%
                                        is.null))
    {
      errormessage <- "Gene name should be the rownames of the matrix and Sample name the column name"
            stop(errormessage)
        } else
        {
            MCXpress <- list()
            X <- as.matrix((X[apply(X, 1, var) > 0, ]))
            X <- X[str_length(X %>% rownames) > 0, ]
            X <- X %>% subset(X %>% rownames %>% duplicated %>% not)
            colnames(X) <- gsub(pattern = " ", replacement = "_",
                x = colnames(X))
            colnames(X) <- gsub(pattern = "-", replacement = "_",
                x = colnames(X))
            colnames(X) <- gsub(pattern = "\\.", replacement = "_",
                x = colnames(X))
            MCXpress$ExpressionMatrix <- X
            class(MCXpress) <- "MCXpress"
            return(MCXpress)
        }
    } else
    {
        errormessage <- "The input is not a matrix"
        stop(errormessage)
    }
}



Df_To_Mat <- function(x) {
  Dframe <- x[,-1] %>% as.matrix %>% set_rownames(x[,1])
  return(Dframe)
}

Mat_To_Df <- function(x, namecol="Sample") {
  Mat <- x %>%  as.data.frame() %>% tibble::rownames_to_column(var = namecol) %>%  as_tibble
  return(Mat)
}

Vec_To_Df <- function(x, namecol=c("X1","X2")){
  Dframe <- x %>%  as.data.frame() %>% tibble::rownames_to_column() %>%  as_tibble
  return(Dframe)
}

Num_Axis_Broken_Stick <- function(x){
  vari<- x$MCA$explained_eigen_variance$Explained_Variance
  vari <- vari/100
  broken <- sum(1/1:(vari %>% length))/(vari %>% length)
  (vari>broken) %>%  which %>%  length
}

SC_Pseudo_Cluster <- function(x){
  x<- x$ExpressionMatrix %>%  colnames() %>%  set_names(.,.)
  return(x)
}


Num_Axis_Eigen_Distance <- function(X){
  eig <- X$MCA$eigen_value
  A <- (eig[1]-eig[length(eig)])/(1-eig %>% length)
  B <- eig[1]-A
  Affine <- function(x){(A*x+B-eig) %>% which.max()}
  return(Affine(1:(eig %>% length)))
}

