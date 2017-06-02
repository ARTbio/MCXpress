#' Scale the expression matrix using splines
#'
#' @param X MCXpress object containing an ExpressionMatrix
#' @return Disjunctive Matrix with 2 bins for each gene.
#' @export
#'
#' @examples
#' MCX64553 <- Initialise_MCXpress(GSE64553)
#' MCX64553 <- discretisation_bsplines(MCX64553)
#' #MCX64553$Disjunctive_Matrix
discretisation_bsplines <- function(X)
{
    if (X %>% class %>% equals("MCXpress_object"))
    {
        exp_matrix <- X$ExpressionMatrix
        nc = ncol(exp_matrix)
        ng = nrow(exp_matrix)
        TEST <- apply(X = exp_matrix, FUN = function(x)
        {
            Bsplines <- bs((x), df = 3, degree = 2, intercept = TRUE)
            cbind(Bsplines[, 1:2] %>% rowSums, Bsplines[, 3])
        }, MARGIN = 1)
        bin2 <- TEST[1:nc, ] %>% set_colnames(paste0(TEST[1:nc,
            ] %>% colnames, "-bin2")) %>% set_rownames(exp_matrix %>%
            colnames())
        bin1 <- TEST[(nc + 1):(2 * nc), ] %>% set_colnames(paste0(TEST[1:nc,
            ] %>% colnames, "-bin1")) %>% set_rownames(exp_matrix %>%
            colnames())
        Disjunctive_Matrix <- cbind(bin1, bin2)
        Disjunctive_Matrix <- Disjunctive_Matrix[, Disjunctive_Matrix %>%
            colnames %>% sort]
        X$Disjunctive_Matrix <- Disjunctive_Matrix
        return(X)
    } else
    {
        errormessage <- "The input is not a MCXpress object, apply the function Initialise_MCXpress on your expression matrix first."
        stop(errormessage)
    }
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
#' MCX64553 <- discretisation_01(MCX64553)
#' #MCX64553$Disjunctive_Matrix
discretisation_01 <- function(X, scaled = TRUE)
{
    if (scaled == TRUE)
    {
        Scale01 <- X$ExpressionMatrix %>% apply(MARGIN = 1, FUN = function(x)
        {
            (x - min(x))/(max(x) - min(x))
        }) %>% t  #Scale from 0 to 1
        Scale01bin1 <- Scale01 %>% set_rownames(Scale01 %>% rownames %>%
            paste0("-bin1"))
        Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%
            rownames %>% paste0("-bin2"))
    } else
    {
        maximum <- X$ExpressionMatrix %>% max
        minimum <- X$ExpressionMatrix %>% min
        Scale01 <- (X$ExpressionMatrix - minimum)/(maximum -
            minimum)
        Scale01bin1 <- Scale01 %>% set_rownames(Scale01 %>% rownames %>%
            paste0("-bin1"))
        Scale01bin2 <- (1 - Scale01) %>% set_rownames(Scale01 %>%
            rownames %>% paste0("-bin2"))
    }
    Disjunctive_Matrix <- Scale01bin1 %>% rbind(Scale01bin2)
    Disjunctive_Matrix <- Disjunctive_Matrix[Disjunctive_Matrix %>%
        rownames %>% sort, ] %>% t
    X$Disjunctive_Matrix <- Disjunctive_Matrix
    return(X)
}
