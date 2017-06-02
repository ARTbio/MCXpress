
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
Initialise_MCXpress <- function(X, min_reads = NULL)
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
            X <- X %>% subset(X %>% rownames %>% duplicated %>%
                not)
            colnames(X) <- gsub(pattern = " ", replacement = "_",
                x = colnames(X))
            colnames(X) <- gsub(pattern = "-", replacement = "_",
                x = colnames(X))
            colnames(X) <- gsub(pattern = "\\.", replacement = "_",
                x = colnames(X))
            if (min_reads %>% is.numeric())
            {
                X <- X[rowSums(X) > min_reads, ]
            }
            MCXpress$ExpressionMatrix <- X
            class(MCXpress) <- "MCXpress_object"
            return(MCXpress)
        }
    } else
    {
        errormessage <- "The input is not a matrix"
        stop(errormessage)
    }
}


