% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Discretisation_Methods.R
\name{discretisation_01}
\alias{discretisation_01}
\title{Scale the expression matrix from 0 to 1}
\usage{
discretisation_01(X, scaled = TRUE)
}
\arguments{
\item{X}{MCXpress object containing an ExpressionMatrix}

\item{scaled}{if FALSE scaling is done using as maximum value the highest value in the entire matrix, if TRUE the scaling is done gene wise using the maximum expression of each gene}
}
\value{
Disjunctive Matrix with 2 bins for each gene.
}
\description{
Scale the expression matrix from 0 to 1
}
\examples{
MCX64553 <- Initialise_MCXpress(GSE64553)
MCX64553 <- discretisation_01(MCX64553)
#MCX64553$disjunctive_matrix
}
