% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Discretisation_Methods.R
\name{discretisation_bsplines}
\alias{discretisation_bsplines}
\title{Scale the expression matrix using splines}
\usage{
discretisation_bsplines(X)
}
\arguments{
\item{X}{MCXpress object containing an ExpressionMatrix}
}
\value{
Disjunctive Matrix with 2 bins for each gene.
}
\description{
Scale the expression matrix using splines
}
\examples{
MCX64553 <- Initialise_MCXpress(GSE64553)
MCX64553 <- discretisation_bsplines(MCX64553)
#MCX64553$disjunctive_matrix
}
