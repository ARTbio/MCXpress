% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Selection_Methods.R
\name{select_diptest}
\alias{select_diptest}
\title{Hartigan's diptest for gene selection}
\usage{
select_diptest(X, pval = 0.05)
}
\arguments{
\item{X}{MCXpress Object with an Expression Matrix}

\item{percentage}{numeric between 0 and 1 indicating the percentages of cells that must express the genes.}

\item{threshold}{gene expression lower threshold.}
}
\value{
A filtered Expression Matrix in the MCXpressObject.
}
\description{
Hartigan's diptest for gene selection
}
\examples{
MCX64553 <- Initialise_MCXpress(GSE64553)
MCX64553 <- filter_outlier(MCX64553, percentage = 0.05, threshold = 3)
}
