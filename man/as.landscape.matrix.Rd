% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landscape.R
\name{as.landscape.matrix}
\alias{as.landscape.matrix}
\title{transfer matrix to landscape.}
\usage{
\method{as.landscape}{matrix}(x, states = NA)
}
\arguments{
\item{x}{A matrix object with factorial content.}

\item{levels}{A character vector giving the original sequence of states.}
}
\value{
A landscape object.
}
\description{
transfer matrix to landscape.
}
\details{
If no vector levels is specified the levels are coerced from the 
  unique states of the matrix, which might cause loss of states (e.g. if not 
  present on this landscape) or a wrong ordering of states. This is not
  recommended! Always provide a levels vector!
}

