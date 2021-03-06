% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caspr.R
\name{mapping}
\alias{mapping}
\title{Create vector maps for cellular automata updating.}
\usage{
mapping(width, height = width, boundary = "periodic",
  i_matrix = matrix(c(0, 1, 0, 1, NA, 1, 0, 1, 0), ncol = 3, byrow = TRUE))
}
\arguments{
\item{width}{An integer number giving the width of the landscape or a landscape
object.}

\item{height}{An integer number. The height of the landscape. Defaults to
value given in width.}

\item{boundary}{A character value representing the type of boundary
conditions. Defaults to \code{"periodic"}. No other boundaries implemented,
yet.}

\item{i_matrix}{A matrix defining the interaction matrix to be evaluated.
Defaults to 4-cell neighborhood.}
}
\value{
Creates a map of the landscape object to translate it into a vector
  with boundaries (\code{x_with_border}) and another one to back-translate it
  to a vector without boundaries  (\code{x_to_evaluate}) into the global
  environment. This will automatically be called within \code{ca()}.
}
\description{
Create vector maps for cellular automata updating.
}
\details{
This function and the objects generated by it are necessary for the
  function \code{\link{neighbors}} to vectorise the calculation of
  neighbors. The interaction matrix can be flexibly defined as a matrix that
  has precisely one cell with value \code{NA}, marking the focal cell. Any
  number larger than 0 marks a cell that is to be taken into account for
  counting.
}

