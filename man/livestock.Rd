% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/livestock.R
\docType{data}
\name{livestock}
\alias{livestock}
\title{Livestock resilience model}
\format{An object of class \code{ca_model} of length 6.}
\usage{
ca(l, livestock, parms = p)
}
\arguments{
\item{b}{environmental quality}

\item{m}{intrinsic plant mortality, i.e. inverse of av. lifespan}

\item{r}{max. regeneration rate of plants}

\item{K}{carrying capacity of the landscape (i.e. global competition)}

\item{f}{local facilitation}

\item{c}{local competition}

\item{L}{livestock density}

\item{a}{attack rate of livestock}

\item{h}{handling time of livestock}

\item{q}{hill-exponent of livestock}

\item{p}{associational resistance against livestock grazing}

\item{v}{attractant-decoy effect of plants to livestock}
}
\description{
A spatially-explicit model of resilience in arid rangelands.
}
\details{
An unpublished, generalized model of positive and negative local feedbacks in arid rangelands, including mechanisms such as local facilitation, competition, associational resistance and attractant-decoy.
}
\author{
F.D. Schneider
}
\seealso{
Other models: \code{\link{forestgap}},
  \code{\link{grazing}}, \code{\link{life}},
  \code{\link{musselbed}}, \code{\link{predprey}}
}
\keyword{datasets}

