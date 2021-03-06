% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/musselbed.R
\docType{data}
\name{musselbed}
\alias{musselbed}
\title{Musselbed Disturbance Model}
\format{An object of class \code{ca_model} of length 7.}
\usage{
ca(l, musselbed, parms = p)
}
\arguments{
\item{r}{A numerical value. recolonisation of empty sites dependent on local 
density.}

\item{d}{A numerical value. probability of disturbance of occupied sites if 
at least one disturbed site is in the direct 4-cell neighborhood.}

\item{delta}{A numerical value. intrinsic disturbance rate.}
}
\description{
Model for spatial pattern of intertidal mussel beds with wave 
  disturbance. With three cell states: occupied by mussels (\code{"+"}), 
  empty sites (\code{"0"} and disturbed sites (\code{"-"}).
}
\details{
The model represents the spatial dynamics in mussel cover of rock substrate 
  in intertidal systems. The stochastic wave disturbances will most likely 
  remove mussels that are located next to a gap because of the losened byssal
  threads in their proximity. This causes a dynamic gap growth.

The model describes the process by simplifying the system into three 
  potential cell states: occupied by mussel ("+"),  empty but undisturbed 
  ("0"), and disturbed, bare rock with loose byssal threads ("-"). 
 
Mussel growth on empty cells is defined by parameter \code{r} multiplied by 
  the local density of mussels in the direct 4-cell neighborhood.  
 
Any cell occupied by mussels has an intrinsic chance of \code{delta} to 
  become disturbed from intrinsic cause, e.g. natural death or predation. 
  Additionally, wave disturbance will remove mussels and leave only bare 
  rock, i.e. disturbed sites, with probability \code{d} if at least one 
  disturbed cell is in the direct 4-cell neighborhood. This causes 
  disturbances to cascade through colonies of mussels.  
 
Disturbed sites will recover into empty sites with a constant rate of 1 per 
  year, i.e. on average a disturbed site becomes recolonisable within one 
  year after the disturbance happened.
}
\examples{

l <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = 50) # create initial landscape
p <- list(delta = 0.01, d = 09, r = 0.4)   # set parameters
r <- ca(l, musselbed, p, t_max = 100)    # run simulation 


}
\author{
Guichard, Halpin, et al. (2003)
}
\references{
Guichard, F., Halpin, P.M., Allison, G.W., Lubchenco, J. & Menge,
  B.A. (2003). Mussel disturbance dynamics: signatures of oceanographic 
  forcing from local interactions. The American Naturalist, 161, 889–904.
}
\seealso{
Other models: \code{\link{forestgap}},
  \code{\link{grazing}}, \code{\link{life}},
  \code{\link{livestock}}, \code{\link{predprey}}
}
\keyword{datasets}

