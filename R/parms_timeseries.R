
#' Create sequence of parameters over time.
#'
#' @param parms A list of parameters. Each list entry can be either a single value or a vector representing the sequence of parameter values over time (in years). This is inflated by the requested time in t_max.  
#' @param t_max An integer number. The number of years to evaluate the timeseries. 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' p <- list(  
#'   r = 1.0,  # max. regeneration rate of plants
#'   b = seq(0.6,0.1, length = 50),  # environmental quality
#'   sigma = 0.1, # random annual variation of environmental quality
#'   f = 0.6,  # local facilitation
#'   alpha = 0, # water runoff
#'   K = 0.9, # carrying capacity of the system
#'   c = 0.2, # local competition
#'   m = 0.05, # intrinsic mortality of plants (inverse of av. lifespan)
#'   v = 0.0, # attractant-decoy
#'   p = 0.99, # associational resistance
#'   L = rep(c(0,10), each = 2), # Livestock density
#'   q = 0, # hill exponent of functional response
#'   h = 50, # handling time 
#'   a = 0.2 # attack rate of livestock
#' ) 
#' 
#' 
#' plist <- parms_timeseries(p, 50)
#' 
#' l <- init_landscape(c("1","0"), cover = c(0.5,0.5), width = 100)
#' 
#' run <- ca(l, livestock, plist, t_max = 50, saveeach = 5)
#' plot(run)

parms_timeseries <- function(parms, t_max) {
  
  parms$i <- 1:t_max
  parms_list <- t(do.call(rbind.data.frame, parms ))
  row.names(parms_list) <- parms$i
  class(parms_list) <- "parms_timeseries"
  return(parms_list)
  
}
