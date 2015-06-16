#' Run ca() over parameter array.
#' 
#' @param x A landscape object.
#' @param parms A list of parameters with one or several parameters containing a vector of parameter values. Those will be combined full-factorially using \code{expand.grid()}
#' @param model A valid object of class 'ca_model'. Defaults to musselbed. Valid values are: \code{musselbed}, \code{grazing}.
#' @param ... parameters handed over to function \code{ca()}.
#' @return Returns a dataframe with global and local cover for each state for each parameter value or combination of parameter values given in \code{parms}. 
#' @details 
#'   The function is used to create gradients along one parameter value or an array of parameter values. 
#' @note  \strong{ToDo:} 
#' 
#'   - initialise with different landscape objects
#'   
#'   - implement selector for output, e.g. indicators = TRUE
#' @examples 
#' \dontrun{
#'   p <- list(
#'     r = 0.4, # recolonisation of empty sites dependent on local density
#'     d = seq(0,1,0.1), # wave disturbance 
#'     delta = 0.01 # intrinsic disturbance rate
#'   ) 
#'   
#' # provides parallel backend
#'  library(foreach)
#'  library(doSNOW)
#' 
#'  workerlist <- c(rep("localhost", times = 7)) 
#' 
#'  cl <- makeSOCKcluster(workerlist, outfile='out_messages.txt')
#' 
#'  registerDoSNOW(cl)
#' 
#'  musselgradient <- ca_array(initial, p, model = musselbed)
#' 
#'  stopCluster(cl)
#' 
#'  musselgradient
#' }


carray <- function(x, parms, ...) {
  require(foreach)
  iterations <- expand.grid(parms)
  iterations <- cbind(ID = 1:nrow(iterations),iterations)
  
  foreach(i = iterations$ID, .combine = rbind) %dopar% {
    
    source("R/functions.r")
    
    run <- ca(x, iterations[i,, drop = TRUE], ...)
    
    return(summary(run)$mean_cover)
  } -> out
  
  
  return(cbind(iterations , out) )
}
