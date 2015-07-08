#' Run ca() over parameter array.
#' 
#' @param model A valid object of class 'ca_model'. Defaults to musselbed. Valid
#'   values are: \code{musselbed}, \code{grazing}.
#' @param init A valid numerical vector of initial cover. Length must be equal 
#'   to \code{length(model$states)}.
#' @param parms A list of parameters with one or several parameters containing a
#'   vector of parameter values. Those will be combined full-factorially using 
#'   \code{expand.grid()}
#' @param width An integer number. Defaults to 50 and is ignored if a landscape 
#'   object is provided in \code{init}.
#' @param height An integer number. Defaults to 50 or to \code{width} if 
#'   provided and is ignored if a landscape object is provided in \code{init}.
#' @param salt An integer number. Used as seed when generating seeds for each 
#'   single simulation run. Those generated seeds are returned in the output of 
#'   the function.
#' @param save Logical. If TRUE, each individual run is saved to a file. They 
#'   can be loaded into the current working environment using 
#'   \code{load("filename")}.
#' @param filename A character vector specifying the root of the filename of the
#'   saved files, including a relative path. It will be extended by an
#'   individual iteration ID and the fileending ".Rd". Note: If running in
#'   parallel on a cluster, output files will be saved in the workers home
#'   directory!
#' @param ... parameters handed over to function \code{ca()}.
#'   
#' @return Returns a dataframe with global and local cover for each state for 
#'   each parameter value or combination of parameter values given in 
#'   \code{parms}.
#' @details The function is used to create gradients along one parameter value 
#'   or an array of parameter values. It runs the simulation for each parameter 
#'   value or combination of parameter values while making use of a parallel 
#'   backend provided by \code{\link[foreach]{foreach}}.
#'   
#' @export
#' @import foreach
#' @examples 
#' \dontrun{
#'   p <- list(
#'     r = 0.4, # recolonisation of empty sites dependent on local density
#'     d = seq(0,1,0.1), # wave disturbance 
#'     delta = 0.01, # intrinsic disturbance rate
#'     replicates = 1:2  # repeat the same parameter set twice
#'   ) 
#'   
#' # provides parallel backend
#'  library(foreach)
#'  library(doSNOW)
#' 
#'  workerlist <- c(rep("localhost", times = 7)) #adjust to your computers' capacities
#' 
#'  cl <- makeSOCKcluster(workerlist, outfile='out_messages.txt')
#' 
#'  registerDoSNOW(cl)
#'  
#'  musselgradient <- carray(musselbed, parms = p, init = c(0.6,0.2,0.2), t_max = 400, save = TRUE, file = "out/musselbed")
#'  
#'  stopCluster(cl)
#' 
#'  musselgradient
#'  plot(musselgradient[,"mean_cover.+"] ~ musselgradient$d, pch = 20)
#' }


carray <- function(model, init, parms = model$parms, save = FALSE, filename = "sim", directory = getwd(), width = 50, height = width, salt = 345678, ...) {
  
  iterations <- expand.grid(parms)
  iterations <- cbind(ID = 1:nrow(iterations),iterations)
  set.seed(salt)
  iterations$seed <- as.integer(runif(length(iterations$ID),100000,999999))
  
  if(!getDoParRegistered() & length(iterations$ID) > 10 ) {
    stop("No parallel backend registered! This would take too much time!")
  }
  
  foreach(i = iterations$ID, .combine = rbind, .packages = c("caspr")) %dopar% {
    
    set.seed(iterations$seed[i])
    
    # get initial landscape
    if(length(init) == length(model$states))  {
      l <- init_landscape(model$states, init, width, height) 
    } else { 
      stop("Provide valid initial cover vector!")
    }
    
    # running the simulation for this iteration
    run <- ca(l, model, parms = iterations[i,, drop = TRUE], seed = iterations$seed[i], ...)
    
    # get summary for output
    out <- unlist(summary(run)[-1])
    
    if(save) {
      obj <- paste0(filename, "_",    # define an unambiguous output object name
                    formatC(i, width= ceiling(log10(length(iterations$ID))), flag="0")
      )
      
      assign(obj, run) 
      eval(parse(text = paste0("save(",obj, ", file = '", paste0(obj, ".Rd"), "')")   ) )
      
      out <- c(out, saved_in =  paste0(obj, ".Rd"))
    }

    return(out)
  } -> out
  
  
  return(cbind(iterations, out) )
}
