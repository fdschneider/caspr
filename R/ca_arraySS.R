#' Run ca_snapsSS() over a parameter array.
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
#' @param height An integer number. Defaults to \code{width} if 
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
#' @param ... parameters handed over to function \code{ca_snapsSS()}. Recommended the use of 
#' nsnaps in \code{ca_snapsSS()} to specify the number of landscape replicate to be saved per 
#' combination of parameters in \code{parms}.
#'   
#' @return Returns a list containing two elements. Firt one is a dataframe with global and 
#' local cover for each state for each parameter value or combination of parameter values 
#' given in \code{parms}. Second one is a list containing the specified number of snapshots 
#' in steady state for each parameter values comination given in \code{parms}.
#' @details The function is used to create gradients along one parameter value 
#'   or an array of parameter values. It runs the simulation for each parameter 
#'   value or combination of parameter values while making use of a parallel 
#'   backend provided by \code{\link[foreach]{foreach}}.
#'   
#' Make sure to adapt the type of parallel backend to your computer
#' infrastructure (see
#' \href{https://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf}{package
#' vignette of the foreach -package}).
#' 
#' @return The function returns a dataframe with global and local cover for each
#'   state for each parameter value or combination of parameter values given in
#'   `parms`, as well as the seed used for the individual simulation run.
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
#'  library(doMC)
#' 
#'  registerDoMC(cores=7)
#'  
#'  musselgradient <- ca_array(musselbed, parms = p, init = c(0.6,0.2,0.2), t_max = 400, save = TRUE, file = "out/musselbed")
#'  
#'  registerDoSEQ()
#' 
#'  musselgradient
#'  plot(musselgradient[,"mean_cover.+"] ~ musselgradient$d, pch = 20)
#' }


ca_arraySS <- function(model, 
                     init, 
                     width = 50, height = width,
                     parms = model$parms, 
                     save = FALSE, 
                     filename = "sim", directory = "", 
                     salt = 345678, 
                     length.stat = 10,
                     ... ) {
  
  iterations <- expand.grid(parms)
  iterations <- cbind(ID = 1:nrow(iterations),iterations)
  set.seed(salt)
  iterations$seed <- as.integer(runif(length(iterations$ID),100000,999999))
  
  if(!getDoParRegistered() & length(iterations$ID) > 10 ) {
    stop("No parallel backend registered! This would take too much time!")
  }
  
  
  
  foreach(i = iterations$ID,
          .export = c("model", "parms", "init"), 
          .packages = c("caspr")) %dopar% {
            
            set.seed(iterations$seed[i])
            
            # get initial landscape
            if("landscape" %in% class(init) ) { 
              l <- init 
            } else {
              if(length(init) == length(model$states))  {
                l <- init_landscape(model$states, init, width, height) 
              } else { 
                stop("Provide valid initial cover vector!")
              }
            }
            
            
            # running the simulation for this iteration
            run <- ca_snapsSS(l, model, parms = iterations[i,, drop = TRUE], 
                              seed = iterations$seed[i], 
                              length.stat = length.stat, ...)
            
            # get summary for output
            out <- c(summary(run, length.stat = length.stat)$mean_cover, 
                     summary(run, length.stat = length.stat)$sd_cover)
            names(out) <- paste0(rep(c("mean_cover", "sd_cover"), each = length(model$states) ), "_",names(out))
            out <- as.data.frame(t(out))
            out$t_start <- as.integer(summary(run)$time[1])
            out$t_end <- as.integer(summary(run)$time[2])
            # ca() does not return a steady_state component anymore
            out <- data.frame(out, is_steady = run[['issteady']][summary(run)$time[2]]) 
            
            # Get snapshots
            Snaps <- run$landscapes
            
            if(save) {
              dir.create(file.path(directory, sub(basename(filename), "", filename)), showWarnings = FALSE)
              path <- normalizePath(file.path(directory, sub(basename(filename), "", filename))) 
              obj <- paste0(basename(filename), "_",    # define an unambiguous output object name
                            formatC(i, width= ceiling(log10(length(iterations$ID))), flag="0")
              )
              
              assign(obj, run) 
              eval(parse(text = paste0("save(", obj, ", file = '", file.path(path, obj) ,".Rd' )" )   ) )
              
              out$saved_in <- paste0(obj, ".Rd")
            }
            #saving snapshots
            
            return(list(out = out, Snaps = Snaps))
          } -> outR
  out<- lapply(outR, `[[`, "out")
  names<-names(out[[1]])
  out<- data.frame(matrix(unlist(out), nrow=length(out), byrow=T),stringsAsFactors=FALSE)
  names(out) <- names
  Snaps <- lapply(outR, `[[`, "Snaps")
  
  return(list (DatBif = cbind(iterations, out),Snaps) )
}


#' @export
carray <- ca_array
