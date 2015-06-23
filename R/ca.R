#' Run cellular automata simulation.
#' 
#' @param x A landscape object.
#' @param model A valid object of class 'ca_model'. Defaults to musselbed. Valid
#'   values are: \code{musselbed}, \code{grazing}.
#' @param parms A list of parameters with one or several parameters containing a
#'   vector of parameter values. Those will be combined full-factorially using 
#'   \code{expand.grid()}. Will be checked against template parameters in 
#'   \code{model}. If not provided, template parameters will be used.
#' @param t_max Maximal number of timesteps. Model will be terminated even if 
#'   unstable.
#' @param t_min Minimal number of timesteps before stability is evaluated and 
#'   model might be terminated if transitory dynamics are surpassed.
#' @param t_eval Timespan that is evaluated for stability.
#' @param isstable Tolerance level for stability. Stability is reached if the 
#'   difference in mean cover of the primary cell state (i.e. the first in the 
#'   vector provided in \code{model$states}) over two subsequent timespans of 
#'   length \code{t_eval} is smaller than \code{isstable}.
#' @param saveeach Timespan between timesteps at which a full snapshot of the 
#'   landscape is saved into the output of the simulation.
#' @param ... Parameters handed over to update function in \code{model$update}.
#'   
#' @return The output is returned as a list object of class \code{ca_result},
#'   containing a full timeseries of global and local cover as well as snapshots
#'   of the entire landscape (stored in \code{out$landscapes}).
#'   
#' @details Runs iterations of the update function \code{model$update()} on the
#'   initial landscape \code{x} until a steady state is reached (as defined by
#'   the tolerance level \code{isstable}), but max \code{t_max} timesteps. The
#'   function saves the full timeseries, i.e. a value for each timestep, of the
#'   global cover of each state as well as the average local cover of each
#'   state. Only for every \code{saveeach}th timestep, the full lattice is saved
#'   in a list within the output file (\code{result$snapshots} ).
#'   
#' 

ca <- function(x, model = musselbed, parms = "default", t_max = 1000, t_min = 500, t_eval = 200, isstable = 0.00001, saveeach = 50, ... )  {
  
  if(parms[1] != "default") { 
    if(!all(names(model$parms) %in% names(parms))) {
      stop(paste("missing parameter(s): ", names(parms)[all(names(model$parms) %in% names(parms))] , "! please specify!"    ) )
    }
  } else {
    parms <- model$parms
    warning("you did not specify 'parms'! using default parameters of model!")
  }
  
  mapping(x$dim[1], x$dim[2])
  xstats <- summary(x)
  states <- levels(x$cells)
  
  # calculate timesteps for snapshots:
  n_snaps <- length(seq(t_min-2*t_eval,t_min,saveeach)) # number of snapshots over t_eval
  t_snaps <-  ceiling(t_max/(saveeach*n_snaps))*(n_snaps*saveeach)
  snapshots <- data.frame(time = seq(t_min-2*t_eval, t_max, saveeach), i= seq(t_min-2*t_eval, t_max, saveeach)+1, pos = c(rep(1:n_snaps,2),1) )  
  ### BUG: this is not dynamically adjusting to t_min, t_max, t_eval, saveeach !!!
  
  
  # initialise result object:
  result <- list()  # generate list object
  result$model <- model
  result$model$parms <- parms
  result$time <- seq(0, t_min) # add vector of realized timesteps
  result$evaluate <- c(t_min, t_min+2*t_eval)+1
  
  result$cover <- as.data.frame(t(xstats$cover))
  result$cover <- result$cover[rep(1, t_min+1),] # preallocating memory
  
  result$local <- as.data.frame(t(xstats$local))
  result$local <- result$local[rep(1, t_min+1),] # preallocating memory
  
  result$snapshots <- snapshots
  result$landscapes <- list()
  result$landscapes <- lapply(1:n_snaps, function(i) x) # preallocating memory
  
  # --------------- simulation -----------------
  
  parms_temp <- parms   # set temporary parameter object
  
  # initialise simulation variables: 
  x_old <- x  # ghost matrix at t_i
  x_new <- x_old
  stability <- 1  # check value for stability
  i = 0  # iterator for simulation timesteps
  
  # starting iterations:
  while(stability > isstable & i <= t_max ) {
    
    i <- i +1  # increase iterator
    
    model$update(x_old, parms, ...) -> x_new
    
    xstats <- summary(x_new)
    result$cover[i,] <- xstats$cover
    result$local[i,] <- xstats$local
    
    x_old <- x_new # replace ghost matrix for next iteration
    
    if(i %in% snapshots$i) {  
      result$landscapes[[match(i, snapshots$i)]] <- x_new
    }
    
    
    if(i > t_min) { # if we are over the minimal timespan 
      t_1 <- (i-2*t_eval):(i-t_eval)-1 # vector of t_eval timesteps previous to the last t_eval timesteps
      t_2 <- (i-t_eval):(i) # vector of the last t_eval timesteps 
      
      if(result$cover[[1]][i] > 0) { 
        stability <- (abs(mean(result$cover[[1]][t_1]) - mean(result$cover[[1]][t_2])))/(mean(result$cover[[1]][t_1])) # calculate stability, i.e. difference between the mean cover in the two evaluation periods t_1 and t_2 
      } else {
        stability <- 0 # set stability to 0 if cover is 0, immediate stop of simulation
      }
      result$evaluate <- c(min(t_1), max(t_2))
      result$time[i] <- i # save timestep to results
      
    }
    
  } 
  
  class(result) <- "ca_result"
  return(result)
}



print.ca_result <- function(x) {
  cat("Model run of", x$model$name, " over ", tail(x$time,1)," timesteps. \n")
  cat("final cover:")
}


#' plot method for 'ca_result' objects
#' 
#' @param x An output object obtained from running function \code{ca()}.
#' @param plotstates A logical vector of the same length as 
#'   \code{x$model$states}. TRUE activates plotting of other cell state populations. 
#' @param cols A vector of valid colors to replace the default color vector.
#' @param lwd Line width for timeseries plot.
#' @param ... Parameters handed to standard plot function, e.g. `bty`, `xlab`.
#'   See \link{par}.
#' @details This prompts a (series of) summary plot for the simulation run. By
#' default this is a timeseries of the primary cell state (i.e. the first state
#' given in \code{x$model$states}).

plot.ca_result <- function(x, plotstates = c(TRUE, rep(FALSE, length(x$model$states)-1)), snapshots = FALSE, cols = x$model$cols , lwd = 1, ...) {
  
  if(snapshots) {
    
  }
  
  plot(NA,NA,  
       type = "l", 
       col = x$model$cols[1], 
       xlab = "time", ylab = paste("cover of", x$model$states[1]) , 
       xlim = c(1, max(x$time)), ylim = c(0,1), 
       ...)
  
  if(length(plotstates) == length(x$model$states)  ) {
    for(i in (1:length(plotstates))[plotstates]) {
      lines(x$time,x$cover[[i]], col = cols[i], lwd = lwd)
      
    }
    
  }
  
}



#' summary method for `ca_result`
#'
#' @param x 
#'
#' @return Returns a list \code{out} containing the model name, the final time, the mean cover  and standard deviation after transitory dynamics. 
#'    
#' 

summary.ca_result <- function(x) {
  out <- list()
  class(out) <- "ca_summary"
  eval <- x$evaluate[1]:x$evaluate[2]
  out$name <- x$model$name
  out$time <- c(min(x$time), max(x$time))
  out$mean_cover <- colMeans(x$cover[eval,])
  out$sd_cover <- sapply(x$cover[eval,], sd)
  
  return(out)
}


#' Transfer output of ca() into a list of matrices. 
#'
#' @param x 
#'
#' @return a list of matrices.
#' @export
#'

as.list.ca_result <- function(x) {
  lapply(x$landscapes, as.matrix)
}


#' Transfer output of ca() into an array. 
#'
#' @param x 
#'
#' @return an array with dimensions 'width', 'height' and 'snaps'.
#' @export
#'

as.array.ca_result <- function(x) {
  width <- x$landscapes[[1]]$dim[1]
  height <- x$landscapes[[1]]$dim[2]
  snaps <- length(x$landscapes)
  array(unlist(lapply(x$landscapes, as.matrix)), dim = c(width, height, snaps), dimnames = c("width", "height", "snaps"))
}
