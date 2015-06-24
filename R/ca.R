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
#'   still in transient dynamics (i.e. model did not reach steady state)
#' @param t_min Minimal number of timesteps before steadyness is evaluated and 
#'   model might be terminated if transitory dynamics are surpassed.
#' @param t_eval Timespan of moving window that is evaluated for the end of 
#'   transient dynamics.
#' @param steady Tolerance level for steadyness. Steady state is reached if the 
#'   difference in mean cover of the primary cell state (i.e. the first in the 
#'   vector provided in \code{model$states}) over two subsequent timespans of 
#'   length \code{t_eval} is smaller than \code{steady}.
#' @param saveeach Timespan between timesteps at which a full snapshot of the 
#'   landscape is saved into the output of the simulation.
#' @param ... Parameters handed over to update function in \code{model$update}.
#'   
#' @return The output is returned as a list object of class \code{ca_result},
#'   containing a full timeseries of global and local cover as well as snapshots
#'   of the landscape.
#'   \describe{
#'     \item{\code{$model}}{The entire model object used to generate this simulation 
#'        run, including the parameters at \code{$model$parms}}
#'     \item{\code{$time}}{Vector of timesteps.}
#'     \item{\code{$evaluate}}{Start and endpoint of steady-state evaluation period.}
#'     \item{\code{$cover}}{A list of cover timeseries for each state of the model.}
#'     \item{\code{$local}}{A list of local cover timeseries for each state of the 
#'        model.}
#'     \item{\code{$snaps}}{A vector of indices of saved snapshots.}
#'     \item{\code{$landscapes}}{A list of landscape objects at each point in 
#'        \code{$snaps}}
#'     \item{\code{$steadyness}}{Steadyness value: difference in mean cover of primary 
#'        state between first and second half of the evaluation period. }
#'   }
#' @details Runs iterations of the update function \code{model$update()} on the
#'   initial landscape \code{x} until a steady state is reached (as defined by
#'   the tolerance level \code{steady}), but max \code{t_max} timesteps. The
#'   function saves the full timeseries, i.e. a value for each timestep, of the
#'   global cover of each state as well as the average local cover of each
#'   state. Only for every \code{saveeach}th timestep, the full lattice is saved
#'   in a list within the output file (\code{result$snapshots} ).
#'   
#' @export
#' 
#' @examples 
#' 
#' # 1. run simulation and save a snapshot each 50 timesteps. plot timeseries and snapshots.
#' 
#' l <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2)) # create initial landscape
#' p <- list(r = 0.4, d = 0.9, delta = 0.01)   # set parameters
#' r <- ca(l, model = musselbed, parms = p, t_max = 500)    # run simulation 
#' plot(r)
#' 
#' par(mfrow= c(2,3))
#' lapply(r$landscapes, plot)
#' 
#' # 2. run simulation and save full landsape at each timestep. create animated gif.
#' 
#' l <- init_landscape(c("1","0"), c(0.6,0.4), 100)
#' r <- ca(l, model = life, t_max = 500, t_eval = 500, saveeach = 1)
#' animate(r, "life01.gif")

ca <- function(x, model = musselbed, parms = "default", t_min = 500, t_eval = 200, saveeach = 50, steady = 0.00001, t_max = 1000, ... )  {
  
  # checking fo valid input
  ## parms
  if(parms[1] != "default") { 
    if(!all(names(model$parms) %in% names(parms))) {
      stop(paste("missing parameter(s): ", names(parms)[all(names(model$parms) %in% names(parms))] , "! please specify!"    ) )
    }
  } else {
    parms <- model$parms
    warning("you did not specify 'parms'! using default parameters of model!")
  }
  
  ## states
  if(!all(x$states %in% model$states)) stop("invalid cell states specified in landscape object!")
 #notworking!  if(any(!model$states %in% x$states)) warning("one or several of the model's default cell states are not present in initial landscape object!") 
  
  ## is interaction matrix provided? if not defaults to four cell neighborhood.   
  if(is.null(model$interact)) {
    I <- matrix(c(0, 1, 0,  1, NA, 1, 0, 1, 0), ncol = 3, byrow = TRUE) 
  } else {
    I <- model$interact 
  }
    
  mapping(x$dim[1], x$dim[2], i_matrix = I )
  xstats <- summary(x)
  states <- levels(x$cells)
  
  # calculate timesteps for snapshots:
  if(t_max <= t_min) {t_eval = t_max; t_min = t_max}
  if(t_eval+1 > t_min) { t_min <- t_eval+1 } # increase t_min if t_eval is too long 

  snaps <- as.integer(t_min) - rev(seq(0, as.integer(t_eval),as.integer(saveeach)))
  
  # initialise result object:
  
  result <- list()  # generate list object
  result$model <- model
  result$model$parms <- parms
  result$time <- seq(0, t_min) # add vector of realized timesteps
  result$evaluate <- c(t_min, t_min+2*t_eval)+1
  result$steadyness <- 1
  
  result$cover <- as.data.frame(t(xstats$cover))
  result$cover <- result$cover[rep(1, t_min+1),] # preallocating memory
  
  result$local <- as.data.frame(t(xstats$local))
  result$local <- result$local[rep(1, t_min+1),] # preallocating memory
  
  result$snaps <- snaps
  result$landscapes <- list()
  result$landscapes <- lapply(1:length(snaps), function(i) x) # preallocating memory
  
  # --------------- start simulation -----------------
  
  parms_temp <- parms   # set temporary parameter object
  
  # initialise simulation variables: 
  x_old <- x  # ghost matrix at t_i
  x_new <- x_old
  steadyness <- 1  # check value for stability
  i = 0  # iterator for simulation timesteps
  
  # starting iterations:
  while(steadyness > steady & i <= as.integer(t_max) ) {
    
    i <- i +1  # increase iterator
    
    # call update function:
    
    model$update(x_old, parms, ...) -> x_new
    
    # save stats of new landscape
    
    xstats <- summary(x_new)
    result$cover[i,] <- xstats$cover
    result$local[i,] <- xstats$local
    
    # replace ghost matrix for next iteration
    
    x_old <- x_new 
    
    # save landscape if snapshot
    
    if(i %in% snaps) {  
      result$landscapes[[match(i, snaps)]] <- x_new
    }
    
    # calculate steadyness and update snaps if not steady
    
    if(i %in% snaps & i >= t_min) { # if we are over the minimal timespan 
      # get vector of t_eval timesteps previous to the last t_eval timesteps
      t_1 <- (i-as.integer(t_eval)):(i-0.5*as.integer(t_eval)) 
      # get vector of the last t_eval timesteps 
      t_2 <- (i-0.5*as.integer(t_eval)+1):(i) 
      
      if(result$cover[[1]][i] > 0) { 
        # calculate steadyness, i.e. difference between the mean cover 
        #    in the two evaluation periods t_1 and t_2
        steadyness <- (abs(mean(result$cover[[1]][t_1]) - mean(result$cover[[1]][t_2]))) /
          (mean(result$cover[[1]][t_1])) 
      } else {
        # set steadyness to 0 if cover is 0, immediate stop of simulation
        steadyness <- 0 
      }
      
      # update snaps vector if not steady
      
      if(steadyness > steady & i == max(snaps) & i+saveeach <= t_max ) {
        snaps[match(min(snaps), snaps)]  <- max(snaps)+as.integer(saveeach)
      }
      
      # save timesteps of last evaluation
      
      result$evaluate <- c(min(t_1), max(t_2))
     
    }
    result$time[i] <- i # save timestep to results
    
  } 
  # ------------ end simulation -----------------
  
  result$steadyness <- steadyness
  result$snaps <- snaps
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
#' 
#' @export

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
#' @export

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


