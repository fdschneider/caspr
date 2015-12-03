#' Run cellular automata simulation.
#' 
#' @param x A landscape object.
#' @param model A valid object of class 'ca_model'. Defaults to musselbed. Valid
#'   values are: \code{musselbed}, \code{grazing}.
#' @param parms A list of parameters with one or several parameters containing a
#'   vector of parameter values. Those will be combined full-factorially using 
#'   \code{expand.grid()}. Will be checked against template parameters in 
#'   \code{model}. If not provided, template parameters will be used.
#'  @param t_min Minimum timespan for starting model evaluation of the steady state
#' @param t_max Maximal number of timesteps. Model will be terminated even if 
#'   still in transient dynamics (i.e. model did not reach steady state)
#' @param t_eval Timespan of moving window that is evaluated for the end of 
#'   transient dynamics.
#' @param nsnaps number of independent snapshots to be saved after steady state is achieved.
#' Minimum is 3. Default is 10.
#' @param stopifsteady Binary parameter, defaults to FALSE. If TRUE, the
#'   function provided in parameter \code{steady} will be applied to test in
#'   each timestep if steady state is reached.
#' @param steady A function returning TRUE or FALSE, taking exactly the
#'   parameters \code{i}, \code{result}, \code{steadyparms}. By default the
#'   function returns TRUE if the difference in mean cover of the primary cell
#'   state (i.e. the first in the vector provided in \code{model$states}) over
#'   two subsequent timespans of length \code{steadyparms$t_eval} is smaller
#'   than \code{steadyparms$accept}.
#' @param steadyparms a list of parameters that are required by the function
#'   provided in \code{steady}.
#' @param seed An integer number serving as seed for random number generation.
#' @param plotting A binary variable. If TRUE, simulation is plotted into an animated gif.
#' @param filename A character string. Filename of animated gif (defaults to "modelrun.gif") which will be placed in current working directory.
#' 
#'   
#'   If not provided global seeds of R apply.
#' @param ... Parameters handed over to update function in \code{model$update}.
#'   
#' @return The output is returned as a list object of class \code{ca_result}, 
#'   containing a full timeseries of global and local cover as well as snapshots
#'   of the landscape.
#'   
#'   \describe{ 
#'   \item{\code{$model}}{The entire model object
#'   used to generate this simulation run, including the parameters at
#'   \code{$model$parms}}
#'   \item{\code{$time}}{Vector of timesteps.} 
#'   \item{\code{$issteady}}{Binary vector reporting for each timestep if the criterion for steady state was fulfilled. The criterion can be customized by adjusting parameter \code{steady}.} 
#'   \item{\code{$cover}}{A list of cover timeseries for each state of the model.} 
#'   \item{\code{$local}}{A list of local cover timeseries for each
#'   state of the model.} 
#'   \item{\code{$snaps}}{A vector of indices of saved snapshots.} 
#'   \item{\code{$landscapes}}{A list of landscape objects at each
#'   point in \code{$snaps}} \item{\code{$issteady}}{A binary vector of the
#'   returned values of function \code{steady}} 
#'   }
#'   
#' @details Runs iterations of the update function \code{model$update()} on the 
#'   initial landscape \code{x} until a \code{t_max} is reached. The function
#'   saves the full timeseries, i.e. a value for each timestep, of the global
#'   cover of each state as well as the average local cover of each state. The
#'   full landscape object of each timestep is stored in a list
#'   \code{result$landscapes} of the output object, but frequency of these
#'   snapshots can be altered by increasing the parameter \code{saveeach}.
#'   
#' @export
#' 
#' @examples 
#' 
#' # 1. run simulation and save a snapshot each 50 timesteps. plot timeseries and snapshots.
#' 
#  # create initial landscape
#' l <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = 100) 
#' p <- list(r = 0.4, d = 0.9, delta = 0.01)   # set parameters
#' r <- ca(l, model = musselbed, parms = p, t_max = 200)    # run simulation
#' plot(r)
#' 
#' par(mfrow= c(2,3))
#' sapply(c(0,25,50,100,150,200)+1, function(i) plot(r$landscapes[[i]]) )
#' 
#' # 2. run simulation and save full landsape at each timestep. create animated gif.
#' 
#' l <- init_landscape(c("1","0"), c(0.6,0.4), 100)
#' r <- ca(l, model = life, t_max = 400)
#' animate(r, "life01.gif")


ca_snapsSS <- function(x, model = grazing, parms = "default", 
               t_min = 50,
               t_max = 200, nsnaps = 10,
               stopifsteady = FALSE, 
               steady = caspr::steady, 
               steadyparms = list(t_eval = 200, accept = 0.001),
               plotting = FALSE,
               filename = "modelrun",
               seed = NULL, ... )  {
  
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
  #if(stopifsteady && identical(steadyparms$dim, x$dim)  ) {
  #  warning(paste0("steadyness criterion was defined for lattice of dimensions: ",steadyparms$dim[1], "x",steadyparms$dim[1], " and might be invalid! Please adjust parameter 'steadyparms'!" ))
  #}
  if (t_min>t_max) stop("t_min is greater than t_max of the simulation. Please, define t_min < t_max")
  ## states
  if(!all(x$states %in% model$states)) stop("invalid cell states specified in landscape object!")
  #notworking!  if(!identical(model$states, x$states)) warning("one or several of the model's default cell states are not present in initial landscape object!") 
  ## is interaction matrix provided? if not defaults to four cell neighborhood. 
  if(is.null(model$interact)) {
    I <- matrix(c(0, 1, 0,  1, NA, 1, 0, 1, 0), ncol = 3, byrow = TRUE) 
  } else {
    I <- model$interact 
  }
  if (nsnaps < 3) stop("please don't specify less than 3 snapshots")
  
  mapping(x$dim[1], x$dim[2], i_matrix = I )
  xstats <- summary(x)
  states <- levels(x$cells)
  
  # initialise result object:
  
  result <- list()  # generate list object
  result$model <- model
  result$model$parms <- parms
  
  #result$evaluate <- c(t_min, t_min+2*t_eval)+1
  result$issteady <- logical(length = t_max + 1 )
  #result$steadyval <- numeric( length = t_max + 1 )
  
  result$cover <- as.data.frame(t(xstats$cover))
  result$cover <- result$cover[rep(1, t_max+1),] # preallocating memory
  
  result$local <- as.data.frame(t(xstats$local))
  result$local <- result$local[rep(1, t_max+1),] # preallocating memory
  
  result$seed <- seed       # save seed 
  result$ini_landscape <- x # save initial landscape object at t=0
  
  result$snaps <- 1:nsnaps
  result$landscapes <- list()
  result$landscapes <- lapply(1:nsnaps, function(i) x) # preallocating memory
  
  # --------------- start simulation -----------------
  
  parms_temp <- parms   # set temporary parameter object
  
  # initialise simulation variables: 
  x_old <- x  # ghost matrix at t_i
  x_new <- x_old
  i = 0  # iterator for simulation timesteps
  if(!is.null(seed)) set.seed(seed)  # get seed from function call
  
  # running untill tmin
  while ( i <= t_min ){

        # call update function:
    
    model$update(x_old, parms, ...) -> x_new
    
    # replace ghost matrix for next iteration
    
    x_old <- x_new 
    i = i+1
    
    # save stats of new landscape
    
    xstats <- summary(x_new)
    result$cover[i,] <- xstats$cover
    result$local[i,] <- xstats$local
    
    #result$steadyval[i] <- steady(i, result, steadyparms, returnvalue = TRUE)
    result$issteady[i] <- steady(i, result, steadyparms)
    
  }
  
  # starting iterations:
  while ( i <= t_max && ( !stopifsteady || !steady(i, result, steadyparms) ) ) {
    
    # call update function:
    
    model$update(x_old, parms, ...) -> x_new
    
    # replace ghost matrix for next iteration
    
    x_old <- x_new 
    i = i+1
    
    # save stats of new landscape
    
    xstats <- summary(x_new)
    result$cover[i,] <- xstats$cover
    result$local[i,] <- xstats$local
    
    
    #result$steadyval[i] <- steady(i, result, steadyparms, returnvalue = TRUE)
    result$issteady[i] <- steady(i, result, steadyparms)
    
  }
  
  result$landscapes[[1]] <- x_new
  xcomp<-as.integer(x_new$cells)
  c<-cor(xcomp,as.integer(x_new$cells))
  tsnaps<-0
  
  while (c >= 0.1) {
    
    # call update function:
    
    model$update(x_old, parms, ...) -> x_new
    
    # replace ghost matrix for next iteration
    
    x_old <- x_new 
    i = i+1
    c = cor(xcomp, as.integer(x_new$cells))
    # save stats of new landscape
    
    xstats <- summary(x_new)
    result$cover[i,] <- xstats$cover
    result$local[i,] <- xstats$local
    
    
    #result$steadyval[i] <- steady(i, result, steadyparms, returnvalue = TRUE)
    result$issteady[i] <- steady(i, result, steadyparms)
    tsnaps <- tsnaps + 1
  }
  result$landscapes[[2]]<-x_new
  for (k in 3:nsnaps){
    for (j in 1:tsnaps){
      
      # call update function:
      
      model$update(x_old, parms, ...) -> x_new
      
      # replace ghost matrix for next iteration
      
      x_old <- x_new 
      i = i+1
      
      # save stats of new landscape
      
      xstats <- summary(x_new)
      result$cover[i,] <- xstats$cover
      result$local[i,] <- xstats$local
      
      # save landscape if snapshot
      
      #if(i %in% snaps) {  
      #  result$landscapes[[match(i, snaps)]] <- x_new
      #}
      
      #result$steadyval[i] <- steady(i, result, steadyparms, returnvalue = TRUE)
      result$issteady[i] <- steady(i, result, steadyparms)
      
    }
    result$landscapes[[k]]<-x_new
  }
  
  result$time <- seq(0, i-1) # add vector of realized timesteps
  
  # ------------ end simulation -----------------
  
  #result$steady_state$issteady <- steadiness <= steady
  #result$steady_state$steadiness <- steadiness 
  class(result) <- "ca_result"
  
  if(plotting) ca_animate(result, filename = filename)
  return(result)
}



#' @export

print.ca_result <- function(x) {
  
  obj <- deparse(substitute(x))
  
  cat(" \n ")
  cat("Model run of ", x$model$name, " over ", max(x$time),
      " timesteps. \n", sep = "")
  
  #if(x$steady_state$issteady == TRUE) {
  # cat("\n average cover reached steady state ('steadiness' <= ", 
  #    x$steady_state$steady ,"): \n", sep = "")
  #} else {
  #  cat("\n average cover did not reach steady state ('steadiness' = ", 
  #      round(x$steady_state$steadiness, 6) ,"): \n", sep = "")
  #}
  cat("\n average cover: \n")
  cat( "  ", formatC(names(summary(x)$mean_cover), width = 6, flag = " " ) , "\n")
  cat( "  ", formatC(round(summary(x)$mean_cover, digits = 3), width = 6, flag = " " ) )
  cat(" \n \n")
  cat(" for more details call 'summary(", obj,")' or 'plot(", obj,")'! \n", sep = "")
  
  cat( " access simulation results: \n")
  cat( "   '", obj,"$model$parms' : simulation parameters \n", sep = "")
  cat( "   '", obj,"$time' : a vector of timesteps \n", sep = "")
  cat( "   '", obj,"$cover' : a dataframe of the states' timeseries \n", sep = "")
  cat( "   '", obj,"$ini_landscape' : the initial landscape object \n", sep = "")
  cat( "   '", obj,"$snaps' : an index table of saved snapshots \n", sep = "")
  cat( "   '", obj,"$landscapes[[i]]' : extract snapshot from list of snapshots \n", sep = "")
  cat(" \n ")
  
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
#' @return Returns a list \code{out} containing the model name, the final time,
#'   the mean cover  and standard deviation in the last 10 timesteps.
#'   
#' @note The mean and sd are supposed to depend on the model being in steady
#'   state dynamics, i.e. beyond transitory. Lacking a universal criterion for
#'   that, we report the final 10 timesteps.
#'   
#' @export

summary.ca_result <- function(x) {
  out <- list()
  class(out) <- "ca_summary"
  eval <- tail(seq_along(x$time),10)
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


#' Transfer output of runca() into an array. 
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
