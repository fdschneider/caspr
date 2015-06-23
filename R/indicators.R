
#' Generic function for calculating indicators.
#'
#' @param ... 
#' @details with methods for objects of class 'landscape' and 'ca_result'.
#' 
#' @export 

indicators <- function(...) UseMethod("indicators")

#' Get spatial indicators for landscape object.
#'
#' @param x A landscape object.
#' @param state The cell state of interest. Defaults to \code{x$state[1]}.
#' @param ... 
#'
#' @return a list of spatial indicators. 
#' @export
#'

indicators.landscape <- function(x, state = levels(x$cells)[1], ...) {
  
  x_mat <- as.matrix(x)
  
  out <- spatialwarnings::spatial_ews(x_mat==state, ...)
  return(out)
}

#' Get spatial and temporal indicators for a ca_result object
#' 
#' @param x A \code{ca_result} object, output of function \code{ca()}. 
#' @param spatial Logical. If true, function returns all spatial indicators.
#' @param temporal Logical. If true, function returns all temporal indicators.
#'   
#' @return An object of class \code{ca_indicators} that contains a detailed
#'   report on all kinds of spatial and temporal indicators on the timeseries
#'   and snapshots.
#'   
#' @export

indicators.ca_result <- function(x, spatial = TRUE, temporal = TRUE) {
  require(moments)
  i_out <- list()
  class(i_out) <- "ca_indicators"
  i_out$model <- x$model
  eval <- x$evaluate[1]:x$evaluate[2]
  i_out$mean_cover <- colMeans(x$cover[eval,])
  i_out$sd_cover <- sapply(x$cover[eval,], sd)
  i_out$skewness_cover <- sapply(x$cover[eval,], moments::skewness)
  i_out$mean_local <- colMeans(x$local[eval,])
  i_out$sd_local <- sapply(x$local[eval,], sd)  
  i_out$skewness_local <- sapply(x$local[eval,], moments::skewness)
  i_out$mean_clustering <- colMeans(x$local[eval,]/x$cover[eval,])
  i_out$sd_clustering <- sapply(x$local[eval,]/x$cover[eval,], sd)
  #i_out$autocorrelation
  i_out$patches <- lapply(x$landscapes, function(i) SDMTools::PatchStat(as.matrix(label(i)))  ) 
  #i_out$cpd <- data.frame(n = unique())
  i_out$largest <- sapply(i_out$patches, max)
  
  #put all landscapes into one long matrix to be handled by spatial_ews()
    x_all <- lapply(1:length(x$landscapes), function(i) x$landscapes[[i]]$cells )
    x_all <- unlist(x_all)
    x_all <- matrix(x_all, ncol = x$landscapes[[1]]$dim[1], byrow = TRUE)
    x_all <- as.matrix(x_all)
    
  i_out$spatialstats <- lapply(x$landscapes, function(i) spatialwarnings::spatial_ews(x_all)) 
  
  return(i_out)
}
