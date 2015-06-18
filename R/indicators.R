
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

indicators <- function(x, spatial = TRUE, temporal = TRUE) {
  require(moments)
  i_out <- list()
  class(i_out) <- "ca_indicators"
  i_out$model <- x$model
  eval <- x$evaluate[1]:x$evaluate[2]
  i_out$mean_cover <- colMeans(x$cover[eval,])
  i_out$sd_cover <- sapply(x$cover[eval,], sd)
  i_out$skewness_cover <- sapply(x$cover[eval,], skewness)
  i_out$mean_local <- colMeans(x$local[eval,])
  i_out$sd_local <- sapply(x$local[eval,], sd)  
  i_out$skewness_local <- sapply(x$local[eval,], skewness)
  i_out$mean_clustering <- colMeans(x$local[eval,]/x$cover[eval,])
  i_out$sd_clustering <- sapply(x$local[eval,]/x$cover[eval,], sd)
  #i_out$autocorrelation
  i_out$patches <- lapply(x$landscapes, patches)
  #i_out$cpd <- data.frame(n = unique())
  i_out$largest <- sapply(i_out$patches, max)
  #i_out$...
  
  return(i_out)
}
