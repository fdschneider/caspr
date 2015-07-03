# This is deprecated and should not be used any more. 

#' @export

count <- function(...) {
  UseMethod("count")
}

count.integer <- function(y, x, neighbor) {
  neighbors <- numeric(length = length(interact))
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  
  neighbors <- neighbors + x_logical_with_border[x_to_evaluate[y]+interact]
  
  return(neighbors) 
}

#' Count neighbors.
#'
#' @param x A landscape object. 
#' @param neighbor A character value. The state to count. 
#'
#' @return Returns a vector of the counts in the neighborhood specified by the \code{\link{mapping}}, by default the 4-cell neighborhood.
#' 
#' @export

count.landscape <- function(x, neighbor) {
  
  neighbors <- numeric(length = prod(x$dim))
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}
