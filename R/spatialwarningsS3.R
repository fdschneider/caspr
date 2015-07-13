#' @title Wrapper functions for package 'spatialwarnings' 
#' 
#' @description A set of wrapper functions that apply the spatial analysis tools provided by 'spatialwarnings' to caspr landscape objects. 
#' 
#' @usage  
#' 
#' indicators(l)
#' indicators_fitpsd(l)
#' patchsizes(l)
#' 
#' @import spatialwarnings
#' @export

patchsizes <- function(...) UseMethod("patchsizes")

#' @export
#' 
patchsizes.matrix <- function(mat)  spatialwarnings::patchsizes(mat)

#' @export
#' 
patchsizes.landscape <- function(l, state = levels(l$cells)[1]) {
  binary <- as.matrix(l) == state
  spatialwarnings::patchsizes(binary)
}


#' @export
patchsizes.list <- function(l, state = levels(l[[1]]$cells)[1]) {
  binary <- lapply(l, function(x) as.matrix(x) == state)
  spatialwarnings::patchsizes(binary)
}


#' @export
#' 
patchsizes.ca_result <- function(x, state = x$model$states[1]) {
  patchsizes(x$landscape)
}

#' @export
indicator_cumpsd <- function(...) UseMethod("indicator_cumpsd")

#' @export
#' 
indicator_cumpsd.matrix <- function(mat) spatialwarnings::indicator_cumpsd(mat)


#' @export
#' 
indicator_cumpsd.landscape <- function(l, state = levels(l$cells)[1]) {
  binary <- as.matrix(l) == state
  spatialwarnings::indicator_cumpsd(binary)
}

#' @export
indicator_cumpsd.list <- function(l, state = levels(l[[1]]$cells)[1]) {
  binary <- lapply(l, function(x) as.matrix(x) == state)
  spatialwarnings::indicator_cumpsd(binary)
}
