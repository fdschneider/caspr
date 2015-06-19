# core functions
# --------------

#' Create vector maps for cellular automata updating
#'
#' @param width An integer number. The width of the lattice. 
#' @param height An integer number. The height of the lattice. 
#' @param boundary A character value representing the type of boundary conditions. Defaults to \code{"periodic"}.  
#' @param i_matrix A matrix defining the interaction matrix to be evaluated. Defaults to 4-cell neighborhood. 
#'
#' @return Creates a map of the landscape object to translate it into a vector with boundaries (\code{x_with_border}) and another one to back-translate it to a vector without boundaries  (\code{x_to_evaluate}) into the global environment. This will automatically be called within \code{ca()}. 
#' @export
#'
#' 

mapping <- function(width, height, boundary = "periodic", i_matrix = matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)) {
  
  # derive helper vectors for counting: 
  # transformation vector for evaluation at the border of the grid
  # set evaluation matrix 
  X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
  # setting the border of the evaluation matrix X
  X <- cbind(X[,width], X, X[,1] )  
  X <- rbind(X[height,], X, X[1,] ) 
  # transformation vector which adds the border to the lattice:
  assign("x_with_border", as.integer(t(X))  , envir = .GlobalEnv )
  
  # from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
  #x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )    
  
  
  assign("x_to_evaluate", sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )	, envir = .GlobalEnv )
  
  
  # defining the neighborhood which is to be evaluated	
  # set interaction matrix
  I <- i_matrix	
  # coordinates of neighbours in Interaction matrix I: 
  neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.ind = TRUE)
  # coordinates relative to the evaluated cell (=  which(is.na(I) ) 
  relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
  relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
  
  # relative position of the four direct neighbours of a cell
  #interact <- (relrow * dim(X)[2] + relcol)
  
  assign("interact", relrow * dim(X)[2] + relcol, envir = .GlobalEnv )
  
}



####################
## count function ##
####################

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
#' @return Returns a vector of the counts in the neighborhood specified by the \link{mapping}, by default the 4-cell neighborhood.
#' 

count.landscape <- function(x, neighbor) {
  
  neighbors <- numeric(length = prod(x$dim))
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}



# helper functions
# ----------------

grayscale <- colorRampPalette(c("black", "white"), space = "rgb")
