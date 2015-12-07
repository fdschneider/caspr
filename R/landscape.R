
#' Create a new landscape object, i.e. a grid for cellular automata processing.
#' 
#' @param states A character vector containing the potential cell states.
#' @param cover A number vector containing the cover.
#' @param width An integer number.
#' @param height An integer number.
#' @return A landscape object of dimensions \code{width} x \code{height} with 
#'   random distribution of \code{states}, in the relative ratio given in 
#'   \code{cover}.
#' @examples 
#' init_landscape(c("+","0","-"), c(0.5,0.25,0.25))
#' @export


init_landscape <- function(states, cover, width = 50, height = width) {
  if(height != width) warning("the landscape is not a square. Spatial analysis will not be possible!")
  if(length(states) != length(cover))
    warning("length of vector 'states' and vector 'cover' differs. I am distributing cover at random!"
  )
  cover <- cover/sum(cover) #normalizing total cover to 1
  cells <- factor(rep(states[1], times = width*height), levels = states) # vector of empty cells in state "0"
  for(i in 2:length(states)) {
    s <- floor(width*height*cover[i]) #how many cells are in state[i]
    cells[sample(which(cells %in% states[1]), s, replace = FALSE)] <- states[i]  # replace s cells, randomly drawn, with state[i]
    }
  #ALT: sample(states, width*height, replace = TRUE, prob = cover )
  # wrap landscape object:
  initial <- list(  
    dim = c(width = as.integer(width), height = as.integer(height)),  # first element contains the dimensions of the landscape 
    cells = cells  #contains a random row-wise, factorial vector to fill the grid 
  )
  levels(initial$cells) <- states  #assign cell states 
  # [Alex] Landscape class should be first so we can override methods for lists 
  # (e.g. as.data.frame, see methods(class='list') ).
  class(initial) <- c("landscape","list") # set class of object (required for plotting)
  return(initial)
}

#' transfer landscape to matrix.
#'
#' @param x A landscape object.
#'
#' @return A matrix. 
#' @export

as.matrix.landscape <- function(x, as = "character") {
  if (as == "character") {
    matrix(x$cells, nrow = x$dim[1], byrow = TRUE)
  } else if (as == 'integer') {
      matrix(as.integer(x$cells), nrow = x$dim[1], byrow = TRUE)
# Cannot make a factorial matrix ?!
#   } else if (as =="factor") { 
#     matrix(x$cells, nrow = x$dim[1], byrow = TRUE)
  } else { 
    stop(paste0('Cannot convert landscape object to class ', as, "."))
  }
}

#' @export
as.binary_matrix <-  function (...) UseMethod("as.binary_matrix")

#' @export
as.binary_matrix.landscape <- function(x, is = levels(x$cells)[1]) { 
  # The use of %in% does not keep the matrix format => use a == instead (/!\ 
  #   this means we cannot consider two states as 1).
  if ( length(is) > 1 ) stop('Please specify only one state')
  tmp <- as.matrix(x) == is
  class(tmp) <- c("binary_matrix", "matrix")
  tmp
}
#' @export
as.binary_matrix.ca_result <- function(x, is = x$model$states[1]) lapply(x$landscapes, function(y) as.matrix(y) %in% is ) 



#' generic S3 function
#' @export

as.landscape <- function (...) UseMethod("as.landscape")

#' transfer matrix to landscape.
#' 
#' @param x A matrix object with factorial content.
#' @param levels A character vector giving the original sequence of states.
#'   
#' @return A landscape object.
#'   
#' @details If no vector levels is specified the levels are coerced from the 
#'   unique states of the matrix, which might cause loss of states (e.g. if not 
#'   present on this landscape) or a wrong ordering of states. This is not
#'   recommended! Always provide a levels vector!
#'   
#' @export

as.landscape.matrix <- function(x, states = NA)  {
  
  # Assumes the matrix x is not a factor, which seems like an impossible 
  # situation [Alex]
  if (is.na(states)) { 
    states <- unique(as.vector(x))
    warning('Importing states from the unique elements of the matrix.',
            'This is not recommended !')
  }
  
  structure(
    list(
      dim = c(width = dim(x)[2] , height = dim(x)[1]), 
      cells = factor(matrix(t(x), nrow = 1 ), levels = states )
    ),
    class = "landscape"
    )
  
}

#' Summary of landscape object.
#' 
#' @param x A landscape object
#' @return A list containing the object dimensions, the number of cells in the different cell states, the global cover of the different cell states, the average local cover of the different cell states.
#' @examples 
#' # obj <- init_landscape(c("+","0","-"), c(0.5,0.25,0.25)) 
#' # summary(obj)
#' @export


summary.landscape <- function(x) {
  if(!exists("x_with_border")) mapping(x$dim[1],x$dim[2], i_matrix = I)
  out <- list()
  out$dim <- x$dim
  out$n <- sapply(levels(x$cells), function(y) {sum(x$cells == y)})
  out$cover <- out$n/sum(out$n)
  out$local <- sapply(levels(x$cells), 
                      function(y) {mean(  (neighbors(x,y)/4)[x$cells == y]  )})
  return(out)
}


#' Plotting an objects of class "landscape"
#' 
#' @param x A landscape object.
#' @param cols A color vector. If \code{"auto"}, then a grayscale vector will be
#'   used.
#' @param grid If TRUE, plot a grid over the cells. Defaults to FALSE.
#' @param axis If TRUE, plot x and y axis with coordinates.  Defaults to FALSE.
#' @param add If TRUE, no primary plot will be called. The change in cells will
#'   be plotted over an existing plot. Used for animated plotting in the screen
#'   plotting device.
#' @param ani If TRUE, an adjustment constant is added when producing a pixel
#'   accurate png or gif file. Required when the function is used to plot
#'   animated figures.
#' @return A landscape object of dimensions \code{width} x \code{height} with 
#'   random distribution of \code{states}, in the relative ratio given in 
#'   \code{cover}.
#' 
#' 
#' @export
#' @examples 
#' obj <- init_landscape(c("+","0","-"), c(0.5,0.25,0.25)) 
#' plot(obj)
#' 

plot.landscape <- function(x, cols = "auto", grid = FALSE, axis = FALSE, 
                           add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") { 
    cols = grayscale(nlev)  # default color value
    warning('No colors provided for landscape plotting: choosing ad-hoc ', 
            'grey levels')
  }
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, x$dim[1]+0.5+adj), ylim = c( x$dim[2]+0.5+adj, 0+0.5+adj), bty = c("n", "o")[grid+1], xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1, ... ) 
  
  if(axis && !add) axis(3) 
  if(axis && !add) axis(2)
  
  if(grid) border = "grey80" else border = cols[as.numeric(x$cells)]
  
  rect(rep(1:x$dim[1], times = x$dim[2])-.5, rep(1:x$dim[2], each = x$dim[1])-.5, rep(1:x$dim[1], times = x$dim[2])+.5, rep(1:x$dim[2], each = x$dim[1])+.5, col = cols[as.numeric(x$cells)], border = border)
  
  if(grid) box()
}



#' @export
#' 
print.landscape <- function(x) {
  
  print(list(dim = x$dim, cells = head(x$cells, 25) )) 
  
}

