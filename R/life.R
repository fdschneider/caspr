#' @title Conway's Game of Life
#'   
#' @description Model implementation for Conway's Game of Life.
#' @author John H. Conway
#' @usage 
#' l <- init_landscape(c("1","0"), c(0.35,0.65), 150,150)
#' run <- ca(l, life)
#' @references Gardner, Martin (October 1970). Mathematical Games - The
#'   fantastic combinations of John Conway's new solitaire game "life".
#'   Scientific American 223. pp. 120–123. ISBN 0-89454-001-7.
#' @family models
#' @details  
#' The Game of Life, also known simply as Life, is a cellular automaton devised
#' by the British mathematician John Horton Conway in 1970.
#' 
#' The "game" is a zero-player game, meaning that its evolution is determined by
#' its initial state, requiring no further input. One interacts with the Game of
#' Life by creating an initial configuration and observing how it evolves or,
#' for advanced players, by creating patterns with particular properties.
#' 
#' \strong{Rules}
#' 
#' The universe of the Game of Life is an infinite two-dimensional orthogonal
#' grid of square cells, each of which is in one of two possible states, alive
#' or dead. Every cell interacts with its eight neighbours, which are the cells
#' that are horizontally, vertically, or diagonally adjacent. At each step in
#' time, the following transitions occur:
#' 
#' \itemize{
#'   \item Any live cell with fewer than two live neighbours dies, as if caused 
#'     by under-population.
#'   \item Any live cell with two or three live neighbours lives on to the next 
#'     generation.
#'   \item Any live cell with more than three live neighbours dies, as if by 
#'     overcrowding.
#'   \item Any dead cell with exactly three live neighbours becomes a live cell, 
#'     as if by reproduction.
#' }
#' 
#' The initial pattern constitutes the seed of the system. The first generation
#' is created by applying the above rules simultaneously to every cell in the
#' seed - births and deaths occur simultaneously, and the discrete moment at which
#' this happens is sometimes called a tick (in other words, each generation is a
#' pure function of the preceding one). The rules continue to be applied
#' repeatedly to create further generations.
#' 
#' 
#' Source: \href{https://en.wikipedia.org/wiki/Conway's_Game_of_Life}{Wikipedia}
#' 
#' 
#' @export

"life"

life <- list()
class(life) <- "ca_model"  
life$name <- "Random model"
life$ref <- NA  # a bibliographic reference
life$states <- c("1", "0")
life$cols <- c("black", "white")
# a list of default model parameters, used to validate input parameters. 
life$interact <- matrix(c(1,2,3, 4,NA,6, 7,8,9), ncol = 3, byrow = TRUE) 
life$parms <- list(  
  B = 3,
  S = c(2,3)
) 
# an update function
life$update <- function(x_old, parms) {
  
  x_new <- x_old
  
  neighb <- neighbors(x_old, "1")
  # define update procedures depending on parms 
  
  x_new$cells[x_old$cells == "0" & neighb %in% parms$B] <- "1"
  x_new$cells[x_old$cells == "1" & !neighb %in% parms$S] <- "0"
  
  return(x_new)
}


