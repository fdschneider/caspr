#' Game of Life
#' 
#' @details Model implementation for Conway's Game of Life.
#' 
#' @examples
#' 
#' l <- init_landscape(c("1","0"), c(0.35,0.65), 150,150)
#' 
#' liferun <- ca(l, life)
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


