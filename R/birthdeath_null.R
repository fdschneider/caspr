#' @title Birth-death model
#'
#' @description Model with two states and a probability of switching between the
#'   two. 
#'
#' @usage 
#' ca(l, birthdeath, parms = p)
#' 
#' @author None
#' @references None
#' @details 
#' This model implements only two states (dead or alive). Each cell can 
#'   transition to the other depending only on a probability transition p (
#'   dead -> alive) or q (alive -> dead). This model serves as a null model 
#'   for the others provided for caspr. 
#'
#' @export
#' 
"birthdeath"

birthdeath <- list()
class(birthdeath) <- "ca_model"  

birthdeath$name <- "Random birth-death model"

birthdeath$ref <- NA  # a bibliographic reference

birthdeath$states <- c("+", "0")

birthdeath$cols <- c("black", "white")

# a list of default model parameters, used to validate input parameters. 
birthdeath$parms <- list(  
  death = 0.5,
  birth = 0.5
) 

# an update function
birthdeath$update <- function(x_old, parms) {
  
  x_new <- x_old
  
  # define update procedures depending on parms 
  
  rnum <- runif(length(x_old$cells))
  x_new$cells[x_old$cells == "+" & rnum <= parms$death] <- "0"
  x_new$cells[x_old$cells == "0" & rnum <= parms$birth] <- "+"
  
  return(x_new)
}


