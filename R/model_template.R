#' @title Template for model objects
#' @description One sentence description. Potential cell states.
#' @usage 
#' ca(l, template, parms = p)
#' @author Author of original model
#' @references Original reference and further references
#' @details 
#' Wordy description about the model details, transition probabilities, biological assumptions.
#' @export
#' 


"template"

template <- list()
class(template) <- "ca_model"  
template$name <- "Random model"
template$ref <- NA  # a bibliographic reference
template$states <- c("1", "0")
template$cols <- c("black", "grey50")
# a list of default model parameters, used to validate input parameters. 
template$parms <- list(  
  death = 0.5,
  birth = 0.5
) 
# an update function
template$update <- function(x_old, parms) {
  
  x_new <- x_old
  
  # define update procedures depending on parms 
  
  rnum <- runif(length(x_old$cells))
  x_new$cells[x_old$cells == "1" & rnum <= parms$death] <- "0"
  x_new$cells[x_old$cells == "0" & rnum <= parms$birth] <- "1"
  
  return(x_new)
}


