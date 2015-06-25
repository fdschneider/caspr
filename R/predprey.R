# Predator-prey model
# 
# 
# These tags are necessary for Rcpp.
# 
#' @useDynLib caspr
#' @importFrom Rcpp sourceCpp
# 

predprey <- list()
class(predprey) <- "ca_model"
predprey$name   <- "Predator-prey Gap Model"
predprey$ref    <- paste0("todo")
predprey$states <- c("f", "s", "0") # fish/shark/empty /!\ ORDER MATTERS
predprey$parms  <- list(
  # Default values are taken from the original publication and produce 
  # bistability (untested) for \delta within 0.15-0.2.
  betaf = 1/3,  # growth rate of prey
  betas = 1/10, # growth rate of predator
  delta = 1/3 # death rate of predator (if starved)
  # The original publication refers to a nu parameter that fixes the mixing 
  #   rate (?). However, it is not used in the c code provided but maybe that is 
  #   because publication always show results where it is set to 1.
#   nu  = 1
)

predprey$update <- function(x,               # landscape object
                            parms,           # set of parms 
                            subs = 10) {     # number of iterations/time step ?
  
  # Adjust data (makes copy)
  x_old <- matrix(as.integer(x$cells), nrow=x$dim[1], byrow=TRUE)
  
  x_new <- predprey_core(x_old, 
                         subs, 
                         length(interact), # number of neighbors
                         parms$betaf, 
                         parms$betas,
                         parms$delta)
  
  x$cells <- as.factor(predprey$states[as.vector(t(x_new))])
  
  ## end of single update call
  return(x)
}
