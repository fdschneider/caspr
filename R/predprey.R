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
predprey$ref    <- paste0("Pascual, M., Roy, M., Guichard, F., & Flierl, G. ",
                          "(2002). Cluster size distributions: signatures of ",
                          "self-organization in spatial ecologies. ",
                          "Philosophical Transactions of the Royal Society of ",
                          "London. Series B, Biological Sciences, 357(1421), ",
                          "657–666. http://doi.org/10.1098/rstb.2001.0983")
predprey$states <- c("f", "s", "0") # fish/shark/empty /!\ ORDER MATTERS
predprey$parms  <- list(
  betaf = 1/3,  # growth rate of preypredprey
  betas = 1/10, # growth rate of predator
  delta = 1/3 # death rate of predator (if starved)
  # The original publication refers to a nu parameter that fixes the mixing 
  #   rate (?). However, it is not used in the c code provided but maybe that is 
  #   because publication always show results where it is set to 1.
#   nu  = 1
)

predprey$update <- function(x,               # landscape object
                            parms,           # set of parms 
                            subs = 1000) {     # number of iterations/time step ?
  
  # Adjust data (makes copy)
  x_old <- matrix(as.integer(x$cells), nrow=x$dim[1], byrow=TRUE)
  
  x_new <- predprey_core(x_old, 
                         subs, 
                         length(interact), # global variable
                         parms$betaf, 
                         parms$betas,
                         parms$delta)
  
  x$cells <- as.factor(predprey$states[as.vector(t(x_new))])
  ## end of single update call
  return(x)
}

