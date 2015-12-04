#' @title Predator-prey model
#'   
#' @description A spatially-explicit predator prey model.
#'   
#' @usage ca(l, predprey, parms = p)
#'   
#' @param betaf A numerical value. The reproduction rate of prey.
#' @param betas A numerical value. reproduction rate of predators.
#' @param delta A numerical value. Probability of starvation for predators that
#'   do not find prey.
#'   
#' @author Mercedes Pascual, M. Roy, F. Guichard & G. Flierl
#'   
#' @references Pascual, M., Roy, M., Guichard, F., & Flierl, G. (2002). Cluster 
#'   size distributions: signatures of self-organization in spatial ecologies. 
#'   Philosophical Transactions of the Royal Society of London. Series B, 
#'   Biological Sciences, 357(1421), 657–666. 
#'   doi: \href{https://doi.org/10.1098/rstb.2001.0983}{10.1098/rstb.2001.0983}
#' @family models
#' @details The particular neighbourhood we consider in the simulations consists
#'   of the four nearest sites. Prey growth occurs as a contact process: a prey 
#'   chooses a neighbouring site at random and gives birth onto it only if this 
#'   site is empty at rate \code{betas}. Predators hunt for prey by inspecting 
#'   their neighbourhood for the presence of prey at rate 1. If prey are 
#'   present, the predator selects one at random and eats it, moving to this 
#'   neighbouring site. Only predators that find a prey can reproduce, and do so
#'   with a specified probability, \code{betas}. The offspring is placed in the 
#'   original site of the predator. Predators that do not find prey are 
#'   susceptible to starvation and die with a probability \code{delta}. Random 
#'   movement occurs through mixing: neighbouring sites exchange state at a 
#'   constant rate \code{nu}.
#'   
#'   
#' @useDynLib caspr
#' @importFrom Rcpp sourceCpp
#' @export

"predprey"

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
predprey$cols <- c("gray50", "black", "white" )
predprey$parms  <- list(
  betaf = 0.01,  # growth rate of fish
  betas = 0.1,   # growth rate of predator
  delta = 0.2,   # death rate of predator (if starved)
  m     = 0.1    # death rate of prey 
  # The original publication refers to a nu parameter that fixes the mixing 
  #   rate (?). However, it is not used in the c code provided but maybe that is 
  #   because publication always show results where it is set to 1.
#   nu  = 1
)

predprey$update <- function(x,               # landscape object
                            parms,           # set of parms 
                            subs = 10) {     # temporal resolution
  
  # Adjust data (makes copy)
#   x_old <- matrix(as.integer(x$cells), nrow = x$dim[1], byrow = TRUE)
  x_old <- as.matrix(x, as = "integer")

  x_new <- predprey_core(x_old, 
                         1,  # DO NOT USE THE SUBS PARAMETER !
                         parms$betaf, 
                         parms$betas,
                         parms$delta,
                         parms$m)
  
  x$cells <- factor(predprey$states[x_new], levels = predprey$states)

  ## end of single update call
  return(x)
}

# Internal function/not exported
.to_vector <- function(xmat, levels) { 
  factor(levels[t(xmat)], levels = levels)
}
