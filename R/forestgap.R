#' Forest gap model
#'
#' @details We implement the version used by Kubo et al. (1996).  More details on the implementation are available at: \link{https://github.com/skefi/SpatialStress/wiki/forestgap}
#' 
#' @export

"forestgap"

forestgap <- list()
class(forestgap) <- "ca_model"
forestgap$name   <- "Forest Gap Model"
forestgap$ref    <- paste0("Kubo, T., Y. Iwasa and N. Furumoto. 1996. Forest ",
                           "spatial dynamics with gap expansion: Total gap ",
                           "area and gap size distribution. Journal of ",
                           "Theoretical Biology 180:229–246.")
forestgap$states <- c("+", "0")
forestgap$cols   <- grayscale(2)
forestgap$parms  <- list(
  # Default values are taken from the original publication and produce 
  # bistability (untested) for \delta within 0.15-0.2.
  alpha = 0.20, # recovery strength 
  delta = 0.17, # strength of gap expansion
  d = 0.01      # intrinsic death rate
)

forestgap$update <- function(x,               # landscape object
                             parms,           # set of parms 
                             subs = 10) {     # number of iterations/time step ?

  for (s in 1:subs) {
    
    # Get the local proportion of empty neighbors 
    localempty <- count(x, "0") / length(interact)
    rho_plus <- sum(x$cells == "+") / prod(x$dim)
    
    # 2 - drawing random numbers
    # one random number between 0 and 1 for each cell
    rnum <- runif(prod(x$dim), 0, 1)
    
    # 3 - set transition probabilities
    p_to_vegetated <- with(parms, alpha * rho_plus)
    p_to_empty     <- with(parms, d + delta * localempty)
    
    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(p_to_empty, p_to_vegetated) > 1 )) {
      warning("a set probability is exceeding 1, please balance parameters!")
    }
    
    # 4 - apply transition probabilities
    x_new <- x
    x_new$cells[x$cells == "+" & rnum <= p_to_empty]     <- "0"
    x_new$cells[x$cells == "0" & rnum <= p_to_vegetated] <- "+"
    
    # 5- store x_new as next x
    x <- x_new
    
  }
  
  ## end of single update call
  return(x)
}
