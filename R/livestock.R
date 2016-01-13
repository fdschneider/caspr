#' @title Livestock resilience model
#' 
#' @description A spatially-explicit model of resilience in arid rangelands.
#' @usage 
#' ca(l, livestock, parms = p)
#' 
#' @param b environmental quality
#' @param m intrinsic plant mortality, i.e. inverse of av. lifespan
#' @param r max. regeneration rate of plants
#' @param K carrying capacity of the landscape (i.e. global competition)
#' @param f local facilitation
#' @param c local competition
#' @param L livestock density
#' @param a attack rate of livestock
#' @param h handling time of livestock
#' @param q hill-exponent of livestock
#' @param p associational resistance against livestock grazing
#' @param v attractant-decoy effect of plants to livestock
#' 
#' @author F.D. Schneider
#' @family models
#' @details 
#' An unpublished, generalized model of positive and negative local feedbacks in arid rangelands, including mechanisms such as local facilitation, competition, associational resistance and attractant-decoy. 
#' 
#' @export
"livestock"

livestock <- list()
class(livestock) <- "ca_model"  
livestock$name <- "Livestock resilience model"
livestock$ref <- NA  # a bibliographic reference
livestock$states <- c("1", "0")
livestock$cols <- c("black", "white")
# a list of default model parameters, used to validate input parameters. 
livestock$parms <- list(  
  r = 1.0,  # max. regeneration rate of plants
  b = 0.5,  # environmental quality
  sigma = 0, #
  f = 0.8,  # local facilitation
  alpha = 0, # water runoff
  K = 0.9, # carrying capacity of the system
  c = 0.2, # local competition
  m = 0.05, # intrinsic mortality of plants (inverse of av. lifespan)
  v = 0.8, # attractant-decoy
  p = 0.99, # associational resistance
  L = 5, # Livestock density
  q = 0, # hill exponent of functional response
  h = 50, # handling time 
  a = 0.2 # attack rate of livestock
) 
# an update function
livestock$update <- function(x_old, parms, subs = 12) {
  
  x_new <- x_old
  
  climate <- parms$b 
  if(parms$sigma != 0)  climate <- rnorm(1, climate, parms$sigma)
  if(climate < 0) climate <- 0
  
  for(s in 1:subs) {
    
    # define update procedures depending on parms 
    
    # model specific part:
    # 1 - setting time-step parameters
    rho_one <- sum(x_old$cells == "1")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    q_one_one <- neighbors(x_old, "1")/4  # count local density of occupied fields for each cell
    
    # 2 - drawing random numbers
    rnum <- runif(x_old$dim[1]*x_old$dim[2]) # one random number between 0 and 1 for each cell
    
    # 3 - setting transition probabilities
    growth <- with(parms, (r * (climate + (1-climate)*f*q_one_one) * rho_one^(1 + alpha) * ( 1 - (rho_one / (K * (1-c*q_one_one) ))) / (1 - rho_one))  *1/subs)  # recolonisation rates of all cells 
    
    growth[growth < 0] <- 0
    
    death <- with(parms,       (m + ( (a+ v*q_one_one) * (1-p*q_one_one) * L * rho_one^(1+q) )/( 1 + (a+ v*q_one_one) * (1-p*q_one_one) * h * rho_one^(1+q) )) *1/subs)   # set probability of death for each cell
    
    death[death < 0] <- 0
    
    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(growth, death) > 1 )) warning(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
    #if(any(c(growth, death) < 0)) warning(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
    
    # 4 - apply transition probabilities  
    
    x_new$cells[which(x_old$cells == "0" & rnum <= growth)] <- "1"
    x_new$cells[which(x_old$cells == "1" & rnum <= death)] <- "0"
    
    # 5- store x_new as next x_old
    
    x_old <- x_new
    
  }
  
  return(x_new)
  
}


