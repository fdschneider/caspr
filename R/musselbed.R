#' @title Musselbed Disturbance Model 
#' @description Model for spatial pattern of intertidal mussel beds with wave 
#'   disturbance. With three cell states: occupied by mussels (\code{"+"}), 
#'   empty sites (\code{"0"} and disturbed sites (\code{"-"}).
#'   
#' @usage 
#' ca(l, musselbed, parms = p)
#' 
#' @param r A numerical value. recolonisation of empty sites dependent on local 
#'   density.
#' 
#' @param d A numerical value. probability of disturbance of occupied sites if 
#'   at least one disturbed site is in the direct 4-cell neighborhood. 
#' 
#' @param delta A numerical value. intrinsic disturbance rate.
#' 
#' @author Guichard, Halpin, et al. (2003)
#' 
#' @references Guichard, F., Halpin, P.M., Allison, G.W., Lubchenco, J. & Menge,
#'   B.A. (2003). Mussel disturbance dynamics: signatures of oceanographic 
#'   forcing from local interactions. The American Naturalist, 161, 889–904.
#' 
#' @family models
#' 
#' @details 
#' 
#' The model represents the spatial dynamics in mussel cover of rock substrate 
#'   in intertidal systems. The stochastic wave disturbances will most likely 
#'   remove mussels that are located next to a gap because of the losened byssal
#'   threads in their proximity. This causes a dynamic gap growth.
#' 
#' The model describes the process by simplifying the system into three 
#'   potential cell states: occupied by mussel ("+"),  empty but undisturbed 
#'   ("0"), and disturbed, bare rock with loose byssal threads ("-"). 
#'  
#' Mussel growth on empty cells is defined by parameter \code{r} multiplied by 
#'   the local density of mussels in the direct 4-cell neighborhood.  
#'  
#' Any cell occupied by mussels has an intrinsic chance of \code{delta} to 
#'   become disturbed from intrinsic cause, e.g. natural death or predation. 
#'   Additionally, wave disturbance will remove mussels and leave only bare 
#'   rock, i.e. disturbed sites, with probability \code{d} if at least one 
#'   disturbed cell is in the direct 4-cell neighborhood. This causes 
#'   disturbances to cascade through colonies of mussels.  
#'  
#' Disturbed sites will recover into empty sites with a constant rate of 1 per 
#'   year, i.e. on average a disturbed site becomes recolonisable within one 
#'   year after the disturbance happened. 
#'
#' @examples 
#' 
#' l <- init_landscape(c("+","0","-"), c(0.6,0.2,0.2), width = 50) # create initial landscape
#' p <- list(delta = 0.01, d = 09, r = 0.4)   # set parameters
#' r <- ca(l, musselbed, p, t_max = 100)    # run simulation 
#' 
#' 
#' @export
#' 


"musselbed"

musselbed <- list()
class(musselbed) <- "ca_model"
musselbed$name <- "Musselbed Disturbance Model"
musselbed$ref <- "Guichard, F., Halpin, P.M., Allison, G.W., Lubchenco, J. & Menge, B.A. (2003). Mussel disturbance dynamics: signatures of oceanographic forcing from local interactions. The American Naturalist, 161, 889–904. doi: 10.1086/375300"
musselbed$states <- c("+", "0", "-")
musselbed$cols <- grayscale(3)
musselbed$parms <- list(
  r = 0.4, # recolonisation of empty sites dependent on local density
  d = 0.9, # probability of disturbance of occupied sites if at least one disturbed site
  delta = 0.01 # intrinsic disturbance rate
)
musselbed$interact<-matrix(c(1,1,1,1,NA,1,1,1,1), ncol = 3, byrow = 3)
musselbed$update <- function(x_old, parms_temp, delta = 0.2, subs = 10, timestep = NA) {
  
  x_new <- x_old
  
  for(s in 1:subs) {
    
    parms_temp$rho_plus <- sum(x_old$cells == "+")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    parms_temp$localdisturbance <- neighbors(x_old, "-") > 0  # any disturbance in neighborhood?
    parms_temp$localcover <- neighbors(x_old, "+")/8 # any occupied in neighborhood? 
    
    # 2 - drawing random numbers
    rnum <- runif(x_old$dim[1]*x_old$dim[2]) # one random number between 0 and 1 for each cell
    
    # 3 - setting transition probabilities
    
    if(parms_temp$rho_plus > 0) {
      recolonisation <- with(parms_temp, (r*localcover)*1/subs)
      
      disturbance <- with(parms_temp, (delta+d*localdisturbance)*1/subs)
      disturbance[disturbance > 1] <- 1 
    } else {
      recolonisation <- 0
      disturbance <- 1
    }
    
    regeneration <- 1*1/subs
    
    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(recolonisation, disturbance, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in  time step", timestep, "! decrease number of substeps!!!")) 
    if(any(recolonisation < 0)) warning(paste("recolonisation falls below 0 in time step",timestep, "! balance parameters!!!")) 
    if(any( disturbance < 0)) warning(paste("disturbance falls below 0 in time step",timestep, "! balance parameters!!!")) 
    if(any(regeneration < 0)) warning(paste("regeneration falls below 0 in time step",timestep, "! balance parameters!!!")) 
    
    # 4 - apply transition probabilities  
    
    x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
    x_new$cells[which(x_old$cells == "+"  & rnum <= disturbance)] <- "-"
    x_new$cells[which(x_old$cells == "-"   & rnum <= regeneration)] <- "0"
    
    # 5- store x_new as next x_old
    
    x_old <- x_new
    
  }
  
  ## end of single update call
  return(x_new)
}
