#' @title Grazing model
#'   
#' @description A spatially-explicit model of grazing pressure in drylands.
#'   
#' @usage ca(l, grazing)
#'   
#' @param del local seed dispersal
#' @param b environmental quality
#' @param c_ global competition
#' @param m0 intrinsic mortality
#' @param g grazing pressure
#' @param r regeneration rate of degraded cells
#' @param f local facilitation
#' @param d intrinsic degradation rate
#' @param p associational resistance against grazing
#'   
#' @author Florian D. Schneider and Sonia Kéfi (2015, in revision)
#'   
#' @details The model builds upon a published model by Kéfi et al. 2007. Spatial
#'   models of vegetation cover so far have considered grazing mortality a 
#'   rather constant pressure, affecting all plants equally, regardless of their
#'   position in space. In the known models it usually adds as a constant to the
#'   individual plant risk (Kéfi et al 2007 TPB). However, grazing has a strong 
#'   spatial component: Many plants in rangelands invest in protective 
#'   structures such as thorns or spines, or develop growth forms that reduce 
#'   their vulnerability to grazing. Therefore, plants growing next to each 
#'   other benefit from the protection of their neighbors.
#'   
#'   Such \strong{associational resistance} is widely acknowledged in vegetation
#'   ecology but hardly integrated in models as a cause for spatially 
#'   heterogenous grazing pressure. It also renders the plant mortality density 
#'   dependent, which has important impacts on the bistability of the system.
#'   
#'   The model investigates how the assumption of spatially heterogeneous 
#'   pressure alters the bistability properties and the response of spatial 
#'   indicators of catastrophic shifts.
#'   
#'   The model knows three different cell states: occupied by vegetation 
#'   \code{"+"}, empty but fertile \code{"0"} and degraded \code{"-"}. 
#'   Transitions between cell states are only possible between vegetated and 
#'   empty (by the processes of plant 'death' and 'recolonization') and between 
#'   empty and degraded (by 'degradation' and 'regeneration').
#'   
#'   To account for the spatially heterogeneous impacts of grazing due to
#'   associational resistance, we assumed that a plant's vulnerability to
#'   grazers decreases with the proportion of occupied neighbors,  $q_{+|+}$.
#'   The individual probability of dying is therefore defined as
#'   
#'   \deqn{	w_{ \left\{ +,0 \right\} }  = m_0 + g_0 \left( 1 - q_{+|+} \right)}
#'   
#'   where the additional mortality due to grazing is maximized to \eqn{g_0} if
#'   a plant has no vegetated neighbor (i.e., \eqn{q_{+|+} = 0}) and gradually
#'   reduces to 0 with an increasing fraction of occupied neighbors,
#'   \eqn{q_{+|+}}.

#' 
#' @seealso \href{https://github.com/cascade-wp6/2015_schneider_kefi}{project on GitHub}
#' 
#' @export

"grazing"

grazing <- list()
class(grazing) <- "ca_model"
grazing$name <- "Spatial Grazing Model"
grazing$ref <- "Schneider and Kéfi 2015, in review"
grazing$states <- c("+", "0", "-")
grazing$cols <- grayscale(3)
grazing$parms <- list(
  del = 0.9, # local seed dispersal
  b = 0.2, # environmental quality
  c_ = 0.2, # global competition
  m0 = 0.05, # intrinsic mortality
  g = 0.2, # grazing pressure
  r = 0.01, # regeneration rate of degraded cells
  f = 0.9, # local facilitation
  d = 0.1, # intrinsic degradation rate 
  p = 1 # associational resistance against grazing
) 
grazing$update <- function(x_old, parms_temp, subs = 10, timestep = NA) {
  
  x_new <- x_old
  
  for(s in 1:subs) {
    
    
    parms_temp$rho_plus <- sum(x_old$cells == "+")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    parms_temp$Q_plus <- neighbors(x_old, "+")/4  # count local density of occupied fields for each cell:
    
    # 2 - drawing random numbers
    rnum <- runif(x_old$dim[1]*x_old$dim[2]) # one random number between 0 and 1 for each cell
    
    # 3 - setting transition probabilities
    
    if(parms_temp$rho_plus > 0) {
      recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*1/subs)
      
      death <- with(parms_temp, (m0+g*(1-p*Q_plus))*1/subs)
      death[death > 1] <- 1 
    } else {
      recolonisation <- 0
      death <- 1
    }
    
    regeneration <- with(parms_temp, (r + f*Q_plus)*1/subs)
    
    degradation <- with(parms_temp, (d*1/subs))
    
    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in  time step", timestep, "! decrease number of substeps!!!")) 
    if(any(recolonisation < 0)) warning(paste("recolonisation falls below 0 in time step",timestep, "! balance parameters!!!")) 
    if(any(degradation < 0)) warning(paste("degradation falls below 0 in time step",timestep, "! balance parameters!!!"))
    if(any( death < 0)) warning(paste("death falls below 0 in time step",timestep, "! balance parameters!!!")) 
    if(any(regeneration < 0)) warning(paste("regeneration falls below 0 in time step",timestep, "! balance parameters!!!")) 
    
    # 4 - apply transition probabilities  
    
    x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
    x_new$cells[which(x_old$cells == "+"  & rnum <= death)] <- "0"
    x_new$cells[which(x_old$cells == "0"  & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"  
    x_new$cells[which(x_old$cells == "-"   & rnum <= regeneration)] <- "0"
    
    # 5- store x_new as next x_old
    
    x_old <- x_new
    
  }
  
  ## end of single update call
  return(x_new)
}


