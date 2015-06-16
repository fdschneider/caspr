########################################################
# The MIT License (MIT)
#
# Copyright (c) 2014 Florian D. Schneider & Sonia Kéfi
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
########################################################

### mussel bed model (Guichard et al 2003)


musselbed <- list()
class(musselbed) <- "ca_model"
musselbed$name <- "Mussel Disturbance Model"
musselbed$ref <- "Guichard et al 2003, American Naturalist, Vol. 161, pp. 889–904"
musselbed$states <- c("+", "0", "-")
musselbed$cols <- grayscale(3)
musselbed$parms <- list(
  r = 0.4, # recolonisation of empty sites dependent on local density
  d = 0.9, # probability of disturbance of occupied sites if at least one disturbed site
  delta = 0.01 # intrinsic disturbance rate
) 
musselbed$update <- function(x_old, parms_temp, delta = 0.2, subs = 10, timestep = NA) {
  
  x_new <- x_old
  
  for(s in 1:subs) {
    
    parms_temp$rho_plus <- sum(x_old$cells == "+")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    parms_temp$localdisturbance <- count(x_old, "-") > 0  # any disturbance in neighborhood?
    parms_temp$localcover <- count(x_old, "+")/4 # any occupied in neighborhood? 
    
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
