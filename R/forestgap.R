########################################################
# The MIT License (MIT)
#
# Copyright (c) 2014 Florian D. Schneider, Sonia Kéfi & Alexandre Génin
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

# Forest gap model
# 
# We implement the version used by Kubo et al. (1996).
# More details on the implementation are available at: 
# 
# https://github.com/skefi/SpatialStress/wiki/forestgap
# 

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
                             subs = 10,       # number of iterations/time step ?
                             timestep = NA) { # ?
  
  # x_new <- x # useless as R makes a local copy of x
  
  for (s in 1:subs) {
    
    # Get global average vegetation cover 
    rho_plus <- sum(x$cells == "+")/(x$dim[1]*x$dim[2]) 
    
    # Get the local proportion of disturbed neighbors 
    localempty <- count(x, "0") > 0 
    
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
