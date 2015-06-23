#Notes on 23/06/2015
# Pred-prey model not working yet
# calling the function indicators does not work
# pb with the arguments of the function ca (tmin and tmax)


library(devtools)
document()#genere doc et recharge les functions

# To install the spatial_warnings package
#install_github("fdschneider/spatial_warnings")
#library(spatialwarnings) # not useful since already installed packages are automatically loaded by devtools

initial <- init_landscape(states = c("+","0","-"), cover = c(0.4,0.1,0.5), width = 25, height = 25)

plot(initial, cols = c("darkgreen", "grey70", "white")) 

#-------------------------------------------------------------------
##### Run Flo's model
#-------------------------------------------------------------------

parms_grazing <- list(
  del = 0.9,
  b = 0.3,
  c_ = 0.2,
  m0 = 0.05,
  g = 0.1,
  r = 0.01,
  f = 0.9,
  d = 0.2,
  protect = 0.9
) 

grazingrun <- ca(initial, parms_grazing, model = grazing)

plot(grazingrun)
plot(grazingrun, plotstates = c(TRUE, TRUE, TRUE), cols = c("darkgreen", "grey50", "grey80"), lwd = 2)

## a simple summary
summary(grazingrun)

## a function that returns  indicators / early warning signs
#indicators(grazingrun)

#-------------------------------------------------------------------
##### Run mussel bed's model
#-------------------------------------------------------------------
parms_mussel  <- list(
  r = 0.4, # recolonisation of empty sites dependent on local density
  d = 0.9, # probability of disturbance of occupied sites if at least one disturbed site
  delta = 0.01 # intrinsic disturbance rate
) 

musselrun <- ca(initial, parms_mussel, model = musselbed)

## specify what to plot 
plot(musselrun, plotstates = c(TRUE, FALSE, TRUE), cols = c("black", "grey50", "grey90"), lwd = 2)

#-------------------------------------------------------------------
##### Run forest gap's model
#-------------------------------------------------------------------

parms_forestgap  <- list(
  alpha = 0.20, # recovery strength 
  delta = 0.17, # strength of gap expansion
  d = 0.01      # intrinsic death rate
) 

forestgaprun <- ca(initial, parms_forestgap, model = forestgap)

## specify what to plot 
plot(forestgaprun, plotstates = c(TRUE, TRUE), cols = c("black","grey90"), lwd = 2)

#-------------------------------------------------------------------
##### Run pred-prey model
#-------------------------------------------------------------------