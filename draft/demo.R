## Example

document()

# creating an initial landscape object

initial <- init_landscape(states = c("+","0","-"), cover = c(0.4,0.1,0.5))

plot(initial) 
plot(initial, cols = c("darkgreen", "grey70", "white")) 

summary(initial)


# change size of lattice

initial <- init_landscape(states = c("+","0","-"), cover = c(0.4,0.1,0.5), width = 25, height = 25)
 
plot(initial, cols = c("darkgreen", "grey70", "white")) 

summary(initial)




# run simulation

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

## specify what to plot 
plot(grazingrun, plotstates = c(TRUE, TRUE, TRUE), cols = c("darkgreen", "grey50", "grey80"), lwd = 2)


## get plots of the snapshots
par(mfrow = c(2,5), mar = c(0,0,0,0)+0.4)
for(i in 1:9) {
  plot(grazingrun$timeseries[[i]])
}

dev.off()
# analyse output

## a simple summary
summary(grazingrun)

## a function that returns  indicators / early warning signs
indicators(grazingrun)



# select model

parms_mussel  <- list(
  r = 0.4, # recolonisation of empty sites dependent on local density
  d = 0.9, # probability of disturbance of occupied sites if at least one disturbed site
  delta = 0.01 # intrinsic disturbance rate
) 

musselrun <- ca(initial, parms_mussel, model = musselbed)


## specify what to plot 
plot(musselrun, plotstates = c(TRUE, FALSE, TRUE), cols = c("black", "grey50", "grey90"), lwd = 2)


## get plots of the snapshots
par(mfrow = c(2,5), mar = c(0,0,0,0)+0.4)
for(i in 1:9) {
  plot(musselrun$timeseries[[i]])
  mtext(musselrun$timeseries[[i]], outer = TRUE)
}
dev.off()


# simulate gradients



parms_mussel_gradient  <- list(
  r = 0.4, # recolonisation of empty sites dependent on local density
  d = seq(0,1,0.1), # probability of disturbance of occupied sites if at least one disturbed site
  delta = 0.01 # intrinsic disturbance rate
) 


# provides parallel backend
library(foreach)
library(doSNOW)

workerlist <- c(rep("localhost", times = 7)) 

cl <- makeSOCKcluster(workerlist, outfile='out_messages.txt')

registerDoSNOW(cl)

musselgradient <- ca_array(initial, parms_mussel_gradient, model = musselbed)

stopCluster(cl)

musselgradient

plot(`+`~ d, data = musselgradient )
