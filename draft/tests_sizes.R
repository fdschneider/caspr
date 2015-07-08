#
# Tests for
# 

# Setup parallel backend
# -------------------------------------------------
library(doParallel)
cl <- makeCluster(5) # 5 cores
registerDoParallel(cl)

grid.size <- 25

# Loads all package functions
library(devtools)
load_all()

# One run example
# -------------------------------------------------

# Model : grazing
initial <- init_landscape(states = c("+","0","-"), 
                          cover = c(0.4,0.1,0.5), 
                          width = grid.size, 
                          height = grid.size)

parms <- grazing$parms

# Change if needed

system.time(
  grazingrun <- ca(initial, parms, model = grazing)
)

plot(grazingrun)
plot(grazingrun, plotstates = c(TRUE, TRUE, TRUE), 
     cols = c("darkgreen", "grey50", "grey80"), lwd = 2)

## a simple summary
summary(grazingrun)


# Example for a bifurcation graph
# -------------------------------------------------

covers.prob <- sample(list(high.branch = c(0.8,  0.1,  0.1),
                           low.branch  = c(0.05, 0.40, 0.40)), 1)

high.branch <- c(.8,.1,.1)

parms <- grazing$parms
# b is the bifurcation parameter
parms$b <- seq(0, 1, length.out = 10)

grazing.bifurc <- carray(grazing, high.branch, parms)
plot(grazing.bifurc[, 'b'], grazing.bifurc[, 'mean_cover.+'])


