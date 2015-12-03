# 
# This file describes a possible workflow when constructing bifurcation 
#   diagrams for a model (here, predprey)
# 

library(caspr)   # use install_github('fdschneider/caspr') if not installed
library(ggplot2) # plots
library(tidyr)   # reformat data frames
library(doParallel) # parallel computing

# Global variables
XRES <- 20   # Number of points along the gradient of parameter
NREP <- 10   # Number of replicates for each parameter value
SIZE <- 128  # Size (width = height) of the lattice
STAT_NTSTEP <- 100 # Number of samples to take to compute the means/sds 
                   #   of covers at the end of simulations
TMIN <- 500 # Minimal number of iterations to consider

# Register a parallel backend for ca_array
#   ca_array will not run if there are too many simulations to be done without 
#   having a parallel backend registered.
registerDoParallel(cores = 2) # set you number of cores to use here (
                              #   hint: put this to 23 on the workstation)
# registerDoSEQ() # do not use parallel computing (turns off the above setting)





# 
# Example 1: run a single simulation, access time series and saved 
#   snapshots.
# ------------------------------------

# Create a lattice of 100x100
init <- init_landscape(states = predprey$states, 
                       cover = c(.01, .25, .74), 
                       width = 100,
                       height = 100)

# Explicitely build a list of parameters
parms <- list(betaf = 1/3, 
              betas = .2, # .1 -> coexistence; extinction pour .3 (0-0.3 ? )
              delta = .01, # delta = .03 => 0% pred
              m     = .01) #  m = 0 -> 0.2

# Compute simulation
# 
result <- ca(init, 
             model = predprey, 
             parms = parms, 
             t_max = 1000, 
             stopifsteady = FALSE)

# Format data, build and display the plot object
# We can directly call the plot function on the result (plot(result)), 
#  but this way is more flexible at the cost of being more complex
# 
# plot(result) # check out ?plot.ca_result for more info
#
plot_df <- data.frame(time = result$time, result$cover) # extract time series
plot_df <- gather_(plot_df, "state", "cover", # reformat data in ggplot format
                   gather_cols = c("f", "s", "X0"))

plot <- ggplot(plot_df) +  # construct plot object
  geom_line(aes(time, cover, color = state)) 

print(plot) # display it

# Visualize a single landscape 
# /!\ The number in brackets must be < t_max
plot(result[['landscapes']][[30]], cols = predprey$cols)






# Example 2: run a bunch of simulations to create a bifurcation diagram
# ------------------------------------

# 2.1 over the well-mixed stressor
# --------

# We create a bifurcation diagram for the m values
parms   <- list(betaf = 1/3,   # aka alpha1 or beta1 (fish repro)
                betas = .2,    # aka alpha2 or beta2 (shark/tiger repro)
                delta = .01,   # starving shark proba of death
                # m: 20 values over the range between 0 and .2, with 10 
                # replicates at each value
                m     = rep(seq(0, .2, length.out = XRES), each = NREP))

bifdiag_data_varying_m_upper <- ca_array(model = predprey, 
                                         init  = c(.75, .25, 0), 
                                         width = SIZE, height = SIZE,
                                         parms = parms, t_max = 1000,
                                         stat_covers_length = STAT_NTSTEP, 
                                         t_min = TMIN)

bifdiag_data_varying_m_lower <- ca_array(model = predprey, 
                                         init  = c(.01, .25, .74), 
                                         width = SIZE, height = SIZE,
                                         parms = parms, t_max = 1000, 
                                         stat_covers_length = STAT_NTSTEP, 
                                         t_min = TMIN)

# Merge both results in a data frame, with a column specifying which branch 
#   the result belongs to 
bifdiag_data_varying_m <- 
  data.frame(branch = c(rep('lower', nrow(bifdiag_data_varying_m_lower)),
                        rep('upper', nrow(bifdiag_data_varying_m_upper))),
             rbind(bifdiag_data_varying_m_lower, 
                   bifdiag_data_varying_m_upper))

# Save the results
save(bifdiag_data_varying_m, 
     file = './bifdiag_data_varying_m.rda')


# Build a plot 
# Format data frame
plot_df <- gather(bifdiag_data_varying_m, sp, cover, 
                  mean_cover_f, mean_cover_s)

bifplot <- ggplot(plot_df) + 
  geom_point(aes(m, cover, color = sp)) + 
  facet_grid( ~ branch )
print(bifplot)





# 2.1 over the spatial stressor
# --------

# We create a bifurcation diagram for the betas param (alpha2/beta2)
parms   <- list(betaf = 1/3,   # alpha1 
                betas = rep(seq(0, .3, length.out = XRES), each = NREP),
                delta = .01, # delta
                m     = .01)

bifdiag_data_varying_betas_upper <- ca_array(model = predprey, 
                                             init  = c(.75, .25, 0), 
                                             width = SIZE, height = SIZE,
                                             parms = parms, t_max = 1000, 
                                             stat_covers_length = 100, 
                                             t_min = TMIN)

bifdiag_data_varying_betas_lower <- ca_array(model = predprey, 
                                             init  = c(.01, .25, .74), 
                                             width = SIZE, height = SIZE,
                                             parms = parms, t_max = 1000,
                                             stat_covers_length = 100, 
                                             t_min = TMIN)

bifdiag_data_varying_betas <- 
  data.frame(branch = c(rep('lower', nrow(bifdiag_data_varying_betas_lower)),
                        rep('upper', nrow(bifdiag_data_varying_betas_upper))),
             rbind(bifdiag_data_varying_betas_lower, 
                   bifdiag_data_varying_betas_upper))

save(bifdiag_data_varying_betas, 
     file = './bifdiag_data_varying_betas.rda')

# Format data 
plot_df <- gather(bifdiag_data_varying_betas, sp, cover, 
                  mean_cover_f, mean_cover_s)

# Make plot
bifplot <- ggplot(plot_df) + 
  geom_point(aes(betas, cover, color = sp)) + 
  facet_grid( ~ branch )
print(bifplot)




# Change the two bifurcation parameters and create a state diagram
# ------------------------------------


# We create a bifurcation diagram for the betas param (alpha2/beta2)
parms   <- list(betaf = 1/3,   # alpha1 
                betas = rep(seq(0, .3, length.out = 2), each = 1),
                delta = .01, # delta
                m     = rep(seq(0, .2, length.out = 2), each = 1))

# Check the number of simulations with : 
# nrow(expand.grid(parms))

statediag_data_upper <- ca_array(model = predprey, 
                                 init  = c(.75, .25, 0), 
                                 width = SIZE, height = SIZE,
                                 parms = parms, t_max = 1000, 
                                 stat_covers_length = 100, 
                                 t_min = TMIN)

statediag_data_lower <- ca_array(model = predprey, 
                                 init  = c(.01, .25, .74), 
                                 width = SIZE, height = SIZE,
                                 parms = parms, t_max = 1000,
                                 stat_covers_length = 100, 
                                 t_min = TMIN)

statediag_data <- 
  data.frame(branch = c(rep('lower', nrow(statediag_data_lower)),
                        rep('upper', nrow(statediag_data_upper))),
             rbind(statediag_data_lower, 
                   statediag_data_upper))

save(statediag_data, file = './statediag_data.rda')

# Format data
plot_df <- gather(statediag_data, sp, cover, 
                  mean_cover_f, mean_cover_s)

# Make plot
stateplot <- ggplot(plot_df) + 
  geom_point(aes(betas, cover, color = sp)) + 
  facet_grid( ~ branch )
print(bifplot)

