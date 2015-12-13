# 
# This is an example file that illustrates the way to load, extract and 
#   process the output from the caspr simulations
# 

# Adjust to the folder containing the data on _your_ computer 
setwd("/home/alex/work/2014-2015/SpatialStress/caspr/draft/ifcam_results/")

# Use this to use the proxy at IISC
Sys.setenv(https_proxy = "proxy.iisc.ernet.in:3128")
Sys.setenv(http_proxy  = "proxy.iisc.ernet.in:3128")

# We will need the following packages
library(devtools)
library(ggplot2)


# Install the spatialwarnings package if it is not already installed (this 
# sometimes require restarting R).
if ( ! "spatialwarnings" %in% row.names(installed.packages())) { 
  install_github('fdschneider/spatial_warnings') # mind the "_"
}
library(spatialwarnings)

# Here we'll take the forestgap model as example
# 
# Let's load the data. /!\ It takes about ~2 to 3G of memory
load('./result_forestgap_processed.rda', verbose = TRUE)

# This loads two variables: upper_branch and lower_branch, corresponding 
#   to different initial conditions. 
# 
# Each of those variables is a list of two components : 
#   
#   1. DatBif: A data.frame of 1600 rows ( = 40*40) with the following columns
#     
#     ID            the simulation number
#     alpha         |
#     delta         | Simulation parameters 
#     d             |
#     seed          seed of RNG
#     mean_cover_.  cover of +          (note that the + is replace by a . as 
#     mean_cover_0  cover of 0             data.frame column names cannot have 
#     sd_cover_.    sd cover of +          special characters)
#     sd_cover_0    sd cover of 0
#     t_start     
#     t_end       
#     is_steady     steadyness
# 
# For example, we can plot the density of trees `+` in the space parameter : 
ggplot(upper_branch$DatBif) + 
  geom_raster(aes(x = d, y = delta, fill = mean_cover_.))
# 
#   2. A list of 1600 components containing 10 snapshots for each simulation
# 
# For example, let's plot the simulation with ID 367, we can access its parameters
#   by subsetting the summary data frame.
simulation_n <- 367
with(upper_branch, DatBif[DatBif$ID == simulation_n, ])

# Let's access the output of the first simulation
snapshot_n <- 1
extracted_matrix <- upper_branch[[2]][[simulation_n]][[snapshot_n]] 
image(extracted_matrix)

# Let's plot all the snapshots for that simulation
par(mfrow=c(2,5))
for (i in seq.int(10)) { 
  image(upper_branch[[2]][[simulation_n]][[i]])
}




# We can now use functions from the spatialwarnings package. As they work 
# on lists we can directly feed them all the snapshots of one simulation. 
# 
# e.g. for variance
var_indic <- indicator_variance(upper_branch[[2]][[simulation_n]], 
                                nreplicates = 0, # turns off shuffling 
                                                 #   (makes no sense here)
                                subsize = 2)
# 
# e.g. for all generic ews (requires spatialwarnings)
ews_indic <- generic_spews(upper_branch[[2]][[simulation_n]], 
                           subsize = 2)
plot(ews_indic)


# Let's merge the two branches together in one table
merge_two_branches <- function(up, low) { 
  data.frame(branch = c(rep("upper", nrow(up$DatBif)), 
                        rep("lower", nrow(low$DatBif))), 
             rbind(up$DatBif, low$DatBif))
} 
merged_results <- merge_two_branches(upper_branch, lower_branch)


# Let's extract a 1D bifurcation diagram
closest_to <- function(val, X) abs(X - val) == min(abs(X - val))
data_subset <- subset(merged_results, closest_to(.0, d))

ggplot(data_subset) + 
  geom_point(aes(delta, mean_cover_., color = branch)) 

# Let's plot a 2D state diagram with both branches 
ggplot(merged_results) + 
  geom_raster(aes(x = delta, 
                  y = d, 
                  fill = mean_cover_.)) + 
  facet_grid( ~ branch) + 
  scale_fill_gradient(low = "#564918" , high = "#20BB20")


# Do a graph "à la flo" (hard to read for now, needs more tweaking)
ggplot() + 
  geom_raster(aes(x = delta, y = d), fill = "grey70",
              data = subset(lower_branch$DatBif, mean_cover_. < .25)) + 
  geom_contour(aes(x = delta, y = d, z = mean_cover_.), 
               data = subset(upper_branch$DatBif), 
               color = "black", linesize = .1) 

# Plot a surface (not very useful)
library(lattice) 
wireframe(mean_cover_. ~ d * delta, 
          data = upper_branch$DatBif, 
          screen = list(z = 700, x = -75))
