# 
# Generate the data needed for the analyses 

# Setwd
setwd('/home/alex/all-data/ifcam_simulations/')
# setwd('/home/alex/system/tmp')

# 
# # Use this to use the proxy at IISC
# Sys.setenv(https_proxy = "proxy.iisc.ernet.in:3128")
# Sys.setenv(http_proxy  = "proxy.iisc.ernet.in:3128")


# Load package/install if necessary
library(devtools)
install_github('fdschneider/caspr') # might need a restart of R
library(caspr)
library(ggplot2)
library(doParallel)


# Global variables
SIZE   <- 100
RES    <- 4 # Number of points on each dimension
NSNAPS <- 10
REDO_COMPUTATIONS <- TRUE
DO_PREDPREY <- FALSE
LENGTH.STAT <- 100 # Number of time steps to consider when computing mean covers

# Computation variables
PARALLEL <- FALSE
NCORES   <- 15

# TODO: 
#   - add time series 
#   - set t_min to 500 (done)
#   - set t_max to 3000 (done)
#   - forestgap: d btw 0 and 0.25, delta between 0 and 0.25 (done)
#   - musselbed: removed the value \delta = 0 from plots (nothing to do here)
#   - decide on whether to use b or m in the grazing model
#   - musselbed: 0 to 0.3 (done)
#   - musselbed: catastrophic shifts along the b parameter ? 
#   - predprey: problem with oscillations: 
#       - take more snapshots ? 
#       - more time steps ?
#   - grazing: use g and m (done but make a quick simulation before to 
#       check for parameter ranges)
# 
# For the slices : 
#   - forestgap: delta \in 0, 0.05, 0.15; 
#   - 

# Register parallel backend
if (PARALLEL) { 
  registerDoParallel(cores = NCORES)
} else { 
  registerDoSEQ()
}

# registerDoSEQ()

# Produce data needed for the 2D state diagram
# --------------------------------------------------------
if (REDO_COMPUTATIONS) {
  
  # Forestgap model -------------
  
  parms <- list(alpha = 0.2,
                delta = seq(0, .25, length.out = RES),
                d     = seq(0, .25, length.out = RES))
  
  result_forestgap_upper <- ca_arraySS(forestgap, 
                                       init = c(.95, .05), 
                                       width = SIZE, height = SIZE, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT)
  
  result_forestgap_lower <- ca_arraySS(forestgap, 
                                       init = c(.05, .95), 
                                       width = SIZE, height = SIZE, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT)
  
  save(result_forestgap_lower, result_forestgap_upper, 
       file = "./result_forestgap.rda", compress = 'gzip2')
  
  rm(result_forestgap_lower)
  rm(result_forestgap_upper)
  
  # Pred-prey model --------------
  if ( DO_PREDPREY ) { 
  # NOTE: subs is disabled (forced to one) for now in the predator-prey 
  #   model !!!
  # 
  # NOTE: due to oscillations and "slow" behavior we increase t_min, t_max and
  #   length.stat so stochasticity does not affect so much the final mean covers.
  # 
  # NOTE: we can't have low predator and high prey so the upper branch is 
  #   necessarily high prey *no* predator, lower branch could be coexistence
  #  - threes domains: 
  #    - coexistence
  #    - extinction of predators (=> full prey)
  #    - coextinction
  #  - we're interested in the transition coexistence -> coextinction
  # 
  #  - We need an idea of the stability landscape: 
  #     - start from upper branch (prey at eq level without predator) and 
  #         add predators: how much do you need to actually go into coexistence
  #         regime given a set of parameters ? 
  #     - implementation: run from full prey to eq, then introduce predators
  #  - As soon as predator survives you get oscillations (=> Hopf bifurcation ?)
  #  - 50% and 5% (for pred&prey)
  #  - Mean cover over 1000 time steps
  #  - Decide on a cutoff to decide where the critical point is
  #  - Instead of a 2D diagram :
  #      - plot only a few bifurcation diagrams
  #        - Hopf along \delta axis (?): delta= 0.01, betaf = 1/3, betas = .2, m = 0.1
  #        - Still to be found along m axis
  #  - Time-series needed for indicators :
  #   - maybe focus on spectrum rather than var/skew indics ?
    parms <- list(betaf = 1/3,
                  betas = seq(0, .5, length.out = RES),
                  delta = 0.01, # .01 to .02
                  m     = seq(0, .2, length.out = RES))
    
    result_predprey_lower <- ca_arraySS(predprey, 
                                        init  = c(.05, .25, .70), 
                                        width = SIZE, height = SIZE, 
                                        parms = parms, nsnaps = NSNAPS, 
                                        length.stat = LENGTH.STAT * 5, 
                                        t_min = 500, 
                                        t_max = 2000)
    
    result_predprey_upper <- ca_arraySS(predprey, 
                                        init  = c(.75, .25, 0), 
                                        width = SIZE, height = SIZE, 
                                        parms = parms, nsnaps = NSNAPS, 
                                        length.stat = LENGTH.STAT * 5, 
                                        t_min = 500, 
                                        t_max = 2000)
    
    save(result_predprey_lower, result_predprey_upper, 
        file = "./result_predprey.rda", compress = 'gzip2')
    rm(result_predprey_lower)
    rm(result_predprey_lower)
  
  }
  
  # Grazing model -------------
  parms <- list(del = 0.9,
                b   = 0.5, # suggested default value in the model
                c_  = 0.2,
                m0  = seq(0,  1, length.out = RES),
                g   = seq(0, .8, length.out = RES),
                r   = 0.01,
                f   = 0.9,
                d   = 0.1,
                p   = 1)
    
  result_grazing_upper <- ca_arraySS(grazing, 
                                     init = c(.99, .01, 0), 
                                     width = SIZE, height = SIZE, 
                                     parms = parms, nsnaps = NSNAPS,
                                     length.stat = LENGTH.STAT)
  
  result_grazing_lower <- ca_arraySS(grazing, 
                                     init = c(.05, .95, 0), 
                                     width = SIZE, height = SIZE, 
                                     parms = parms, nsnaps = NSNAPS,
                                     length.stat = LENGTH.STAT)
  
  save(result_grazing_lower, result_grazing_upper, 
       file = "./result_grazing.rda", compress = 'gzip2')
  rm(result_grazing_upper)
  rm(result_grazing_lower)
  
  # Musselbed model -------------
  parms <- list(r     = 0.4,
                d     = seq(0, 1,  length.out = RES),
                delta = seq(0, .3, length.out = RES))
                
  result_musselbed_upper <- ca_arraySS(musselbed, 
                                       init = c(.99, .01, 0), 
                                       width = SIZE, height = SIZE, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT)
  
  result_musselbed_lower <- ca_arraySS(musselbed, 
                                       init = c(.05, .95, 0), 
                                       width = SIZE, height = SIZE, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT)
  
  save(result_musselbed_lower, result_musselbed_upper, 
       file = "./result_musselbed.rda", compress = 'gzip2')
  rm(result_musselbed_lower)
  rm(result_musselbed_upper)
  
  
  # Save results --------------
  
} else { # do not redo simulations
  load("./data_simulations.rda")
}

# Merge and plot results
make_df <- function(result_up, result_low) { 
  data.frame(branch = c(rep("upper", nrow(result_up[['DatBif']])), 
                        rep("lower", nrow(result_low[['DatBif']]))), 
             rbind(result_up[['DatBif']], 
                   result_low[["DatBif"]]))
}

ggplot(make_df(result_musselbed_upper, result_musselbed_lower)) + 
  geom_raster(aes(d, delta, fill = mean_cover_.)) + 
  facet_grid( ~ branch)

ggplot(make_df(result_predprey_upper, result_predprey_lower)) + 
  geom_raster(aes(m, betas, fill = mean_cover_f)) + 
  facet_grid( ~ branch)

ggplot(make_df(result_forestgap_upper, result_forestgap_lower)) + 
  geom_raster(aes(d, delta, fill = mean_cover_.)) + 
  facet_grid( ~ branch)

ggplot(make_df(result_grazing_upper, result_grazing_lower)) + 
  geom_raster(aes(b, g, fill = mean_cover_.)) + 
  facet_grid( ~ branch)


