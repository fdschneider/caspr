README for R-package caspr
==========================

The code provided here is running spatial disturbance models in a cellular automata framework. 

The package is based on objects of class `ca_model` (S3 objects in R) that contain all the model specific information, like the original publication, the update functions, and the cell states. 
These objects can be handled by a function called `ca()` which runs the cellular automata over time. 

Other functions allow the plotting of single snapshots or timeseries, the generation of initial lattices/grids/landscapes, the calculation of spatial and temporal indicators.

### object structure

A model object is a list `x` that contains

- `x$name`: the name of the model
- `x$ref`: the original reference
- `x$states`: the potential cell states
- `x$cols`: colors for cell states
- `x$parms`: a template for parameters or default parameters 
- `x$update`: the update function, which takes a landscape object `x_old` and returns a landscape object `x_new`, representing the updating of one single timestep

### functions

#### `init_landscape(states, cover, width, height)`

returns an object `i` of class `landscape` that contains the dimensions (`x$dims`) and the cell contents in a character vector (`i$cells`). Specific method for functions `plot`, `summary` and `print` do exist, which means you can call:

```
i <- init_landscape(c("1", "0"), c(0.2,0.8), 100, 100)
i
plot(i)
summary(i)
```

#### `ca()`

#### `indicators()`

