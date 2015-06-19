R-package caspr
===============

The code provided here is running spatial disturbance models in a cellular automata framework. This is part of a collaborative project between the group of [Sonia Kéfi](http://sonia.kefi.fr/) (Institut de Sciences d'Evolution, CNRS, IRD, Université Montpellier, France) and [Vishwesha Guttal](https://teelabiisc.wordpress.com/) (Center for Ecological Sciences, Indian Institute of Science, Bangalore, India).

The R-package is build around objects of a particular class `ca_model` (S3 objects in R) that contain all the model specific information, like the original publication, the update functions, and the cell states. 
These objects can be handled by a function called `ca()` which runs the cellular automata over time. 

Other functions allow the plotting of single snapshots or timeseries, the generation of initial lattices/grids/landscapes, the calculation of spatial and temporal indicators.

## Install package

To *use* the package, it can be installed directly from GitHub using the `devtools` package. 

```
install.packages("devtools")
devtools::install_github("fdschneider/spatial_warnings")
devtools::install_github("fdschneider/caspr")
```

## Contribute to the package

Just a couple of guidelines for contributing code to the package:

1. The collaborative work on the package is coordinated via the [issues](https://github.com/fdschneider/caspr/issues). No uninvited pull-requests, please! Before you start developing a feature, please create an issue and assign yourself (if it does not exist yet). 
2. Clone the repository to your computer and work locally. If you add a new function or modify an existing one please include a valid documentation using [roxygen2](http://r-pkgs.had.co.nz/man.html).
3. Before you push your commits to GitHub, please be sure that everything works fine by testing the new functionality on your local repository. Test functions that rely on your function, too! If your feature works, create a last commit where you increase the version number in `DESCRIPTION` by `0.0.1` and mention the related issue in your commit message, e.g. `git commit -m "transfer code into R-package #1"`. If you're not listed yet, add your name to the `LICENSE` and `README.md` file. 


## Objects and functions

In the following the objects and functions are described in the sequence of their use (somewhat). 


####  function `init_landscape(states, cover, width, height)`

landscape objects are created using the function `init_landscape()` , e.g. 

```
l <- init_landscape(c("1", "0"), c(0.2,0.8), 100, 100)
``` 

#### `landscape` object class

A landscape object is a list `l` that contains

- `l$dims` : a named vector of the default form `c(width = 50, height = 50)`
- `l$cells` : a factorial vector of length `prod(l$dims)` that contains the state that each cell of the landscape matrix is in (row-wise from top to bottom). 

Specific method for functions `plot`, `summary` and `print` do exist, which means you can call:

```
l
plot(l)
summary(l)
```

The vectorial storage of the grid is allowing for a speedy evaluation than with matrices, using the map vectors provided globally by executing `mapping(l$dims[1], l$dims[2])`. 

For compatibility the landscape object can be translated into a matrix using `as.matrix(l)`. Any matrix of factorial content can be reverted into a landscape object using `as.landscape(l)`. 

#### `ca_model` object class 

A model object is a list `model` that contains

- `model$name`: the name of the model
- `model$ref`: the original reference
- `model$states`: the potential cell states
- `model$cols`: colors for cell states
- `model$parms`: a template for parameters or default parameters 
- `model$update`: the update function, which takes a landscape object `x_old` and returns a landscape object `x_new`, representing the updating of one single timestep

Those objects are not created by a function, but provided manually in separate files. If you add a model, following `model_template.R`, please test the update function for valid output and optimize for speed. 

(A method `print.ca_model` will be developed to quickly review a model's specifications.)

#### function `ca(x, model, parms)` 

The function `ca(x, model, parms)` runs a cellular automata simulation starting from the landscape object provided in `x` and using the update function provided by `model` with the parameter set `parms` (a list of parameters). Before running the model, the function validates the parameter set provided against the template stored in `model$parms`. 

Further parameters can be used to adjust the timespan run and the snapshots of the landscape saved. 

```
l <- init_landscape(c("+", "0", "-"), c(0.2,0.7,0.1))
p <- list(d = 0.4, r = 0.2, delta = 0.4) 
r <- ca(l, musselbed, p)

```

After simulation, the function returns an object of class `ca_results`.  


#### `ca_result` object class

The result object `r` of a simulation contains the full timeseries of cover and selected snapshots of the full landscape. It has the structure:

- `r$model`: the `ca_model` object used to run the simulation, including the provided parameter set.  
- `r$time`: a vector of the distinct timesteps of the simulation
- `r$evaluate`: the begin and end of a steady state period in the simulation run. 
- `r$cover`: a data frame reporting the global cover with one column for each state of the model and one row for each timestep. Thus, `r$cover[1]` returns the timeseries of the primary cell state. 
- `r$local`: a data frame reporting the average local cover of the cell states, i.e. the average probability that a state *i* is found in the neighborhood given that the focal cell is in state *i*. Thus, `r$local[1]` returns the timeseries of the primary cell state. 
- `r$snapshots`: a registry table of the snapshots and at which point in time they were taken. 
- `r$landscapes`: a list of landscape objects for the given snapshots.

#### function `indicators()`

#### `ca_indicators` object class

#### function `carray()`

## Contributors

Alain Danet, Alex Genin, Vishwesha Guttal, Sonia Kefi, Sumithra Sankaran, Florian Schneider (Maintainer), Sabiha Majumder

## License

The MIT License (MIT)

Copyright (c) 2015 the authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
