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
devtools::install_github(fdschneider/caspr)
```

## Contribute to the package

Just a couple of guidelines for contributing code to the package:

1. The collaborative work on the package is coordinated via the [issues](https://github.com/fdschneider/caspr/issues). No uninvited pull-requests, please! Before you start developing a feature, please create an issue and assign yourself (if it does not exist yet). 
2. Clone the repository to your computer and work locally. If you add a new function or modify an existing one please include a valid documentation using [roxygen2](http://r-pkgs.had.co.nz/man.html).
3. Before you push your commits to GitHub, please be sure that everything works fine by testing the new functionality on your local repository. Test functions that rely on your function, too! If your feature works, create a last commit where you increase the version number in `DESCRIPTION` by `0.0.1` and mention the related issue in your commit message, e.g. `git commit -m "transfer code into R-package #1"`. If you're not listed yet, add your name to the `LICENSE` and `README.md` file. 


## Objects and functions

### object structure

#### `landscape`

#### `ca_model` 

A model object is a list `x` that contains

- `x$name`: the name of the model
- `x$ref`: the original reference
- `x$states`: the potential cell states
- `x$cols`: colors for cell states
- `x$parms`: a template for parameters or default parameters 
- `x$update`: the update function, which takes a landscape object `x_old` and returns a landscape object `x_new`, representing the updating of one single timestep

Those objects are not created by a function, but provided manually in separate files. If you add a model, please test the update function for valid output and optimize for speed. 

#### `ca_result`


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

## License

The MIT License (MIT)

Copyright (c) 2015 Florian D. Schneider

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
