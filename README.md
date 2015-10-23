R-package caspr
===============

The code provided here is running spatial disturbance models in a cellular automata framework. This is part of a collaborative project between the group of [Sonia Kéfi](http://sonia.kefi.fr/) (Institut de Sciences d'Evolution, CNRS, IRD, Université Montpellier, France) and [Vishwesha Guttal](https://teelabiisc.wordpress.com/) (Center for Ecological Sciences, Indian Institute of Science, Bangalore, India).

The R-package is build around objects of a particular class `ca_model` (S3 objects in R) that contain all the model specific information, like the original publication, the update functions, and the cell states. 
These objects can be handled by a function called `ca()` which runs the cellular automata over time. 

Other functions allow the plotting of single snapshots or timeseries, the generation of initial lattices/grids/landscapes, the calculation of spatial and temporal indicators as provided by the package ['spatialwarnings'](https://github.com/fdschneider/spatial_warnings).


## Contributors

Alain Danet, Alex Genin, Vishwesha Guttal, Sonia Kefi, Sabiha Majumder, Sumithra Sankaran, [Florian Schneider (Maintainer)](mailto:florian.schneider@univ-montp2.fr)

## Install package

To *use* the package, it can be installed directly from GitHub using the `devtools` package. 

```
install.packages("devtools")
devtools::install_github("fdschneider/caspr")
```

## Details 

Further description of the functions of the package, the models provided and the application of spatial indicators can be found in the [package vignette](http://htmlpreview.github.io/?https://github.com/fdschneider/caspr/blob/master/inst/doc/caspr.html) and the package [manual](https://github.com/fdschneider/caspr/blob/master/inst/doc/caspr-manual.pdf).

## Contribute to the package

Just a couple of guidelines for contributing code to the package:

1. The collaborative work on the package is coordinated via the [issues](https://github.com/fdschneider/caspr/issues). No uninvited pull-requests, please! Before you start developing a feature, please create an issue and assign yourself (if it does not exist yet). 
2. Clone the repository to your computer and work locally. If you add a new function or modify an existing one please include a valid documentation using [roxygen2](http://r-pkgs.had.co.nz/man.html). If you refer to functions of other packages, put them in the `Import` list in DESCRIPTION and use the structure `package::function()` to apply them in the code. 
3. Before you push your commits to GitHub, please be sure that everything works fine by testing the new functionality on your local repository. Test functions that rely on your function, too! If your feature works, create a last commit where you increase the version number in `DESCRIPTION` by `0.0.1` and mention the related issue in your commit message, e.g. `git commit -m "transfer code into R-package #1"`. If you're not listed yet, add your name to the `LICENSE` and `README.md` file. 


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
