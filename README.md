# WatershedTools
R package for watershed analysis

This package wraps up some of the tools in GRASS GIS 7.4 to perform watershed analysis in R on R data structures.
The package relies heavily on data structures from the raster and sp packages for spatial objects on the R side.
On the grass side, it is necessary first to install GRASS GIS 7.4 on your system, then install this package.

To install, the easiest method is via github. You can do this by first installing `devtools` then using it to install this package:

    # install.packages('devtools') ## if needed
    library('devtools')
    install_github("mtalluto/WatershedTools")

The package will create a Grass session for you within your R session (though you can do so manually if you wish with `GrassSession()`; this can save some time).

Most functions work the way you expect other R functions to work; communication with Grass is handled transparently. In some cases, you will need to provide additional information, such as the location of your Grass installation. On some systems, this can be found by opening a command line window and running the command:

    grass74 --config path
    
## Installation Instructions

### Linux
Tested on Mint Cinnamon 19.1.

1. Install Grass GIS 7.4, including the grass-dev package. Verify that it is working with the command `grass74 --config path`.
2. Install the `r-base` and `r-base-dev` packages (and optionally RStudio).

### Windows
Tested on a fresh install of Windows 10.

1. Install [Base R](https://cran.r-project.org/bin/windows/base/) and [RTools](https://cran.r-project.org/bin/windows/Rtools/).
2. Install [Grass GIS 7.4](https://grass.osgeo.org/download/software/ms-windows/). Note that later versions will not work.
3. 


