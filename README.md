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


## Installation Notes

### macOS

Developer tools are necessary to compile the package from source, which presently is the only way to install the package. I suggest you follow the [instructions on CRAN](https://cran.r-project.org/bin/macosx/) for installing tools necessary to compile packages from source. You will also need to install xcode from Apple.

If you are getting compilation errors after installing xcode, be sure that you have installed the developer tools. This can be accomplished by opening a terminal window and running `xcode-select --install`. Sometimes it is necessary to do this step again if you have recently updated xcode or your macOS installaion.