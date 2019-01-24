.onUnload <- function (libpath) {
  library.dynam.unload("WatershedTools", libpath)
}

.onLoad <- function(libname, pkgname) {
	require(deSolve, quietly=TRUE)
}