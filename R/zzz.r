.onUnload <- function (libpath) {
  library.dynam.unload("WatershedTools", libpath)
}

.onLoad <- function(libname, pkgname) {
	require(deSolve, quietly=TRUE)
}

#' @keywords internal
getGISBase <- function(appPath = "/Applications", suffix = "Contents/Resources") {
	grassBase <- list.files(appPath, pattern='GRASS')
	gVersion <- as.numeric(sub(".*(7\\.[0-9\\.]+)\\.app", "\\1", grassBase))
	if(length(grassBase) > 1)
		grassBase <- grassBase[which.max(gVersion)]
	file.path(appPath, grassBase, suffix)
}
