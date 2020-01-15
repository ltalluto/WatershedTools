.onUnload <- function (libpath) {
  library.dynam.unload("WatershedTools", libpath)
}

.onLoad <- function(libname, pkgname) {
	require(deSolve, quietly=TRUE)
	Sys.setenv("GRASS_VERBOSE"=0)
	err <- tryCatch(rgrass7::use_sp(), error = function(e) e)
}

#' @keywords internal
getGISBase <- function(appPath = "/Applications", suffix = "Contents/Resources") {
	if(!grepl("[Dd]arwin", Sys.info()['sysname']))
		stop("Only currently supported on mac")
	grassBase <- list.files(appPath, pattern='GRASS')
	gVersion <- as.numeric(sub(".*(7\\.[0-9\\.]+)\\.app", "\\1", grassBase))
	if(length(grassBase) > 1)
		grassBase <- grassBase[which.max(gVersion)]
	file.path(appPath, grassBase, suffix)
}
