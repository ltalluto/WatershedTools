
#' Produce a [sp::SpatialPixelsDataFrame] from a [raster::RasterLayer]
#'
#' @param x A [raster::RasterLayer] object
#' @return A [sp::SpatialPixelsDataFrame]
#' @keywords internal
rasterToSPDF <- function(x)
{
	coords <- sp::coordinates(x)
	gr <- data.frame(x=coords[,1], y=coords[,2], val=raster::values(x))
	sp::coordinates(gr) <- c(1,2)
	sp::proj4string(gr) <- sp::proj4string(x)
	sp::gridded(gr) <- TRUE
	return(gr)
}


#' Set up a GRASS session
#'
#' @param layer A [raster::raster] object (recommended); or any other object which has `extent()` and `proj4string()` methods defined.
#' @param gisBase character; the location of the GRASS installation (see `details`)
#' @param layer_name NA or character; if NA (the default), the layer will not be added to the grass session, otherwise it will be appended with this name.
#' @param home Location to write GRASS settings, `details`.
#' @param gisDbase Location to write GRASS GIS datasets; see `details`.
#' @param location Grass location name
#' @param mapset Grass mapset name
#' @param override Logical;  see `details`
#'
#' @details if `gisBase` is not provided, it can be automatically deduced in some cases.
#' On some systems, you can run `grass74 --config path` from the command line to get this path.
#' 
#' The extent, projection, and resolution of the grass session will be determined by `layer`.
#' 
#' by default `override` will be TRUE if home and gisDbase are set to their defaults, otherwise FALSE. If TRUE, the new session will override any existing grass session (possibly damaging/overwriting existing files). It is an error if override is FALSE and there is an already running session.
#' @return An S3 [GrassSession] object
GrassSession <- function(layer, gisBase, layer_name = NA, home = tempdir(), gisDbase = home,location = 'NSmetabolism', mapset = 'PERMANENT', override)
{
	if(missing(gisBase)) 
		gisBase <- system2("grass74", args=c("--config path"), stdout=TRUE)

	if(missing(override))
	{
		if(home == gisDbase  && home == tempdir())
		{
			override <- TRUE
		} else
			override <- FALSE
	}

	gs <- list()
	rgrass7::initGRASS(gisBase, home=home, gisDbase = gisDbase, location = location, 
		mapset = mapset, override = override)
	gs$gisBase <- gisBase
	gs$home <- home
	gs$gisDbase <- gisDbase
	gs$location <- location
	gs$mapset <- mapset

	err <- rgrass7::execGRASS("g.proj", flags = "c", proj4 = sp::proj4string(layer), intern=TRUE)
	gs$proj4string <- sp::proj4string(layer)

	ext <- as.character(as.vector(raster::extent(layer)))
	rasres <- as.character(raster::res(layer))
	rgrass7::execGRASS("g.region", n = ext[4], s = ext[3], e = ext[2], w = ext[1], 
			rows=raster::nrow(layer), cols=raster::ncol(layer), nsres = rasres[2], 
			ewres = rasres[1])
	gs$extent <- raster::extent(layer)
	gs$resolution <- raster::res(layer)

	class(gs) <- c("GrassSession", class(gs))

	if(!is.na(layer_name))
		gs <- GSAddRaster(layer, layer_name, gs)

	return(gs)
}

#' Print method for a [GrassSession] object
#'
#' @param x A [GrassSession] object
#'
print.GrassSession <- function(x)
{
	print(rgrass7::gmeta())
}

#' Add a RasterLayer to a grass session and return the modified session
#' @param x A [raster::raster] object or a [sp::SpatialGridDataFrame]
#' @param grass_session A [GrassSession] object
#' @param layer_name character; the name of the layer to add to grass.
#'
#' @return An S3 [GrassSession] object
GSAddRaster <- function(x, layer_name, grass_session, overwrite = TRUE)
{
	if("RasterLayer" %in% class(x))
		xSPDF <- rasterToSPDF(x)
	flags <- NULL
	if(overwrite)
		flags <- c(flags, "overwrite")
	rgrass7::writeRAST(xSPDF, layer_name, flags = flags)
	grass_session$layers <- c(grass_session$layers, layer_name)
	return(grass_session)
}
