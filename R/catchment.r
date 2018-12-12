#' Catchment deliniation
#' 
#' @param x An [sp::SpatialPoints] objector a matrix of x-y coordinates
#' @param drainage Raster or character, Drainage direction raster, from e.g., [accumulate()]. If specified as a character, a layer with that name from the existing [GrassSession()] given by gs will be used.
#' @param gs An optional [GrassSession()]; if missing a new one will be created
#' @param areas Logical; if `TRUE`, returns catchment area for each point, otherwise returns
#'  a [raster::RasterStack], one layer per point, with each layer delimiting the catchment for its
#'  respective point.
#' @param file The file name of the raster to be returned (if `areas` is `FALSE`), see `details`.
#' @param ... Additional parameters to pass to [GrassSession()]
#' @details This is a wrapper for [r.water.outlet](https://grass.osgeo.org/grass74/manuals/r.water.outlet.html)
#' 
#' It is recommended to specify the `file` parameter (including the extension to specify
#' file format; e.g., .tif, .grd). If not specified, a temp file will be created and will be
#' lost at the end of the R session
#'
#' @return A [raster::stack()] with two layers, flow accumulation ('accumulation') and drainage
#' 		direction ('drainage') (if outputName is missing), or a GrassSession otherwise



#' @details If `grass_session` is specified, the existing session will be used. In this case
#'   `drainage` is specified it will be written to the grass session with the name `drain_name`.
#'   If `drainage` is missing, then the grass session will be searched for an existing layer with
#'   a name given by `drain_name` (this is the fastest mode of operation).
#'
#' If no `grass_session` is specified, then a new session will be created, in which case `drainage`
#'   is required and will  be copied to the session.
#' 
#' This tool will perform much better if it is run with an existing grass session with 
#'	drainage direction already added:
#' 		`gs <- GrassSession(layer = drainage, layer_name = "drainage")`
#' 		`cAreas <- catchment(points, drain_name = "drainage", grass_session = gs, areas = TRUE)`
#' 	However note that one raster file is created *per point* in x, thus it is a good idea
#' 	with many points to specify area = TRUE to get area (rather than the raster) for each point
#' @return A vector of catchment areas (if `areas = TRUE`), otherwise a [raster::stack()] of 
#'   delineated catchments
catchment <- function(x, drainage, gs, areas = TRUE, file = NULL, ...)
{
	if(missing(gs)) {
		gs <- GrassSession(drainage, layerName = "drainage", ...)
		drainage <- "drainage"
	} else if(!is.character(drainage)) {
		gs <- GSAddRaster(drainage, layerName = "drainage", gs)
		drainage <- "drainage"
	}

	if(!is.matrix(x) & is.numeric(x)) {
		x <- matrix(x, ncol=2)
	} else if(!is.matrix(x)) {
		x <- sp::coordinates(x)
	}

	result <- if(areas) numeric(nrow(x)) else list()
	catchName <- "catchment"
	for(i in 1:nrow(x))
	{
		rgrass7::execGRASS("r.water.outlet", flags=c("overwrite"), input = drainage, 
			output = catchName, coordinates = x[i,])
		if(areas) {
			result[i] <- catchmentArea(catchName)
		} else {
			ras <- GSGetRaster(catchName, gs)
			ras[ras == 0] <- NA
			result <- c(result, ras)
		}
	}
	if(!areas) {
		result <- stack(result)
		if(!missing(file))
			result <- writeRaster(result, file)
	}
	return(result)
}


#' Compute the area of a delineated catchment
#' @param layer Character; the name of the grass catchment from which to compute area
#' @return Numeric; area of the catchment
#' @keywords internal
catchmentArea <- function(layer)
{
	vname <- paste0('v_', layer)

	rgrass7::execGRASS("r.to.vect", flags=c("overwrite", "quiet", "s"), 
		input = layer, output = vname, type='area', column = 'one')
	res <- rgrass7::execGRASS("v.to.db", flags=c("quiet", "p"), map = vname, 
				option = "area", intern=TRUE)
	idres <- as.numeric(sub("^(-?[0-9])+\\|(.+)", "\\1", res))
	res <- as.numeric(sub(".+?([0-9\\.]+)$", "\\1", res))
	val <- rgrass7::execGRASS("v.to.db", flags=c("quiet", "p"), map = vname, 
				option = "query", intern=TRUE, query_column="one")
	idval <- as.numeric(sub("^(-?[0-9])+\\|(.+)", "\\1", val))
	suppressWarnings(val <- as.numeric(sub(".+?([0-9\\.]+?)$", "\\1", val)))
	keep <- merge(data.frame(id = idres, area = res), 
		data.frame(id = idval, value = val), all.x = TRUE)

	## clean up
	rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), type="vector", name=vname)
	sum(keep$area[keep$value == 1])
}

