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
#' If no `grass_session` is specified, then a new session will be created, in which case `drainage`
#'   is required and will be copied to the session.
#' 
#' This tool will perform much better if it is run with an existing grass session with 
#'	drainage direction already added:
#' 		`gs <- GrassSession(layer = drainage, layer_name = "drainage")`
#' 		`cAreas <- catchment(points, drain_name = "drainage", grass_session = gs, areas = TRUE)`
#' 	However note that one raster file is created *per point* in x, thus it is a good idea
#' 	with many points to specify area = TRUE to get area (rather than the raster) for each point
#' @return A vector of catchment areas (if `areas = TRUE`), otherwise a [raster::stack()] of 
#'   delineated catchments
#' @export
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
			result[i] <- catchmentArea(catchName, gs)
		} else {
			ras <- GSGetRaster(catchName, gs)
			ras[ras == 0] <- NA
			result <- c(result, ras)
		}
	}
	if(!areas) {
		if(length(result) > 1) {
			result <- stack(result)
		} else {
			result <- result[[1]]
		}
		if(!missing(file))
			result <- writeRaster(result, file)
	}
	return(result)
}


#' Compute the area of a delineated catchment
#' @param layer Character; the name of the grass catchment from which to compute area
#' @param gs A [GrassSession()]
#' @return Numeric; area of the catchment
#' @keywords internal
catchmentArea <- function(layer, gs)
{
	vname <- paste0('v_', layer)

	GSRastToPoly(layer, vname, gs)
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
	gs <- GSClean(vname, gs, 'vector')
	sum(keep$area[keep$value == 1])
}


#' Crops a delineated stream to the catchment at a specified point
#' @param x The point at which to compute the catchment
#' @param streamRaster A RasterLayer stream delineation (required), or character giving the 
#' 		name of a Grass object
#' @param streamVector A SpatialLines stream delineation (optional), or character giving the 
#' 		name of a Grass object
#' @param drainage Raster or character, Drainage direction raster, from e.g., [accumulate()]. 
#' 		If specified as a character, a layer with that name from the existing [GrassSession()]
#' 		given by gs will be used.
#' @param gs An optional [GrassSession()]; if missing a new one will be created
#' @param file The file name of the raster to be returned (if `areas` is `FALSE`), see `details`.
#' @param trim Should the resulting layers be trimmed (cropped to the non-NA extent)?
#' @param ... Additional parameters to pass to [GrassSession()]
#' @return As [extractStream()]
#' @export
cropToCatchment <- function(x, streamRaster, streamVector, drainage, gs, file, trim = TRUE, ...)
{

	if(missing(gs))
	 	gs <- GrassSession(drainage, layerName = "cropToCatchment_drainage", ...)

	crCatchment <- catchment(x, drainage, gs, areas = FALSE)
	if(is.character(streamRaster))
		streamRaster <- GSGetRaster(streamRaster, gs)
	
	# for rasters we can do the cropping in R
	output <- raster::mask(streamRaster, crCatchment)
	if(trim)
		output <- raster::trim(output)
	if(!missing(file))
		output <- raster::writeRaster(output, file)

	if(!missing(streamVector)) {
		if(any(grepl('SpatialLines', class(streamVector)))) {
			rgrass7::writeVECT(streamVector, "cropToCatchment_streamVector", 
				v.in.ogr_flags = "overwrite", ignore.stderr=TRUE)
			streamVector <- "cropToCatchment_streamVector"
		}
		catchVname <- "cropToCatchment_catchArea"
		GSRastToPoly(crCatchment, catchVname, gs)
		
		oname <- "cropToCatchment_output"
		rgrass7::execGRASS("v.overlay", flags = c("overwrite", "quiet"), ainput = streamVector,
			atype = "line", binput = catchVname, operator = "and", output = oname)
		vout <- rgrass7::readVECT(oname)
		if(trim)
			vout <- raster::crop(vout, output)

		## clean up
		gs <- GSClean(catchVname, gs, "vector")
		gs <- GSClean(streamVector, gs, "vector")
		gs <- GSClean(oname, gs, "vector")

		## return lines
		output <- list(raster = output, vector = vout)
	}

	return(output)
}
