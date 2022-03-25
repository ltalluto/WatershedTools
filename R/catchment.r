#' Catchment deliniation
#' 
#' @param x An object of class [sp::SpatialPoints], [sf::sf], or a matrix of x-y coordinates
#' @param drainage Raster or character, Drainage direction raster, from e.g., [accumulate()]. If specified as a character, a layer with that name from the existing [GrassSession()] given by gs will be used.
#' @param gs An optional [GrassSession()]; if missing a new one will be created
#' @param areas DEPRECATED, use `output` instead. Logical; if `TRUE`, returns catchment area for each point, 
#'  otherwise returns
#'  a [raster::RasterStack], one layer per point, with each layer delimiting the catchment for its
#'  respective point.
#' @param file The file name of the raster to be returned if `output='raster'`, see `details`.
#' @param output One of 'area', 'raster', 'sf' determining the output format
#' @param overwrite If writing a raster to a file, should it be overwritten?
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
#' @return Depends on the value of the argument `output`:
#'  * 'area': a vector of catchment areas, one per point in x
#'  * 'raster': a [rasterLayer][raster::raster()] (either a single raster or a stack, if length(x) > 1)
#'  * 'sf': a polygon layer of class 'sf'
#' @export
catchment <- function(x, drainage, gs, areas, file = NULL, output = c('area', 'raster', 'sf'), 
		overwrite = FALSE, ...)
{
	output = match.arg(output)
	if(!missing(areas)) {
		warning("Argument areas is deprecated; use output instead")
		if(areas) {
			output = 'area'
		} else {
			output = 'raster'
		}
	}
	
	if(missing(gs)) {
		gs <- GrassSession(drainage, layerName = "drainage", ...)
		drainage <- "drainage"
	} else if(!is.character(drainage)) {
		gs <- GSAddRaster(drainage, layerName = "drainage", gs)
		drainage <- "drainage"
	}

	if(is(x, "Spatial")) {
		x = sp::coordinates(x)
	} else if(is(x, 'sf')) {
		x = sf::st_coordinates(x)
	} else if(!is.matrix(x) & is.numeric(x)) {
		x = matrix(x, ncol=2)
	}

	result <- list()
	catchName <- "catchment"
	for(i in 1:nrow(x))
	{
		rgrass7::execGRASS("r.water.outlet", flags=c("overwrite", "quiet"), input = drainage, 
			output = catchName, coordinates = x[i,])
		if(output == 'area') {
			result[[i]] = catchmentArea(catchName, gs)
		} else if(output == 'raster') {
			ras = GSGetRaster(catchName, gs)
			ras[ras == 0] = NA
			result = c(result, ras)
		} else if(output == 'sf') {
			if(nrow(x) > 1)
				stop("Output to sf is not supported with more than a single input point")
			vname = paste0(catchName, "_v")
			rgrass7::execGRASS("r.to.vect", flags = c('overwrite', 'quiet'), input = catchName, 
				output = vname, type = 'area', column='value')
			vect = sf::st_as_sf(rgrass7::readVECT(vname, ignore.stderr = TRUE))
			result[[i]] = vect
		}
	}
	
	if(output == 'area') {
		result = unlist(result)	
	} else if(output == 'raster') {
		if(length(result) > 1) {
			result = raster::stack(result)
		} else {
			result = result[[1]]
		}
		result = raster::trim(result)
		result = as.integer(result)
		if(!missing(file))
			result = writeRaster(result, file, overwrite = overwrite, 
				datatype="INT2U", options="COMPRESS=LZW")
	} else if(output == 'sf') {
		result = result[[1]]
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
#' @param streamVector A SpatialLines or sf stream delineation (optional), or character giving the 
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

	crCatchment <- catchment(x, drainage, gs, output = "raster")
	if(is.character(streamRaster))
		streamRaster <- GSGetRaster(streamRaster, gs)
	
	# for rasters we can do the cropping in R
	output <- raster::mask(streamRaster, crCatchment)
	if(trim)
		output <- raster::trim(output)
	if(!missing(file))
		output <- raster::writeRaster(output, file)

	if(!missing(streamVector)) {
		crCatchment <- catchment(x, drainage, gs, output = "sf")
		if(methods::is(streamVector, "SpatialLines")) {
			streamVector = sf::st_as_sf(streamVector)
		}
		
		vout = tryCatch(sf::st_intersection(streamVector, crCatchment), 
				error = function(e) {
					suppressWarnings(
						pts <- lapply(1:nrow(streamVector), function(x) sf::st_cast(streamVector[x,], "POINT"))
					)
					badrows = which(sapply(pts, nrow) == 1)
					if(length(badrows > 0)) {
						warning(length(badrows), " lines consisted of a single vertex and were removed")
						streamVector = streamVector[-badrows, ]
						sf::st_intersection(streamVector, crCatchment)
					} else {
						error(e)
		}})

		
		if(trim) {
			ext_poly = sf::st_as_sf(methods::as(raster::extent(output), "SpatialPolygons"))
			sf::st_crs(ext_poly) = sf::st_crs(output)
			ext_poly = sf::st_transform(ext_poly, sf::st_crs(vout))
			vout = sf::st_intersection(vout, ext_poly)
		}

		## return lines
		output <- list(raster = output, vector = methods::as(vout, "Spatial"))
	}

	return(output)
}
