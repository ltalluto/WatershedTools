#' Produce a DEM with basins filled in, as well as a flow direction raster
#'
#' @param dem Raster or character; a digital elevation model. If specified as a character, a layer with that name from the existing [GrassSession()] given by gs will be used.
#' @param gs An optional [GrassSession()]; if missing a new one will be created
#' @param filledDEM Optional character; if missing, a raster stack is returned, otherwise the filled DEM is written to the grass session with this name and the modified GrassSession is returned.
#' @param probs Optional character; as with filled DEM, for problem areas
#' @param file The file name of the raster to be returned (if no `outputName` is provided), see `details`.
#' @param ... Additional parameters to pass to [GrassSession()]
#'
#' @details This is a wrapper for [r.fill.dir](https://grass.osgeo.org/grass74/manuals/r.fill.dir.html)
#' 
#' It is recommended to specify the `file` parameter (including the extension to specify
#' file format; e.g., .tif, .grd). If not specified, a temp file will be created and will be
#' lost at the end of the R session
#'
#' @return A [raster::stack], optionally written to `file` with 2 layers: the filled DEM and a list of problem areas (if outputName is missing), or a GrassSession otherwise
fillDEM <- function(dem, gs, filledDEM, probs, file = NULL, ...)
{
	if(missing(gs)) {
		gs <- GrassSession(dem, layerName = "dem", ...)
		dem <- "dem"
	} else if(!is.character(dem)) {
		gs <- GSAddRaster(dem, layerName = "dem", gs)
		dem <- "dem"
	}
	filledDEMName <- if(missing(filledDEM)) "filledDEM" else filledDEM
	probsName <- if(missing(probs)) "problem_areas" else probs

	flowDirection <- "flow_direction"
	rgrass7::execGRASS("r.fill.dir", flags=c("overwrite", "quiet"), input=dem,
		output = filledDEMName, direction = flowDirection, areas = probsName)
	gs <- GSAppendRasterName(c(probsName, filledDEMName), gs)


	# if both missing, return raster stack; if only one return the other, if neither return session
	if(missing(filledDEM) & missing(probs)) {
		res <- GSGetRaster(c(filledDEMName, probsName), gs)
	} else if(missing(filledDEM)) {
		res <- GSGetRaster(filledDEMName, gs)
	} else if(missing(probs)) {
		res <- GSGetRaster(probsName, gs)
	} else {
		res <- gs
	}
	return(res)
}


#' Watershed analysis tools
#'
#' @param dem Raster or character; a digital elevation model, preferably one filled using [fillDEM()]. If specified as a character, a layer with that name from the existing [GrassSession()] given by gs will be used.
#' @param threshold Minimum size of an exterior watershed, in *cells*
#' @param gs An optional [GrassSession()]; if missing a new one will be created
#' @param outputName Optional character vector; if missing, a raster layer is returned, otherwise the raster is written to the grass session with this name and the modified GrassSession is returned.
#' @param file The file name of the raster to be returned (if no `outputName` is provided), see `details`.
#' @param ... Additional parameters to pass to [GrassSession()]
#' @details #' @details This is a wrapper for [r.watershed](https://grass.osgeo.org/grass74/manuals/r.watershed.html)
#' 
#' It is recommended to specify the `file` parameter (including the extension to specify
#' file format; e.g., .tif, .grd). If not specified, a temp file will be created and will be
#' lost at the end of the R session
#'
#' @return A [raster::stack()] with two layers, flow accumulation ('accumulation') and drainage
#' 		direction ('drainage') (if outputName is missing), or a GrassSession otherwise
watershed <- function(dem, threshold = 250, gs, accumulation, drainage, file = NULL, ...)
{
	if(missing(gs)) {
		gs <- GrassSession(dem, layerName = "dem", ...)
		dem <- "dem"
	} else if(!is.character(dem)) {
		gs <- GSAddRaster(dem, layerName = "dem", gs)
		dem <- "dem"
	}

	accumulationName <- if(missing(accumulation)) "accumulation" else accumulation
	drainageName <- if(missing(drainage)) "drainage" else drainage

	rgrass7::execGRASS("r.watershed", flags=c("overwrite", "quiet"), elevation = dem, 
		accumulation = accumulationName, drainage = drainageName)
	gs <- GSAppendRasterName(c(accumulationName, drainageName), gs)

	# if both missing, return raster stack; if only one return the other, if neither return session
	if(missing(accumulation) & missing(drainage)) {
		res <- GSGetRaster(c(accumulationName, drainageName), gs)
	} else if(missing(accumulation)) {
		res <- GSGetRaster(accumulationName, gs)
	} else if(missing(drainage)) {
		res <- GSGetRaster(drainageName, gs)
	} else {
		res <- gs
	}
	return(res)
}




#' Produce a raster showing which pixels are part of a stream
#' 
#' @param dem Raster or character; a digital elevation model, preferably one filled using [fillDEM()]. If specified as a character, a layer with that name from the existing [GrassSession()] given by gs will be used.
#' @param accumulation Raster or character; flow accumulation later; see details.
#' @param threshold Accumulation threshold; optional, see details
#' @param qthresh Quantile threshold to use if `threshold` is missing
#' @param weights A raster layer specifying the weights to apply to flow accumulation before comparing to `threshold`.
#' @param gs An optional [GrassSession()]; if missing a new one will be created
#' @param outputName Optional character; if missing, a raster layer is returned, otherwise the raster is written to the grass session with this name and the modified GrassSession is returned.
#' @param file The file name of the raster to be returned (if no `outputName` is provided), see `details`.
#' @param ... Additional parameters to pass to [GrassSession()] or r.stream.extract
#' @details A flow accumulation raster can be provided (as a rasterLayer or as the name of a raster in an existing Grass Session). If not provided, it will be computed.
#'
#' If no threshold is used, then one will be set automatically using the quantiles of the accumulation raster
#'
#' This is a wrapper for [r.stream.extract](https://grass.osgeo.org/grass74/manuals/r.stream.extract.html)
#' 
#' It is recommended to specify the `file` parameter (including the extension to specify
#' file format; e.g., .tif, .grd). If not specified, a temp file will be created and will be
#' lost at the end of the R session
#' @return R [raster::raster] with values set to a unique ID for pixels in a stream, NA otherwise (if outputName is missing) or a GrassSession otherwise
extractStream <- function(dem, accumulation, threshold, qthresh = 0.95, weights, gs, 
	outputName, file = NULL, ...)
{
	if(missing(gs)) {
		gs <- GrassSession(dem, layerName = "dem", ...)
		dem <- "dem"
	} else if(!is.character(dem)) {
		gs <- GSAddRaster(dem, layerName = "dem", gs)
		dem <- "dem"
	}

	if(!is.character(accumulation)) {
		if(!missing(weights))
			accumulation <- accumulation * weights
		gs <- GSAddRaster(accumulation, layerName = "accumulation", gs)
		if(missing(threshold))
			accTh <- raster::values(accumulation)
		accumulation <- "accumulation"
	} else if(!missing(weights)) {
		accu <- GSGetRaster(accumulation, gs)
		accu <- accu * weights
		accumulation <- paste0(accumulation, "_weighted")
		gs <- GSAddRaster(accu, accumulation, gs)
	}

	if(missing(threshold)) {
		if(!exists("accTh"))
			accTh <- raster::values(GSGetRaster(accumulation, gs))
		accTh[accTh < 0] <- NA
		threshold <- quantile(accTh, qthresh, na.rm=TRUE)
	}

	streamRaster <- if(missing(outputName)) "streams" else outputName

	rgrass7::execGRASS("r.stream.extract", flags=c("overwrite", "quiet"), elevation=dem, 
		accumulation = accumulation, threshold = threshold, stream_raster = streamRaster)
	gs <- GSAppendRasterName(streamRaster, gs)

	res <- if(missing(outputName)) GSGetRaster(streamRaster, gs) else gs
	return(res)
}