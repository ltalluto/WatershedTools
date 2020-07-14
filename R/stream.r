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
#' @export
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
#' 		direction ('drainage') (if outputName is missing), or a GrassSession otherwise.
#' 		Drainage direction values can be negative (indicating water flowing out of the map region),
#' 		zero (indicating a portion of a natural depression), or positive. Positive values are 
#' 		measured counterclockwise from northeast; 1 is northeast, 2 north,, etc through 8 flowing #'	  due east.
#' @export
drainageAccumulation <- function(dem, threshold = 250, gs, accumulation, drainage, 
	file = NULL, ...)
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




#' Produce a map of a stream network
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
#' @return If `outputName` is missing and `type=='raster'`; a [raster::raster]; if 
#'   `outputName` is missing and `type=='vector'`; a [sp::SpatialLinesDataFrame];
#'   otherwuse a `GrassSession`
#' @export
extractStream <- function(dem, accumulation, threshold, qthresh = 0.95, weights, gs, 
	outputName, type = c('raster', 'vector', 'both'), file = NULL, ...)
{
	type <- match.arg(type)

	if((is.character(dem) | is.character(accumulation)) & missing(gs))
		stop("If gs is missing, dem and accumulation must be a Raster or SpatialPixelsDataFrame")
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

	streamName <- if(missing(outputName)) "streams" else outputName

	rgrass7::execGRASS("r.stream.extract", flags=c("overwrite", "quiet"), elevation=dem, 
		accumulation = accumulation, threshold = threshold, stream_raster = streamName, stream_vector = streamName)
	gs <- GSAppendRasterName(streamName, gs = gs)

	if(!missing(outputName)) {
		res <- gs
	} else if(type == 'raster') {
		res <- GSGetRaster(streamName, gs, file = file)
	} else if(type == 'vector') {
		res <- rgrass7::readVECT(streamName, ignore.stderr=TRUE, type="line")
	} else {
		res <- list(raster = GSGetRaster(streamName, gs, file = file), 
			vector = rgrass7::readVECT(streamName, ignore.stderr=TRUE, type="line"))
	}

	gs <- GSClean(streamName, gs, "vector")
	return(res)
}

#' Snap points to a stream network raster
#' 
#' @param x SpatialPoints or similar object
#' @param stream Raster or character; a stream network raster (e.g., from [extractStream()]). 
#'   If specified as a character, a layer with that name from the existing [GrassSession()] 
#'   given by gs will be used.
#' @param buff The distance (in meters if x is in geographic coordinates, map units otherwise) to restrict the search from points in x
#' @details If buff is too small, a nearby stream may not be found, in which case the original coordinates are returned
#' @return A SpatialPointsDataFrame, with the attributes from `x` and new coordinates
#' @export
snapToStream <- function(x, stream, buff)
{
	if(grepl("longlat", sp::proj4string(x))) {
		warning("x has geographic coordinates; reprojecting to epsg:32632")
		projOrig <- sp::proj4string(x)
		x <- sp::spTransform(x, sp::CRS("+init=epsg:32632"))
		stream <- raster::projectRaster(stream, crs = sp::proj4string(x))
	}

	newCoords <- t(sapply(1:length(x), function(i) findClosest(x[i,], stream, buff = buff)))
	result <- raster::as.data.frame(x)
	result <- cbind(newCoords, result)
	sp::coordinates(result) <- c(1,2)
	sp::proj4string(result) <- sp::proj4string(x)
	if(exists("projOrig"))
		result <- sp::spTransform(result, projOrig)
	return(result)
}

#' Find the closest non NA point in a raster
#' @param x a spatial object
#' @param y a raster
#' @param buff The distance (in map units) to restrict the search from points in x
#' @details Finds the points in y that are closest to x
#' @return The index in y of the closest point
#' @keywords internal
findClosest <- function(x, y, buff)
{
	if(length(x) > 1)
		stop("findClosest is not vectorised at this time")
	smRas <- CropPtBuff(x, y, buff)
	xc <- sp::coordinates(x)
	inds <- which(!is.na(raster::values(smRas)))
	if(length(inds) == 0)
		return(xc)
	yc <- sp::coordinates(smRas)[inds,,drop=FALSE]
	xy <- rbind(xc,yc)
	dists <- dist(xy)[1:nrow(yc)]
	yc[which.min(dists),,drop=FALSE]
}


#' Crop raster around a point
#' @param pt a spatial point
#' @param ras a raster
#' @param buff distance to crop around point, in map units
#' @return Cropped raster
#' @keywords internal
CropPtBuff <- function(pt, ras, buff)
{
	xc <- sp::coordinates(pt)
	newextent <- raster::extent(xc[1] - buff, xc[1] + buff, xc[2] - buff, xc[2] + buff)
	raster::crop(ras, newextent)
}




#' Compute strahler order for pixel x
#' @param ws A watershed
#' @param x A pixel id
#' @param streamOrder Vector of streamOrders already computed (should be NA if not done yet)
#' @return A matrix, giving pixel IDs in the first column and computed stream order in the second
#' @keywords internal
doStrahler <- function(ws, x, streamOrder) {
	upIDs <- which(ws$adjacency[x,] == 1)
	if(length(upIDs) < 1) {
		val <- 1
	} else {
		upOrders <- streamOrder[upIDs]
		if(any(is.na(upOrders))) {
			val <- NA
		} else if(sum(upOrders == max(upOrders)) == 1) {
			val <- max(upOrders)
		}  else
			val <- max(upOrders) + 1
	}

	pixes <- ws$data$id[ws$data$reachID == ws$data$reachID[x]]
	cbind(pixes, rep(val, length(pixes)))
}

#' Compute strahler stream order
#' 
#' @param ws A watershed
#' @param parallel If TRUE, will compute in parallel to speed things up. Defaults to TRUE
#' if parallel package is installed
#' @details Computes Strahler stream order
#' @return a vector, with one element per pixel in the watershed.
strahler <- function(ws, parallel = TRUE) {
	if(parallel && requireNamespace("parallel")) {
		FUN <- parallel::mclapply
	} else
		FUN <- lapply
	streamOrder <- rep(NA, length=nrow(ws$data))
	pix <- c(headwaters(ws)$id, confluences(ws)$id)
	while(any(is.na(streamOrder))) {
		ord <- FUN(pix, function(x) doStrahler(ws, x, streamOrder))
		ord <- do.call(rbind, ord)
		streamOrder[ord[,1]] <- ord[,2]
		prev <- pix
		pix <- pix[is.na(streamOrder[pix])]
		if(length(pix) == length(prev))
			stop("No progress being made, check topology")
	}
	streamOrder
}



