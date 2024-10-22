#' Creates a Watershed object
#' @param stream Stream network raster, required
#' @param drainage Drainage direction raster, required
#' @param elevation Optional elevation raster
#' @param accumulation Optional flow accumulation raster
#' @param catchmentArea Optional catchment area raster
#' @param otherLayers RasterStack of other data layers to add to the Watershed object
#' @details All raster maps will be cropped to the stream network. The values in `stream` will
#' 		be automatically assigned to a reachID field in the Watershed object.
#' @return A watershed object
#' @export
Watershed <- function(stream, drainage, elevation, accumulation, catchmentArea, otherLayers) {
	dataRasters =list()
	if(!missing(drainage)) dataRasters$drainage = drainage
	if(!missing(elevation)) dataRasters$elevation = elevation
	if(!missing(accumulation)) dataRasters$accumulation = accumulation
	if(!missing(catchmentArea)) dataRasters$catchmentArea = catchmentArea
	if(!missing(otherLayers)) dataRasters$otherLayers = otherLayers
	layerStack =lapply(dataRasters, function(x) {
		if(!raster::compareRaster(stream, x, stopiffalse = FALSE))
			x =raster::crop(x, stream)
		raster::mask(x, stream)
	})

	## create pixel IDs and add other layers, if present
	allRasters = raster::stack(stream, stream)
	names(allRasters) = c('reachID', 'id')
	if(length(layerStack) > 0) {
		layerStack = raster::stack(layerStack)
		allRasters = raster::addLayer(allRasters, layerStack)
	}
	maskIndices = which(!is.na(raster::values(stream)))
	allRasters$id[maskIndices] = 1:length(maskIndices)

	allSPDF = as(allRasters, "SpatialPixelsDataFrame")
	allSPDF = allSPDF[complete.cases(allSPDF@data),]
	adjacency = Matrix::t(watershed::pixel_topology(drainage = allRasters$drainage, 
		stream = allRasters$reachID, id = allRasters$id))
	allSPDF$length = WSComputeLength(allSPDF$drainage, raster::res(drainage))
	allSPDF$vReachNumber <- allSPDF$reachID

	wsobj <- list(data = allSPDF, adjacency = adjacency)
	class(wsobj) <- c("Watershed", class(wsobj))
	
	wsobj = .rebuild_reaches(wsobj)
		
	attr(wsobj, "version") <- packageVersion("WatershedTools")
	return(wsobj)
}



#' Convert a Watershed back to a watershed-style stack
#' @param ws A watershed
#' @return A raster stack following the format from the watershed package
#' @keywords internal
.make_ws_stack = function(ws) {
	crs = sp::proj4string(ws$data)
	co = coordinates(ws$data)
	val = raster::stack(raster::rasterFromXYZ(cbind(coordinates(ws$data), 
		ws$data$accumulation, ws$data$drainage, ws$data$reachID, ws$data$catchmentArea, ws$data$id)))
	names(val) = c('accum', 'drainage', 'stream', 'catchment', 'id')
	val
}






#' Compute length to next pixel given drainage direction
#' @param drainage drainage direction vector
#' @param cellsize size of each cell (vector of length 2)
#' @keywords internal
#' @return vector of lengths
WSComputeLength <- function(drainage, cellsize) {
	cellLength <- rep(cellsize[1], length(drainage))
	if(abs(cellsize[1] - cellsize[2]) > 1e-4) {
		vertical <- which(drainage %in% c(2,6))
		cellLength[vertical] <- cellsize[2]
	}
	diagonal <- which(drainage %in% c(1,3,5,7))
	cellLength[diagonal] <- sqrt(cellsize[1]^2 + cellsize[2]^2)
	cellLength
}


#' Get data from all confluences of a watershed
#' 
#' @param ws Watershed object
#' @return A `data.frame` containing data for all confluences
#' @export
confluences <- function(ws) {
	as.data.frame(ws$data[Matrix::rowSums(ws$adjacency) > 1,])
}

#' Get data from all headwaters of a watershed
#' 
#' @param ws Watershed object
#' @return a `data.frame` containing data for all headwaters
#' @export
headwaters <- function(ws) {
	as.data.frame(ws$data[Matrix::rowSums(ws$adjacency) == 0,])
}

#' Get data from all outlets of a watershed
#' 
#' @param ws Watershed object
#' @param rid vector of reach IDs, if missing returns outlet for entire network
#' @param output Output type to return
#' @return a `data.frame` or a `SpatialPixelsDataFrame` containing data for all outlets 
#' @export
outlets <- function(ws, rid, output = c("data.frame", "Spatial")) {
	output = match.arg(output)
	if(!missing(rid)) {
		out_ind = sapply(rid, function(i) {
			ii = which(ws$data$reachID == i)
			mat = ws$adjacency[ii,ii, drop=FALSE]
			pix = which(Matrix::colSums(mat) == 0)
			as.integer(rownames(mat)[pix])
		})
	} else {
		out_ind = which(Matrix::colSums(ws$adjacency) == 0)
	}
	res = ws$data[out_ind,]
	if(output == "data.frame") {
		res = as.data.frame(res)
	}
	res
}


#' Get the pixel ID of the next downstream pixel for each pixel in the watershed
#' @param ws A watershed object
#' @return A vector of pixel IDs
#' @export
downstreamPixelIds <- function(ws) {
	mat <- Matrix::which(ws$adjacency != 0, arr.ind = TRUE)
	endpt <- Matrix::which(Matrix::colSums(ws$adjacency) == 0)
	mat <- rbind(mat, c(NA, endpt))
	# rearrange so the UPSTREAM pixels (second column) indicate the row number
	mat <- mat[order(mat[,2]),]
	if(!all(mat[,2] == 1:nrow(mat)))
		stop("There is an error with the topology")
	mat[,1]
}




#' Extract pixelIDs from a watershed at spatial locations
#' @param ws A watershed object
#' @param x An object inheriting from `sp::SpatialPoints()`
#' @return A vector of pixel IDs
#' @export
extract <- function(ws, x) {
	ras <- raster::rasterFromXYZ(ws[,c('x', 'y', 'id')])
	raster::extract(ras, x)
}



#' Compute a site by pixel accumulation matrix
#' 
#' The default behavior computes distance, where positive numbers indicate downstream
#' distances and negative numbers indicate upstream distances. Other variables can also
#' be used, but in all cases the values will be summed to compute the 'distance'
#' 
#' Upstream distances do NOT include intermediate pixels; they only include pixels in `x`
#' 
#' @param ws A Watershed
#' @param x A vector of pixel ids from which to compute the distance 
#' @param variable The variable to use for the distance
#' @return A matrix with dimensions `length(x)` by `nrow(ws)`
#' @export
siteByPixel <- function(ws, x, variable = 'length') {
	dsPixes <- downstreamPixelIds(ws)
	dm <- dmat(x, dsPixes, nrow(ws$data), ws[,variable])
	rownames(dm) <- x
	colnames(dm) <- ws[,'id']
	dm
}

