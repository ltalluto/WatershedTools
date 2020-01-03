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
	## drainage will be added later, after it is fixed by WSConnectivity
	dataRasters <- list()
	if(!missing(elevation)) dataRasters$elevation <- elevation
	if(!missing(accumulation)) dataRasters$accumulation <- accumulation
	if(!missing(catchmentArea)) dataRasters$catchmentArea <- catchmentArea
	if(!missing(otherLayers)) dataRasters$otherLayers <- otherLayers
	layerStack <- lapply(dataRasters, function(x) {
		if(!raster::compareRaster(stream, x, stopiffalse = FALSE))
			x <- raster::crop(x, stream)
		raster::mask(x, stream)
	})

	## create pixel IDs and add other layers, if present
	allRasters <- raster::stack(stream, stream)
	names(allRasters) <- c('reachID', 'id')
	if(length(layerStack) > 0) {
		layerStack <- raster::stack(layerStack)
		allRasters <- raster::addLayer(allRasters, layerStack)
	}
	maskIndices <- which(!is.na(raster::values(stream)))
	allRasters$id[maskIndices] <- 1:length(maskIndices)

	allSPDF <- rasterToSPDF(allRasters, complete.cases = TRUE)
	if(!raster::compareRaster(allRasters, drainage, stopiffalse = FALSE))
		drainage <- raster::crop(drainage, allRasters)
	adjacency <- WSConnectivity(drainage, allRasters$id)
	allSPDF <- sp::merge(allSPDF, adjacency$drainage, by = 'id', all.x = TRUE)
	allSPDF$length <- WSComputeLength(allSPDF$drainage, raster::res(drainage))
	allSPDF$reachID <- renumberReaches(allSPDF$reachID)

	wsobj <- list(data = allSPDF, adjacency = adjacency$adjacency)
	class(wsobj) <- c("Watershed", class(wsobj))
	wsobj$data$reachID <- renumberReaches(wsobj$data$reachID)
	wsobj$reach_adjacency <- reachByReachAdj(wsobj)
	wsobj$reach_connectivity <- reachByReachConn(wsobj, self = FALSE)
	attr(wsobj, "version") <- packageVersion("WatershedTools")
	return(wsobj)
}


#' Compute connectivity matrix
#'
#' @param drainage Drainage direction raster
#' @param stream Stream network raster; see `details`
#'
#' @details The stream network raster should be NA in all cells not considered a part of the
#'		river network. The pixel values of the raster must be unique IDs representing individual
#'		stream reaches to model. At present, the only supported reach size is a single pixel, thus
#'		each pixel must have a unique value.
#' @return A list with two elements, the first containing corrected drainage directions, the 
#' 		second with A [Matrix::sparseMatrix()] representation of the river network. For a `stream` 
#' 		input raster with `n` non-NA cells, the dimensions of this matrix will be n by n. Dimnames
#'		of the matrix will be the pixel IDs from the `stream` input raster. Values of the
#'		matrix cells are either 0 or 1; a zero indicates no flow, a one in cell i,j indicates
#'		that reach `i` receives water from reach `j`.
#' @keywords internal
WSConnectivity <- function(drainage, stream) {
	ids <- raster::values(stream)
	inds <- which(!is.na(ids))
	ids <- ids[inds]
	if(any(duplicated(ids)))
		stop("Stream IDs must be all unique")

	rowMat <- matrix(1:raster::nrow(drainage), nrow=raster::nrow(drainage), 
		ncol=raster::ncol(drainage))
	colMat <- matrix(1:raster::ncol(drainage), nrow=raster::nrow(drainage), 
		ncol=raster::ncol(drainage), byrow=TRUE)
	coordRas <- raster::stack(list(x = raster::raster(colMat, template = drainage), 
		y = raster::raster(rowMat, template = drainage), drainage = drainage, id = stream))
	coordRas <- raster::mask(coordRas, stream)

	xy <- WSFlowTo(coordRas[inds])
	res <- xy[,c('fromID', 'drainage')]
	colnames(res)[1] <- 'id'
	list(drainage = res, adjacency = Matrix::sparseMatrix(xy[,'toID'], xy[,'fromID'], 
		dims=rep(length(inds), 2), dimnames = list(ids, ids), x = 1))
}




#' Compute which pixels flow into which other pixels
#' @param mat A matrix with minimum three columns, the first being the x-coordinate, second the y,
#'		and third the ID.
#' @return A matrix of IDs, the first column the source, the second column the destination
#' @keywords internal
WSFlowTo <- function(mat) {
	newy <- mat[,2]
	newx <- mat[,1]

	ind <- which(mat[,3] > 0)
	xoffset <- c(1, 0, -1, -1, -1, 0, 1, 1)
	yoffset <- c(-1, -1, -1, 0, 1, 1, 1, 0)
	newx[ind] <- newx[ind] + xoffset[mat[ind,3]]
	newy[ind] <- newy[ind] + yoffset[mat[ind,3]]
	na_ind <- which(mat[,3] < 0 | newx < 1 | newy < 1 | newx > max(mat[,1]) | newy > max(mat[,2]))
	newx[na_ind] <- newy[na_ind] <- mat[na_ind, 'drainage'] <- NA
	resMat <- cbind(mat, newx, newy)
	resMat <- merge(resMat[,c('newx', 'newy', 'id', 'drainage')], resMat[,c('x', 'y', 'id')], by = c(1,2), all.x = TRUE)
	resMat <- resMat[,c(1,2,4,3,5)]
	colnames(resMat)[4:5] <- c('fromID', 'toID')

	resMat <- WSCheckDrainage(resMat, mat)
	resMat <- resMat[order(resMat[,'fromID']),]
	resMat <- resMat[complete.cases(resMat),]
	return(resMat)
}

#' Check and fix problems with drainage direction
#' @param connMat preliminary connectivity matrix
#' @param drainMat Drainage direction matrix
#' @param prevProbs Previous number of problems, to allow stopping if no improvement on 
#' 		subsequent calls
#' @details In some cases, drainage direction rasters don't agree with flow accumulation, resulting
#' 		in a delineated stream that doesn't have the right drainage direction. This function
#' 		attempts to detect and fix this in the adjacency matrix and the drainage layer
#' @keywords internal
#' @return A corrected connectivity matrix
WSCheckDrainage <- function(connMat, drainMat, prevProbs = NA) {
	probs <- which(is.na(connMat[,'toID']) & connMat[,'drainage'] > 0)
	if(length(probs) == 0 | (!is.na(prevProbs) & length(probs) == prevProbs))
		return(connMat)

	prFix <- do.call(rbind, 
		lapply(connMat[probs,'fromID'], WSFixDrainage, drainMat = drainMat, connMat = connMat))
	connMat <- connMat[-probs,]
	connMat <- rbind(connMat, prFix)
	WSCheckDrainage(connMat, drainMat, prevProbs = length(probs))
}


#' Fix drainage direction for a single pixel
#' @param id ID of the problematic pixel
#' @param connMat preliminary connectivity matrix
#' @param drainMat Drainage direction matrix
#' @keywords internal
WSFixDrainage <- function(id, drainMat, connMat) {
	i <- which(drainMat[,'id'] == id) # problem cell index
	j <- which(connMat[,'toID'] == id) # upstream of problem cell
	upID <- connMat[j,'fromID']
	x <- drainMat[i,'x']
	y <- drainMat[i,'y']
	downInd <- which(drainMat[,'x'] >= x-1 & drainMat[,'x'] <= x+1 & drainMat[,'y'] >= y-1 & drainMat[,'y'] <= y+1 & !(drainMat[,'id'] %in% c(id, upID)))
	out <- connMat[connMat[,'fromID'] == id,]
	if(length(downInd) == 1) {
		out[,'newx'] <- drainMat[downInd,'x']
		out[,'newy'] <- drainMat[downInd,'y']
		out[,'toID'] <- drainMat[downInd,'id']
		xo <- out[,'newx'] - drainMat[i,'x']
		yo <- out[,'newy'] - drainMat[i,'y']
		xoffset <- c(1, 0, -1, -1, -1, 0, 1, 1)
		yoffset <- c(-1, -1, -1, 0, 1, 1, 1, 0)
		out[,'drainage'] <- which(xoffset == xo & yoffset == yo)
	}
	out
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
#' @return a `data.frame` containing data for all outlets
#' @export
outlets <- function(ws) {
	as.data.frame(ws$data[Matrix::colSums(ws$adjacency) == 0,])
}


#' Get the pixel ID of the next downstream pixel for each pixel in the watershed
#' @param ws A watershed object
#' @return A vector of pixel IDs
#' @export
downstreamPixelIds <- function(ws) {
	mat <- Matrix::which(ws$adjacency == 1, arr.ind = TRUE)
	endpt <- Matrix::which(Matrix::colSums(ws$adjacency) == 0)
	mat <- rbind(mat, c(NA, endpt))
	# rearrange so the UPSTREAM pixels (second column) indicate the row number
	mat <- mat[order(mat[,2]),]
	if(!all(mat[,2] == 1:nrow(mat)))
		stop("There is an error with the topology")
	mat[,1]
}


#' Splits a reachID at selected points
#' @param ws A watershed object
#' @param points A vector of ids at which to split reaches
#' @param na_ignore Logical, if `TRUE` NAs will be dropped from the data. If `FALSE`, any NAs
#' 		in `points` will cause an error.
#' @return A modified watershed with new reachIDs
#' @export
splitReaches <- function(ws, points, na_ignore = FALSE) {
	if(na_ignore & any(is.na(points))) {
		points <- points[!is.na(points)]
	} else if(any(is.na(points))) {
		stop("NA in points; use na_ignore = TRUE to ignore")
	}
	for(pt in points) {
		ids <- which(ws$data$reachID == ws$data$reachID[pt])
		reachAdj <- ws$adjacency[ids, ids]
		mostUpstream <- ids[which(Matrix::rowSums(reachAdj) == 0)]
		if(pt != mostUpstream) {
			chIds <- connect(ws, mostUpstream, pt)
			ws$data$reachID[chIds] <- max(ws$data$reachID) + 1
		}
	}
	
	# re-create topology
	ws$data$reachID <- renumberReaches(ws$data$reachID)
	ws$reach_adjacency <- reachByReachAdj(ws)
	ws$reach_connectivity <- reachByReachConn(ws, self = FALSE)

	ws
}


#' Produce a reach by reach connectivity matrix
#' 
#' @param ws A watershed
#' @param self If TRUE, a reach is considered connected to itself
#' Nonzero values indicate that a reachID in a row is downstream from that column
#' @return A [Matrix::sparseMatrix()]
#' @keywords internal
reachByReachConn <- function(ws, self = TRUE) {
	adj <- ws$reach_adjacency
	if(self)
		diag(adj) <- 1
	for(i in 1:nrow(adj))
		adj <- adj %*% adj + adj
	adj[adj != 0] <- 1
	adj
}


#' Produce a reach by reach adjacency matrix
#' @keywords internal
reachByReachAdj <- function(ws) {
	if(max(ws$data$reachID) != length(unique(ws$data$reachID)))
		stop("ReachIDs not in 1:length(unique(reachID)); they must be renumbered")
	
	reaches <- sort(unique(ws$data$reachID))
	adjMat <- parallel::mclapply(reaches, function(r) reachAdj(ws, r))
	adjMat <- do.call(rbind, adjMat)
	adjspMat <- Matrix::sparseMatrix(adjMat[,1], adjMat[,2], dims=rep(max(ws$data$reachID), 2),
		dimnames = list(reaches, reaches))
	adjspMat
}

#' Find all reaches directly upstream from a given reach
#' @keywords internal
reachAdj <- function(ws, rch) {
	ids <- which(ws$data$reachID == rch)
	reachAdj <- ws$adjacency[ids, ids, drop = FALSE]
	mostUpstream <- ids[which(Matrix::rowSums(reachAdj) == 0)]
	upRch <- which(ws$adjacency[mostUpstream,] == 1)
	if(length(upRch) > 0) {
		upRch <- ws$data$reachID[upRch]
		return(cbind(rch, upRch))
	} else return(NULL)
}


#' Renumber reachIDs to start at 1 and increase one at a time, preserving ordering
#' @keywords internal
renumberReaches <- function(rIDs) {
	newNums <- 1:length(unique(rIDs))
	oldNums <- sort(unique(rIDs))
	newNums[match(rIDs, oldNums)]
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

