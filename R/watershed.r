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
	dataRasters <- list()
	if(!missing(elevation)) dataRasters$elevation <- elevation
	if(!missing(accumulation)) dataRasters$accumulation <- accumulation
	if(!missing(catchmentArea)) dataRasters$catchmentArea <- catchmentArea
	if(!missing(otherLayers)) dataRasters$otherLayers <- otherLayers
	allRasters <- lapply(dataRasters, function(x) {
		if(!raster::compareRaster(stream, x, stopiffalse = FALSE))
			x <- raster::crop(x, stream)
		raster::mask(x, stream)
	})
	allRasters <- raster::stack(allRasters)

	## create pixel IDs
	allRasters$reachID <- allRasters$id <- stream
	maskIndices <- which(!is.na(raster::values(stream)))
	allRasters$id[maskIndices] <- 1:length(maskIndices)

	allSPDF <- rasterToSPDF(allRasters, complete.cases = TRUE)
	adjacency <- WSConnectivity(drainage, allRasters$id)
	allSPDF <- merge(allSPDF, adjacency$drainage, by = 'id', all.x = TRUE)
	allSPDF$length <- WSComputeLength(allSPDF$drainage, res(drainage))

	wsobj <- list(data = allSPDF, adjacency = adjacency$adjacency)
	class(wsobj) <- c("Watershed", class(wsobj))
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
