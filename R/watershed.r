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
	dataRasters <- list(drainage = drainage)
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
	# print(head(allSPDF))
	# allSPDF <- allSPDF[!is.na(allSPDF$id),]
	wsobj <- list(data = allSPDF)
	wsobj$adjacency <- WSConnectivity(allRasters$drainage, allRasters$id)

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
#' @return A [Matrix::sparseMatrix()] representation of the river network. For a `stream` input
#'		raster with `n` non-NA cells, the dimensions of this matrix will be n by n. Row/column
#'		names of the matrix will be the pixel IDs from the `stream` input raster. Values of the
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
	Matrix::sparseMatrix(xy[,2], xy[,1], dims=rep(length(inds), 2), 
		dimnames = list(ids, ids), x = 1)
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
	yoffset <- c(-1, -1, -1, 0, 1, 1, 1, 1)
	newx[ind] <- newx[ind] + xoffset[mat[ind,3]]
	newy[ind] <- newy[ind] + yoffset[mat[ind,3]]
	na_ind <- which(mat[,3] < 0 | newx < 1 | newy < 1 | newx > max(mat[,1]) | newy > max(mat[,2]))
	newx[na_ind] <- newy[na_ind] <- NA
	resMat <- cbind(mat, newx, newy)
	resMat <- merge(resMat[,c('newx', 'newy', 'id')], resMat[,c('x', 'y', 'id')], by = c(1,2), all.x = TRUE)
	resMat <- resMat[,c('id.x', 'id.y')]
	colnames(resMat) <- c('fromID', 'toID')
	resMat <- resMat[order(resMat[,'fromID']),]
	resMat <- resMat[complete.cases(resMat),]
	return(resMat)
}