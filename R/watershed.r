#' Creates a Watershed object
#' @export
Watershed <- function(stream, drainage, ...) {
	wsobj <- list()

	## need to compute pixel IDs
	wsobj$adjacency <- WSConnectivity(drainage, stream)

	## maybe don't do this; instead store attributes in a spatialGridDataFrame
	attr(wsobj, "projection") <- sp::proj4string(stream)
	attr(wsobj, "resolution") <- raster::res(stream)

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
	adjacency <- Matrix::sparseMatrix(xy[,2], xy[,1], dims=rep(length(inds), 2), 
		dimnames = list(ids, ids), x = 1)
}

#' Run the transport model
#' @param adjacency The adjacency matrix describing the stream network
#' @param initial A vector of initial conditions
#' @param steps How may time steps to run the model
#' @param weights Optional weighting to apply to the adjacency matrix
#' @export
transport <- function(adjacency, initial, steps, weights) {
	if(!missing(weights))
		adjacency <- t(t(adjacency) * weights)

	state <- matrix(NA, nrow=length(initial), ncol = steps+1)
	state[,1] <- initial
	for(i in 2:(steps+1)) {
		state[,i] <- (adjacency %*% state[,i-1])[,1]
	}
	return(state)
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
	resMat <- merge(test, test[,c('x', 'y', 'id')], by.x = c(6,7), by.y = c(1,2), all.x = TRUE)
	resMat <- resMat[,c('id.x', 'id.y')]
	colnames(resMat) <- c('fromID', 'toID')
	resMat <- resMat[order(resMat[,'fromID']),]
	resMat <- resMat[complete.cases(resMat),]
	return(resMat)
}