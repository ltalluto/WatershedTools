

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



