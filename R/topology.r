#' Compute the slope along the river channel
#' 
#' This computes the change in elevation from one watershed unit to the next. Computation by
#' pixel is only recommended if the elevation model used to create the watershed is very
#' accurate; otherwise unrealistic slopes can occur. Reach-based computation uses a much
#' larger range of pixels to compute the slope, and so is less prone to DEM errors.
#' 
#' @param ws A watershed
#' @param by Should the slope be computed for each reach (recommended) or from one pixel to the 
#' next
#' 
#' @return A vector of slopes, one per pixel
#' @examples
#' \donttest{
#' ws_slope(ws)
#' ws_slope(ws, by = "pixel")
#' }
#' @export
ws_slope = function(ws, by = c("reach", "pixel")) {
	by = match.arg(by)
	if(by == "pixel") {
		adj = methods::as(ws$adjacency, "dgTMatrix")
		adj = cbind(adj@i + 1, adj@j + 1, adj@x) ## Matrix is 0-indexed

		## this function assumes the adjacenty matrix is unweighted
		adj[,3] = ws$data$length[adj[,2]]
		delev = ws$data$elevation[adj[,1]] - ws$data$elevation[adj[,2]]

		slope = numeric(nrow(ws$data)) * NA
		slope[adj[,1]] = delev/adj[,3]
		slope[adj[,2]] = delev/adj[,3]
		slope = abs(slope)
	} else {
		delev = tapply(ws$data$elevation, ws$data$reachID, function(x) max(x) - min(x))
		length = tapply(ws$data$length, ws$data$reachID, sum)
		slope = delev / length
		slope = slope[match(ws$data$reachID, names(slope))]		
	}
	return(slope)
}

#' Nearest neighbors of pixels
#' @keywords internal
ds = function(pix, ws) {
	val = which(ws$adjacency[,pix] == 1)
	if(length(val) == 0)
		val = NA
	val
}

#' Nearest neighbors of pixels
#' @keywords internal
us = function(pix, ws) {
	val = which(ws$adjacency[pix,] == 1)
	if(length(val) == 0)
		val = NA
	val
}

#' Produces a site by pixel distance matrix
#'
#' Note that pixels in the watershed that are not connected to any of the sites will be excluded
#' from the result.
#' @param ws A watershed
#' @param x A vector of sites
#' @param variable The variable to use as the distance
#' @return A site by pixel distance matrix
#' @export
wsDistance <- function(ws, x, variable = 'length') {
	downs <- parallel::mclapply(x, function(xx) {
		res <- accumulate(ws, xx, parallel=FALSE, variable = variable)
		if(!is.matrix(res)) res <- matrix(res, nrow=1)
		res
	})
	ups <- parallel::mclapply(x, function(xx) {
		res <- accumulate(ws, upstream=Inf, downstream = xx, parallel = FALSE, 
			direction = 'up', variable = variable)
		if(!is.matrix(res)) res <- matrix(res, nrow=1)
		res
	})

	matTall <- do.call(rbind.data.frame, mapply(function(xx, ds, us) {
		dsMat <- cbind(row=rep(xx, nrow(ds)), col = ds[,1], val = ds[,2])
		usMat <- cbind(row=rep(xx, nrow(us)), col = us[,1], val = us[,2])
		rbind(dsMat, usMat)
	}, xx = x, ds = downs, us = ups, SIMPLIFY = FALSE))
	matWide <- reshape2::acast(matTall, row ~ col, value.var = 'val', 
		fun.aggregate = function(x) ifelse(length(x) == 0, NA, x[1]), fill=NA_real_, drop=FALSE)

	matWide[order(as.integer(rownames(matWide))), order(as.integer(colnames(matWide)))]
}


#' Returns all points in the watershed connecting two points
#' 
#' This is a simple wrapper to accumulate that only returns the pixels
#' @param ws Watershed object
#' @param upstream ID number of upstream point; if Inf ALL upstream pixels will be returned
#' @param downstream ID of downstream point; if Inf ALL downstream pixels will be returned
#' 
#' @return A vector of pixel ids
#' @export
connect <- function(ws, upstream, downstream = Inf) {
	direction <- "down"
	if(is.infinite(upstream))
		direction <- "up"
	accumulate(ws, upstream, downstream, direction)[,1]
}


#' Returns all points in the watershed connecting two points, accumulating a value between
#' 
#' If direction is down, then the accumulated distances will reflect the distance from the 
#' `upstream` point to each downstream point, and values will be positive. If direction is "up",
#' then the distances will be negative and accumulate from `downstream` to `upstream`.
#' @param ws Watershed object
#' @param upstream ID number of upstream point; if Inf ALL upstream pixels will be returned
#' @param downstream ID of downstream point; if Inf ALL downstream pixels will be returned
#' @param direction In which direction should the variable be accumulated; see 'details'
#' @param variable The variable to accumulate
#' @param parallel Boolean, should parallel processing be used?
#' @return A matrix, first column connected pixelIDs, second the accumulated variable
#' @export
accumulate <- function(ws, upstream, downstream = Inf, direction = c("down", "up"), 
			variable = 'length', parallel = FALSE) {
	direction <- match.arg(direction)
	if(parallel && !requireNamespace("parallel"))
		parallel <- FALSE
	dsPixes <- downstreamPixelIds(ws)
	if(length(downstream) == 1 && is.infinite(downstream))
		downstream <- outlets(ws)$id
	if(length(downstream) == 1) {
		if(length(upstream) == 1 && is.infinite(upstream)) {
			rid <- ws[downstream, 'reachID']
			tops <- headwaters(ws)[,'id']
			upReaches <- Matrix::which(ws$reach_connectivity[rid,] == 1)
			if(rid %in% ws[tops, 'reachID'])
				upReaches <- c(upReaches, rid)
			upstream <- tops[ws[tops, "reachID"] %in% upReaches]
		}
		if(length(upstream) == 1) {
			accum <- do.call(cbind, connectCPP(dsPixes, upstream, downstream, ws[, variable]))
			if(direction == 'up') 
				accum[,2] <- accum[,2] - accum[accum[,1] == downstream,2]
		} else {
			if(parallel) {
				accum <- do.call(rbind, 
					parallel::mclapply(upstream, function(x) 
						accumulate(ws, x, downstream, direction, variable)))
			} else {
				accum <- do.call(rbind, lapply(upstream, function(x) 
						accumulate(ws, x, downstream, direction, variable)))
			}
		}
	} else {
		accum <- do.call(rbind, 
				lapply(downstream, function(x) accumulate(ws, upstream, x, direction, 
					variable, parallel)))
	}
	## remove redundancies; it is quite common to traverse the same pixel multiple times
	accum <- accum[!duplicated(accum[,1]),, drop=FALSE]
	return(accum)
}


#' Produce a site by reach connectivity matrix
#' 
#' A value of 1 at indices `[i,j]` indicates that reach `i` is connected to (i.e., downstream of)
#' reach `j`.
#' @param ws A watershed object
#' @param points A vector of pixel id numbers for sites
#' @param names An optional vector of site names
#' @param self If TRUE, a reach is considered connected to itself
#' @return A [Matrix::sparseMatrix()]
#' @export
siteByReach <- function(ws, points, names, self = TRUE) {
	rIDs <- ws[points, 'reachID']
	conn <- ws$reach_connectivity
	if(self)
		diag(conn) <- rep(1, nrow(conn))
	mat <- conn[rIDs,, drop=FALSE]
	if(missing(names)) {
		rownames(mat) <- points
	} else {
		rownames(mat) <- names
	}
	mat
}

#' Construct a matrix identifying the nearest downstream neighbor from a list of sites
#' 
#' @param ws A watershed
#' @param x A vector of pixel IDs
#' @param names Optional vector of site names
#' @return A 2-column matrix; the first column gives the ID of a pixel, the second its nearest
#' downstream neighbor. Pixels in `x` that have no nearest neighbor are excluded.
#' If 'names' is specified, the output matrix will give site names instead of pixelIDs
#' @export
nearestDownstreamNeighbor <- function(ws, x, siteNames) {
	dmat <- downstreamDist(ws, x)
	res <- do.call(c, apply(dmat, 1, function(x) {
		if(all(x == 0)) {
			NULL
		} else {
			x[x==0] <- Inf
			names(x)[which.min(x)]
		}
	}))
	mat <- cbind(from=as.integer(names(res)), to=as.integer(res))
	if(!missing(siteNames)) {
		mat <- cbind(from=siteNames[match(mat[,1], x)], to=siteNames[match(mat[,2], x)])
	}
	mat
}


#' Build a downstream distance matrix for a list of sites
#' 
#' For each site, this function identifies the other sites that are downstream of it and computes
#' a 'distance' between each site and all of its downstream sites. By default, this is river,
#' distance, computed by summing the length of stream between the two sites. For any entry in
#' this distance matrix [i,j], the value is either a nonzero number giving the downstream
#' distance from site i to site j, or zero, indicating that j is not downstream of i.
#' 
#' For a more general (and much slower) distance computation, see [wsDistance()]
#' 
#' @param ws A watershed
#' @param x A vector of pixel ids
#' @param variable The variable to use for the distance metric
#' @param fun The function to apply between each pair of sites
#' @return A site by site distance [Matrix::sparseMatrix()]
#' @export
downstreamDist <- function(ws, x, variable = 'length', fun = sum) {
	res <- parallel::mclapply(x, function(i) {
		dspix <- connect(ws, i, Inf)
		dsSites <- x[which(x %in% dspix & !x==i)]
		if(length(dsSites) == 0) {
			NULL
		} else {
			t(sapply(dsSites, function(j) {
				ijPixes <- connect(ws, i, j)
				c(i, j, fun(ws[ijPixes, variable]))
			}))
		}
	})
	res <- do.call(rbind, res)
	res[,1] <- match(res[,1], x)
	res[,2] <- match(res[,2], x)
	Matrix::sparseMatrix(i = res[,1], j = res[,2], x = res[,3], 
		dims = c(length(x), length(x)), dimnames = list(x, x))
}



#' Find the nearest neighbors for a focal pixel among a list of sites
#' 
#' For a list of `sites`, this function first computes a distance matrix(if not specified),
#' then finds the nearest upstream and downstream neighbors for focal pixels `x`. Note that
#' `x` does not necessarily have to be in the list of sites, they can be anywhere in the watershed.
#' 
#' @param ws A watershed, optional unless `distMatrix` is missing
#' @param x A vector of focal pixels (always required)
#' @param distMatrix A site by pixel distance matrix, will be computed if not specified
#' @param sites A list of sites for the distance matrix, required if distMatrix is missing
#' @param selfAdjacency if TRUE, a site can be downstream of itself
#' 
#' @return A list with a downstream and an upstream element. Each element is a matrix, where
#' the first column corresponds to the focal pixels in `x`, and the second column lists the
#' nearest neighbors in `sites`.
#' @export
nearestNeighbors <- function(ws, x, distMatrix, sites, selfAdjacency = FALSE) {
	if(missing(distMatrix))
		distMatrix <- wsDistance(ws, sites)

	if(length(x) > 1) {
		res <- lapply(x, function(xx) nearestNeighbors(ws, xx, distMatrix, 
			selfAdjacency = selfAdjacency))
		resds <- do.call(rbind, lapply(res, function(xx) xx$downstream))
		resus <- do.call(rbind, lapply(res, function(xx) xx$upstream))
		return(list(downstream = resds, upstream = resus))
	}

	x <- as.character(x)

	# nearest downstream neighbor is easy, there is just one
	if(selfAdjacency & x %in% rownames(distMatrix)) {
		ds <- x
	} else if(! x %in% colnames(distMatrix)) {
		ds <- numeric(0)
	} else {
		ds <- dsNeighbors(x, distMatrix)
		ds <- rownames(distMatrix[ds,x, drop=FALSE])[which.max(distMatrix[ds,x])]
	}
	if(length(ds) > 0) {
		ds <- c(as.integer(x), as.integer(ds))
	} else
		ds <- NULL

	# for upstream, we look for upstream neighbors that have no downstream neighbors
	# within the set of all upstream neighbors (in other words, their only downstream neighbors)
	# are either the site x or sites downstream of x
	if(selfAdjacency & x %in% rownames(distMatrix)) {
		us <- x
	} else if(! x %in% colnames(distMatrix)) {
		us <- numeric(0)
	} else {
		us <- usNeighbors(x, distMatrix)
		us_ds <- lapply(us, dsNeighbors, distMatrix = distMatrix[us,,drop = FALSE])
		if(length(us_ds) == 0) {
			us <- character(0)
		} else
			us <- us[sapply(us_ds, function(xx) all(is.na(xx)))]
	}

	if(length(us) > 0) {
		us <- cbind(as.integer(x), as.integer(us))
	} else 
		us <- NULL


	list(downstream = ds, upstream = us)
}

# returns all downstream neighbors given a distance matrix
#' @keywords internal
dsNeighbors <- function(x, distMatrix) {
	nbs <- distMatrix[,x]
	names(nbs)[which(nbs < 0)]
}

# returns all upstream neighbors given a distance matrix
#' @keywords internal
usNeighbors <- function(x, distMatrix) {
	nbs <- distMatrix[,x]
	names(nbs)[which(nbs > 0)]
}
