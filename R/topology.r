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
	downs <- parallel::mclapply(x, function(xx) accumulate(ws, xx, parallel=FALSE))
	ups <- parallel::mclapply(x, function(xx) 
		accumulate(ws, upstream=Inf, downstream = xx, parallel = FALSE, direction = 'up'))
	matTall <- do.call(rbind.data.frame, mapply(function(xx, ds, us) {
		dsMat <- cbind(row=rep(xx, nrow(ds)), col = ds[,1], val = ds[,2])
		usMat <- cbind(row=rep(xx, nrow(us)), col = us[,1], val = us[,2])
		rbind(dsMat, usMat)
	}, xx = x, ds = downs, us = ups))
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
	accumulate(ws, upstream, downstream)[,1]
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
	if(is.infinite(downstream))
		downstream <- outlets(ws)$id
	if(is.infinite(upstream)) {
		rid <- ws[downstream, 'reachID']
		tops <- headwaters(ws)[,'id']
		upReaches <- which(ws$reach_connectivity[rid,] == 1)
		if(rid %in% ws[tops, 'reachID'])
			upReaches <- c(upReaches, rid)
		upstream <- tops[ws[tops, "reachID"] %in% upReaches]
	}
	if(length(downstream == 1)) {
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
				lapply(downstream, function(x) accumulate(ws, upstream, x, variable, parallel)))
	}
	## remove redundancies; it is quite common to traverse the same pixel multiple times
	accum <- accum[!duplicated(accum[,1]),]
	return(accum)
}

