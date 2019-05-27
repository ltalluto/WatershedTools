# #' Produces a site by pixel distance matrix
wsDistance <- function(ws, x, variable = 'length') {
	downs <- parallel::mclapply(x, function(xx) accumulate(ws, xx, parallel=FALSE))
	ups <- parallel::mclapply(x, function(xx) 
		accumulate(ws, upstream=Inf, downstream = xx, parallel = FALSE))

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
#' @param ws Watershed object
#' @param upstream ID number of upstream point; if Inf ALL upstream pixels will be returned
#' @param downstream ID of downstream point; if Inf ALL downstream pixels will be returned
#' @param variable The variable to accumulate
#' @param parallel Boolean, should parallel processing be used?
#' @return A matrix, first column connected pixelIDs, second the accumulated variable
#' @export
accumulate <- function(ws, upstream, downstream = Inf, variable = 'length', parallel = FALSE) {
	if(parallel && !requireNamespace("parallel"))
		parallel <- FALSE
	dsPixes <- downstreamPixelIds(ws)
	if(is.infinite(downstream))
		downstream <- outlets(ws)$id
	if(is.infinite(upstream)) {
		rid <- ws[downstream, 'reachID']
		tops <- headwaters(ws)[,'id']
		upReaches <- which(ws$reach_connectivity[rid,] == 1)
		upstream <- tops[ws[tops, "reachID"] %in% upReaches]
	}
	if(length(downstream == 1)) {
		if(length(upstream) == 1) {
			accum <- do.call(cbind, connectCPP(dsPixes, upstream, downstream, ws[, variable]))
		} else {
			if(parallel) {
				accum <- do.call(rbind, 
					parallel::mclapply(upstream, function(x) 
						accumulate(ws, x, downstream, variable)))
			} else {
				accum <- do.call(rbind, lapply(upstream, function(x) 
						accumulate(ws, x, downstream, variable)))
			}
		}
	} else {
		accum <- do.call(rbind, 
				lapply(downstream, function(x) accumulate(ws, upstream, x, variable, parallel)))
	}
	return(accum)
}


# #' Accumulate a variable for between pixels
# #' 
# #' @param ws A watershed
# #' @param x A focal site
# #' @param y A second site
# #' @param length The variable to accumulate
# #' @return a vector of accumulated values (or NA for unconnected pixels). #' Accumulated values
# #' will be positive going downstream, negative going upstream
# accumulateBetween <- function(ws, x, y, variable = 'length') {
# 	rid_x <- ws[x, 'reachID']
# 	rid_y <- ws[y, 'reachID']

# 	if(rid_x == rid_y) {
# 		## handle special case when they are in the same reach
# 	} else {
# 		upReaches_x <- which(ws$reach_connectivity[rid_x,] == 1)
# 		upReaches_y <- which(ws$reach_connectivity[rid_y,] == 1)
# 		downReaches_x <- which(ws$reach_connectivity[,rid_x] == 1)
# 		downReaches_y <- which(ws$reach_connectivity[,rid_y] == 1)

# 		## check if they are not connected
# 		if(!(rid_y %in% downReaches_x) & !(rid_x %in% downReaches_y))
# 			return(NA)

# 		inReachLens_x <- accumInReach(ws, x, variable)
# 		inReachLens_y <- accumInReach(ws, y, variable)

# 		## y is downstream of x
# 		if(rid_y %in% downReaches_x) {
# 			sharedReaches <- downReaches_x[downReaches_x %in% upReaches_y]
# 			inReachSum <- sum(inReachLens_x[inReachLens_x > 0]) + 
# 					abs(sum(inReachLens_y[inReachLens_y < 0]))
# 			mult <- 1					
# 		} else if(rid_x %in% downReaches_y) {
# 			## y is upstream of x
# 			sharedReaches <- downReaches_y[downReaches_y %in% upReaches_x]
# 			inReachSum <- abs(sum(inReachLens_x[inReachLens_x < 0])) + 
# 					sum(inReachLens_y[inReachLens_y > 0])
# 			mult <- -1					
# 		}

# 		pixes <- which(ws$data$reachID %in% sharedReaches)
# 		res <- mult * (sum(ws[pixes, variable]) + inReachSum)
# 	}
# 	res
# }



# #' Accumulate a variable for one pixel to all up- and downstream pixels
# #' 
# #' @param ws A watershed
# #' @param x A focal site
# #' @param length The variable to accumulate
# #' @return a vector of accumulated values (or NA for unconnected pixels). #' Accumulated values
# #' will be positive going downstream, negative going upstream
# accumulateOne <- function(ws, x, variable = 'length') {
# 	rid <- ws[x, 'reachID']
# 	upReaches <- which(ws$reach_connectivity[rid,] == 1)
# 	upPixes <- which(ws[,'reachID'] %in% upReaches)
# 	downReaches <- which(ws$reach_connectivity[,rid] == 1)
# 	downPixes <- which(ws[,'reachID'] %in% downReaches)

# 	upSum <- tapply(ws[upPixes, variable], ws[upPixes, 'reachID'], function(x) -1 * sum(x))
# 	downSum <- tapply(ws[downPixes, variable], ws[downPixes, 'reachID'], function(x) sum(x))

# 	adj <- ws$reach_adjacency[c(rid, upReaches), c(rid, upReaches)]
# 	tops <- Matrix::which(Matrix::rowSums(adj) == 0)	
# 	usReaches <- downstreamPixelIds(list(adjacency = adj))
# 	usLengths <- parallel::mclapply(tops, function(x) {
# 		res <- do.call(cbind, connectCPP(usReaches, x, 1, c(0, upSum)))
# 		res[,2] <- res[nrow(res), 2] - res[,2]
# 		res
# 	})
# 	usLengths <- do.call(rbind, usLengths)
# 	usLengths <- usLengths[!duplicated(usLengths[,1]), ]
# 	usLengths[,1] <- c(rid, upReaches)[usLengths[,1]] ## restore reachIDs

# 	adj <- ws$reach_adjacency[c(rid, downReaches), c(rid, downReaches)]
# 	bottom <- Matrix::which(Matrix::colSums(adj) == 0)
# 	dsReaches <- downstreamPixelIds(list(adjacency = adj))
# 	dsLengths <- do.call(cbind, connectCPP(dsReaches, 1, bottom, c(downSum, 0)))
# 	dsLengths <- dsLengths[-1, ]
# 	dsLengths[,1] <- c(rid, downReaches)[dsLengths[,1]] ## restore reachIDs


# 	inReachLens <- accumInReach(ws, x, variable)
# 	usLengths[,2] <- usLengths[,2] + min(inReachLens[,2]) # add the upstream in-reach length
# 	dsLengths[,2] <- dsLengths[,2] + max(inReachLens[,2]) # add the downstream in-reach length

# 	reachLengths <- rbind(usLengths, dsLengths)
# 	pixLens <- reachLengths[match(ws$data$reachID, reachLengths[,1]),2]
# 	pixLens[inReachLens[,1]] <- inReachLens[,2]
# 	pixLens
# }

# #' Accumulates a variable within a reach
# #' @param ws A watershed
# #' @param x A focal site
# #' @param length The variable to accumulate
# #' @return A vector of accumulated values (or NA for unconnected pixels). 
# #' Accumulated values
# #' will be positive going downstream, negative going upstream
# accumInReach <- function(ws, x, variable = 'length') {
# 	rid <- ws[x, 'reachID']
# 	inreach <- ws$data$id[ws[,'reachID'] == rid]

# 	adj <- ws$adjacency[inreach,inreach]
# 	# rownames(adj) <- colnames(adj) <- 1:nrow(adj)
# 	dsPixes <- downstreamPixelIds(list(adjacency = adj))
# 	bottom <- Matrix::which(Matrix::colSums(adj) == 0)
# 	top <- Matrix::which(Matrix::rowSums(adj) == 0)
# 	pix <- which(inreach == x)

# 	downAccum <- do.call(cbind, connectCPP(dsPixes, pix, bottom, ws[inreach,variable]))
# 	upAccum <- do.call(cbind, connectCPP(dsPixes, top, pix, ws[inreach,variable]))
# 	upAccum[,2] <- upAccum[,2] - upAccum[nrow(upAccum), 2]
# 	upAccum <- upAccum[-which(upAccum[,1] == pix),]

# 	# restore pixelIDs
# 	res <- rbind(upAccum, downAccum)
# 	res[,1] <- inreach[res[,1]]
# 	res
# }
