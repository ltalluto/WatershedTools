# #' Produces a site by pixel distance matrix
# ws_distance <- function(ws, x) {

# 	## the approach here will be to use mclapply on connectAll for each point in x, 
# 	## get all pixels connected, then sum their lengths (*-1 for upstream)

# 	### argh, no, this won't work
# 	### the easiest thing to do is to accumulate WHILE extracting with connectCPP
# 	### but this means I need a CPP way to move upstream
# 	### perhaps the easiest way to using connectCPP as rewritten, find all upstream headwaters
# 	### then connect(and compute distance) from them down to the target pixel
# 	### so for pixel x;
# 		1. connectWithDist(x, outlet)
# 		2. connectUpstream(x)
# 		3. hws <- which(headwaters(ws) %in% connectUpstream(x))
# 		4. usDists <- -1 * mclapply(hws, connectWithDist(hw, x))
# 		5. remove redundant info in US dists
# 		6. make into matrix
# }







#' Accumulate a variable for one pixel to all up- and downstream pixels
#' 
#' @param ws A watershed
#' @param x A focal site
#' @param length The variable to accumulate
#' @return a vector of accumulated values (or NA for unconnected pixels). #' Accumulated values
#' will be positive going downstream, negative going upstream
accumulateOne <- function(ws, x, variable = 'length') {
	rid <- ws[x, 'reachID']
	inreach <- ws$data$id[ws[,'reachID'] == rid]

	upReaches <- which(ws$reach_connectivity[rid,] == 1)
	upPixes <- which(ws[,'reachID'] %in% upReaches)
	downReaches <- which(ws$reach_connectivity[,rid] == 1)
	downPixes <- which(ws[,'reachID'] %in% downReaches)

	upSum <- tapply(ws[upPixes, variable], ws[upPixes, 'reachID'], function(x) -1 * sum(x))
	downSum <- tapply(ws[downPixes, variable], ws[downPixes, 'reachID'], function(x) sum(x))

	adj <- ws$reach_adjacency[c(rid, upReaches), c(rid, upReaches)]
	tops <- Matrix::which(Matrix::rowSums(adj) == 0)	
	usReaches <- downstreamPixelIds(list(adjacency = adj))
	usLengths <- parallel::mclapply(tops, function(x) {
		res <- do.call(cbind, connectCPP(usReaches, x, 1, c(0, upSum)))
		res[,2] <- res[nrow(res), 2] - res[,2]
		res
	})
	usLengths <- do.call(rbind, usLengths)
	usLengths <- usLengths[!duplicated(usLengths[,1]), ]
	usLengths[,1] <- c(rid, upReaches)[usLengths[,1]] ## restore reachIDs

	adj <- ws$reach_adjacency[c(rid, downReaches), c(rid, downReaches)]
	bottom <- Matrix::which(Matrix::colSums(adj) == 0)
	dsReaches <- downstreamPixelIds(list(adjacency = adj))
	dsLengths <- do.call(cbind, connectCPP(dsReaches, 1, bottom, c(downSum, 0)))
	dsLengths <- dsLengths[-1, ]
	dsLengths[,1] <- c(rid, downReaches)[dsLengths[,1]] ## restore reachIDs


	inReachLens <- accumInReach(ws, x, variable)
	usLengths[,2] <- usLengths[,2] + min(inReachLens[,2]) # add the upstream in-reach length
	dsLengths[,2] <- dsLengths[,2] + max(inReachLens[,2]) # add the downstream in-reach length

	reachLengths <- rbind(usLengths, dsLengths)
	pixLens <- reachLengths[match(ws$data$reachID, reachLengths[,1]),2]
	pixLens[inReachLens[,1]] <- inReachLens[,2]
	pixLens
}

#' Accumulates a variable within a reach
#' @param ws A watershed
#' @param x A focal site
#' @param length The variable to accumulate
#' @return A vector of accumulated values (or NA for unconnected pixels). 
#' Accumulated values
#' will be positive going downstream, negative going upstream
accumInReach <- function(ws, x, variable = 'length') {
	rid <- ws[x, 'reachID']
	inreach <- ws$data$id[ws[,'reachID'] == rid]

	adj <- ws$adjacency[inreach,inreach]
	# rownames(adj) <- colnames(adj) <- 1:nrow(adj)
	dsPixes <- downstreamPixelIds(list(adjacency = adj))
	bottom <- Matrix::which(Matrix::colSums(adj) == 0)
	top <- Matrix::which(Matrix::rowSums(adj) == 0)
	pix <- which(inreach == x)

	downAccum <- do.call(cbind, connectCPP(dsPixes, pix, bottom, ws[inreach,variable]))
	upAccum <- do.call(cbind, connectCPP(dsPixes, top, pix, ws[inreach,variable]))
	upAccum[,2] <- upAccum[,2] - upAccum[nrow(upAccum), 2]
	upAccum <- upAccum[-which(upAccum[,1] == pix),]

	# restore pixelIDs
	res <- rbind(upAccum, downAccum)
	res[,1] <- inreach[res[,1]]
	res
}
