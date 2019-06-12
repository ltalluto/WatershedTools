#' Extract a complete subcatchment from the watershed
#' 
#' This function returns a complete subset watershed object, including all tributaties upstream
#' from a specified point. PixelIDs and reachIDs will be renumbered.
#' @param ws A [Watershed()].
#' @param x The point to make the new outlet for the extracted watershed.
#' @return A watershed.
#' @export
subcatchment <- function(ws, x) {
	newSet <- sort(connect(ws, downstream = x, upstream = Inf))
	newNums <- 1:length(newSet)
	dat <- ws[newSet,]
	dat$id <- newNums
	oldReaches <- dat$reachID
	dat$reachID <- renumberReaches(dat$reachID)

	adj <- ws$adjacency[newSet, newSet]
	rownames(adj) <- colnames(adj) <- newNums
	reach_adj <- ws$reach_adjacency[oldReaches, oldReaches]
	rownames(reach_adj) <- colnames(reach_adj) <- dat$reachID
	reach_conn <- ws$reach_connectivity[oldReaches, oldReaches]
	rownames(reach_conn) <- colnames(reach_conn) <- dat$reachID
	wsobj <- list(data = dat, adjacency = adj, reach_adjacency = reach_adj, 
			reach_connectivity = reach_conn)
	class(wsobj) <- class(ws)
	wsobj$reach_adjacency <- reachByReachAdj(wsobj)
	wsobj$reach_connectivity <- reachByReachConn(wsobj, self = FALSE)

	wsobj
}
