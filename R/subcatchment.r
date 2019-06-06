#' Extract a complete subcatchment from the watershed
#' 
#' This function returns a complete subset watershed object, including all tributaties upstream
#' from a specified point. Pixel IDs will not be renamed, but reachIDs will.
#' @param ws A [Watershed()].
#' @param x The point to make the new outlet for the extracted watershed.
#' @return A watershed.
#' @export
subcatchment <- function(ws, x) {
	newSet <- connect(ws, downstream = x, upstream = Inf)
	dat <- ws[newSet,]
	oldReaches <- dat$reachID
	dat$reachID <- renumberReaches(dat$reachID)

	adj <- ws$adjacency[newSet, newSet]
	reach_adj <- ws$reach_adjacency[oldReaches, oldReaches]
	reach_conn <- ws$reach_connectivity[oldReaches, oldReaches]
	wsobj <- list(data = dat, adjacency = adj, reach_adjacency = reach_adj, 
			reach_connectivity = reach_conn)
	class(wsobj) <- class(ws)
	wsobj$reach_adjacency <- reachByReachAdj(wsobj)
	wsobj$reach_connectivity <- reachByReachConn(wsobj, self = FALSE)

	wsobj
}
