#' Extract a complete subcatchment from the watershed
#' 
#' This function returns a complete subset watershed object, including all triburaties upstream
#' from a specified point. PixelIDs and reachIDs will be renumbered.
#' @param ws A [Watershed()].
#' @param x The point to make the new outlet for the extracted watershed.
#' @return A watershed.
#' @export
subcatchment = function(ws, x) {
	newSet = sort(connect(ws, downstream = x, upstream = Inf))
	newNums = 1:length(newSet)
	dat = ws$data[newSet,]
	dat$id = newNums

	adj = ws$adjacency[newSet, newSet]
	rownames(adj) = colnames(adj) = newNums
	wsobj = list(data = dat, adjacency = adj)
	class(wsobj) = class(ws)
	
	.rebuild_reaches(wsobj)
}
