
#' Remove short headwater reaches below a certain threshold
#' @param ws A watershed
#' @param size The minimum size of a reach to retain, in map units
#' @param rebuild Boolean; should topology be rebuilt? Normally yes, but can set to FALSE if it will be rebuilt later
#' @return A modified watershed
#' @export
trim_reaches = function(ws, size, rebuild = TRUE) {
	find_short_reaches = function(ws, size) {
		lens = tapply(ws$data$length, ws$data$reachID, sum)
		reaches = data.table::data.table(reachID = as.integer(names(lens)), length = unname(lens))
		reaches$hw = 0
		hw = headwaters(ws)$reachID
		reaches$hw[match(hw, reaches$reachID)] = 1
		reaches = reaches[length < size & hw == 1]
		pix_ids = which(ws$data$reachID %in% reaches$reachID)
		pix_ids
	}
	
	changed = FALSE
	pix_ids = find_short_reaches(ws, size)
	while(length(pix_ids) > 0) {
		changed = TRUE
		ws$data = ws$data[-pix_ids,]
		ws$data$id = 1:nrow(ws$data)
		ws$adjacency = ws$adjacency[-pix_ids, -pix_ids]
		rownames(ws$adjacency) = colnames(ws$adjacency) = ws$data$id
		pix_ids = find_short_reaches(ws, size)
	}
	if(rebuild && changed)
		ws = .rebuild_reaches(ws)
	return(ws)
}







#' Prepare watershed by setting up proper reach structures
#' @param x A watershed
#' @return A modified watershed with reaches properly set
#' @keywords internal
.rebuild_reaches = function(x) {
	# renumber all reaches
	old_reach_ids = unique(x$data$reachID)
	index_reach_ids = 1:max(old_reach_ids)
	if(!(all(index_reach_ids %in% old_reach_ids) && all(old_reach_ids %in% index_reach_ids))) {
		new_ids = rank(old_reach_ids)
		x$data$reachID = match(x$data$reachID, old_reach_ids)		
	}

	# make the topology & connectivity matrices
	ws_st = .make_ws_stack(x)
	x$reach_adjacency = Matrix::t(watershed::reach_topology(ws_st, Matrix::t(x$adjacency)))
	x$reach_connectivity = .create_reach_connectivity(x, self = FALSE)

	x
}


#' Produce a reach by reach connectivity matrix
#' 
#' @param ws A watershed
#' @param self If TRUE, a reach is considered connected to itself
#' Nonzero values indicate that a reachID in a row is downstream from that column
#' @return A [Matrix::sparseMatrix()]
#' @keywords internal
.create_reach_connectivity = function(ws, self = TRUE) {
	adj = ws$reach_adjacency
	if(self)
		diag(adj) = 1
	for(i in 1:nrow(adj))
		adj = adj %*% adj + adj
	adj[adj != 0] = 1
	adj
}





