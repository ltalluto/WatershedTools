
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
		ws = .rebuild_reach_topology(ws)
	return(ws)
}



#' Reconstruct the reach topology
#' 
#' After editing a watershed's pixels or reaches, the reach topology will need to be rebuilt using this function.
#' @param ws A watershed
#' @return A modified watershed
#' @keywords internal
.rebuild_reach_topology = function(ws) {
	ws = .renumber_reaches(ws)
	ws$reach_adjacency = .create_reach_adjacency(ws)
	ws$reach_connectivity = .create_reach_connectivity(ws, self = FALSE)
	ws
}

#' Resize a single reach
#' @details reach numbers will start from 1
#' @param adj subsetted adjacency matrix from a watershed, with just the entries for a single reach
#' @param len length of each pixel
#' @param start The reach number to start with
#' @param size Reach size; the desired reach size, see 'details'
#' @param min_size Minimum reach size
.resize_reach = function(adj, len, size, min_size, start = 1) {
	if(sum(len) <= size) {
		rnums = rep(start, nrow(adj))
	} else {
		rnums = rep(0, nrow(adj))
	}
	c_reach = which(Matrix::rowSums(adj) == 0)
	while(any(rnums == 0)) {
		while(sum(len[c_reach]) < size) {
			nxt = which(adj[,c_reach[length(c_reach)]] == 1)
			c_reach = c(c_reach, nxt)
		}
		rnums[c_reach] = start
		ind = which(rnums == 0)
		if(length(ind) > 0) {
			if(sum(len[ind]) < min_size) {
				rnums[ind] = start
			} else if(sum(len[ind]) < size) {
				rnums[ind] = start + 1
			} else {
				start = start + 1
				c_reach = which(adj[,c_reach[length(c_reach)]] == 1)
			}
		}
	}

	return(rnums)
}





#' Produce a reach by reach connectivity matrix
#' 
#' @param ws A watershed
#' @param self If TRUE, a reach is considered connected to itself
#' Nonzero values indicate that a reachID in a row is downstream from that column
#' @return A [Matrix::sparseMatrix()]
#' @keywords internal
.create_reach_connectivity <- function(ws, self = TRUE) {
	adj <- ws$reach_adjacency
	if(self)
		diag(adj) <- 1
	for(i in 1:nrow(adj))
		adj <- adj %*% adj + adj
	adj[adj != 0] <- 1
	adj
}





