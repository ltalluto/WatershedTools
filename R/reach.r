#' Resize a watershed's reaches
#' @details By default watersheds are created with simple reaches that stretch between confluences. These can be 
#' quite long, however, so this function will resize them, starting at headwaters.
#' 
#' The size parameter only gives a target size, reaches can be smaller or larger in specific circumstances. Resizing
#' starts at headwaters, When a confluence is reached, the final reach before the confluence will often be smaller than
#' `size`. Thus, it is also necessary to specify `min_size`. If the reach upstream of a confluence is smaller than 
#' `min_size`, the reach will be incorporated into the next reach upstream. If there is no reach upstream, the reach will
#' be deleted. If it is larger than `min_size`, it will be left as a separate reach. Thus, the largest reaches will be
#' `size + min_size` in length, and the smallest will be `min_size`
#' 
#' @param ws A watershed
#' @param size Reach size; the desired reach size, see 'details'
#' @param min_size Minimum reach size
#' @param parallel Use the parallel package to speed computation on mac/linux?
#' @return A modified watershed
#' @export
resize_reaches = function(ws, size, min_size, parallel=TRUE) {
	ws = trim_reaches(ws, min_size, rebuild = FALSE)
	vals = rep(0, nrow(ws$data))
	reaches = unique(ws$data$reachID)
	pixes = lapply(reaches, function(x) which(ws$data$reachID == x))
	adj = lapply(pixes, function(x) ws$adjacency[x,x, drop=FALSE])
	len = lapply(pixes, function(x) ws$data$length[x])
	
	if(parallel) {
		nums = parallel::mcmapply(.resize_reach, adj, len, MoreArgs = list(size = size, min_size = min_size), SIMPLIFY=FALSE)
	} else {
		nums = mapply(.resize_reach, adj, len, MoreArgs = list(size = size, min_size = min_size), SIMPLIFY=FALSE)
	}

	if(length(nums) > 1) {
		for(i in 2:length(nums))
			nums[[i]] = nums[[i]] + max(nums[[i-1]])
	}
	for(i in seq_along(pixes)) {
		vals[pixes[[i]]] = nums[[i]]
	}
	ws$data$reachID = vals
	
	.rebuild_reach_topology(ws)
}

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

#' Splits a reachID at selected points
#' @param ws A watershed object
#' @param points A vector of ids at which to split reaches
#' @param na_ignore Logical, if `TRUE` NAs will be dropped from the data. If `FALSE`, any NAs
#' 		in `points` will cause an error.
#' @return A modified watershed with new reachIDs
#' @export
splitReaches <- function(ws, points, na_ignore = FALSE) {
	if(na_ignore & any(is.na(points))) {
		points <- points[!is.na(points)]
	} else if(any(is.na(points))) {
		stop("NA in points; use na_ignore = TRUE to ignore")
	}
	for(pt in points) {
		ids <- which(ws$data$reachID == ws$data$reachID[pt])
		reachAdj <- ws$adjacency[ids, ids]
		mostUpstream <- ids[which(Matrix::rowSums(reachAdj) == 0)]
		if(pt != mostUpstream) {
			chIds <- connect(ws, mostUpstream, pt)
			ws$data$reachID[chIds] <- max(ws$data$reachID) + 1
		}
	}
	
	# re-create topology
	.rebuild_reach_topology(ws)
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

#' Renumber reachIDs to start at 1 and increase one at a time, preserving ordering
#' @keywords internal
renumberReaches <- function(rIDs) {
	stop("Function is deprecated, use .renumber_reaches")
	newNums <- 1:length(unique(rIDs))
	oldNums <- sort(unique(rIDs))
	newNums[match(rIDs, oldNums)]
}

#' Renumber reaches, starting with one at the headwaters
#' @param ws A watershed
.renumber_reaches = function(ws) {
	result = rep(0, nrow(ws$data))
	reaches = headwaters(ws)$reachID
	while(any(result == 0)) {
		start = max(result)
		pixes = lapply(reaches, function(x) which(ws$data$reachID == x))
		vals = seq_along(pixes) + start
		newvals = mapply(function(pp, vv) rep(vv, length(pp)), pixes, vals, SIMPLIFY = FALSE)
		newvals = cbind(unlist(pixes), unlist(newvals))
		result[newvals[,1]] = newvals[,2]
		ids = outlets(ws, reaches)$id
		nxt = unique(unlist(apply(ws$adjacency[, ids, drop=FALSE], 2, function(x) which(x != 0))))
		nxt = nxt[which(result[nxt] == 0)]
		reaches = ws$data$reachID[nxt]
	}
	ws$data$reachID = result
	ws
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




#' Find all reaches directly upstream from a given reach
#' @keywords internal
reachAdj <- function(ws, rch) {
	ids <- which(ws$data$reachID == rch)
	reachAdj <- ws$adjacency[ids,ids, drop = FALSE]
	mostUpstream <- ids[which(Matrix::rowSums(reachAdj) == 0)]
	adjMatUp <- as.matrix(ws$adjacency[mostUpstream,])
	upRch <- which(adjMatUp == 1)
	if(length(upRch) > 0) {
		upRch <- ws$data$reachID[upRch]
		return(cbind(rch, upRch))
	} else return(NULL)
}

