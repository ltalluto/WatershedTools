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
#' @return A modified watershed
resize_reaches = function(ws, size, min_size) {
	ws = .trim_reaches(ws, min_size)
	vals = rep(0, nrow(ws$data))
	reaches = unique(ws$data$reachID)
	pixes = lapply(reaches, function(x) which(ws$data$reachID == x))
	adj = lapply(pixes, function(x) ws$adjacency[x,x, drop=FALSE])
	len = lapply(pixes, function(x) ws$data$length[x])
	nums = mapply(.resize_reach, adj, len, MoreArgs = list(size = size, min_size = min_size), SIMPLIFY=FALSE)
	if(length(nums) > 1) {
		for(i in 2:length(nums))
			nums[[i]] = nums[[i]] + max(nums[[i-1]])
	}
	for(i in seq_along(pixes)) {
		vals[pixes[[i]]] = nums[[i]]
	}
	ws$data$reachID = vals
	
	## now renumber them, starting with the headwaters
	ws = .new_renumber_reaches_headwaters(ws)
	ws$reach_adjacency = reachByReachAdj(ws)
	ws$reach_connectivity = reachByReachConn(ws, self = FALSE)
	return(ws)
}

#' Remove short headwater reaches below a certain threshold
.trim_reaches = function(ws, size) {
	stop("finish me")
	hw = headwaters(ws)$reachID

}

#' Resize a single reach
#' @details reach numbers will start from 1
#' @param adj subsetted adjacency matrix from a watershed, with just the entries for a single reach
#' @param len length of each pixel
#' @param start The reach number to start with
#' @param size Reach size; the desired reach size, see 'details'
#' @param min_size Minimum reach size
.resize_reach = function(adj, len, size, min_size, start = 1) {
	rnums = rep(0, nrow(adj))
	c_reach = which(rowSums(adj) == 0)
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

#' Renumber reaches, starting with one at the headwaters
#' @param ws A watershed
.new_renumber_reaches_headwaters = function(ws) {
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
