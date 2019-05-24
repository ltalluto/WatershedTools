#' Produces a site by pixel distance matrix
ws_distance <- function(ws, x) {

	## the approach here will be to use mclapply on connectAll for each point in x, 
	## get all pixels connected, then sum their lengths (*-1 for upstream)

	### argh, no, this won't work
	### the easiest thing to do is to accumulate WHILE extracting with connectCPP
	### but this means I need a CPP way to move upstream
	### perhaps the easiest way to using connectCPP as rewritten, find all upstream headwaters
	### then connect(and compute distance) from them down to the target pixel
	### so for pixel x;
		1. connectWithDist(x, outlet)
		2. connectUpstream(x)
		3. hws <- which(headwaters(ws) %in% connectUpstream(x))
		4. usDists <- -1 * mclapply(hws, connectWithDist(hw, x))
		5. remove redundant info in US dists
		6. make into matrix
}



#' For a single pixel, find all upstream and downstream pixels
#' Slow version for historical purposes
#' @param ws A watershed
#' @param x The pixel to connect
# connectAll <- function(ws, x) {
# 	rid <- ws[x, 'reachID']
# 	inreach <- ws$data$id[ws[,'reachID'] == rid]

# 	upReaches <- which(ws$connectivity[rid,] == 1)
# 	downReaches <- which(ws$connectivity[,rid] == 1)
# 	direction <- rep(0, length(inreach)) ## 1 is downstream, -1 is up
# 	out <- outlets(ws)$id
# 	pid <- x
# 	curr_rid <- rid
# 	## need to speed this up
# 	while(curr_rid == rid & pid != out) {
# 		pid <- which(ws$adjacency[,pid] == 1)
# 		curr_rid <- ws[pid, 'reachID']
# 		direction[which(inreach == pid)] <- 1
# 	}
# 	direction[direction == 0 & inreach != x] <- -1

# 	upPixes <- c(ws$data$id[ws$data$reachID %in% upReaches], inreach[direction == -1])
# 	downPixes <- c(ws$data$id[ws$data$reachID %in% downReaches], inreach[direction == 1])

# 	list(upstream = upPixes, downstream = downPixes)
# }


#' Yo this works and it's a fast way to get every pixel that's connected to a focal pixel
#' Probably need to use it elsewhere
connectAll <- function(ws, x) {
	rid <- ws[x, 'reachID']
	inreach <- ws$data$id[ws[,'reachID'] == rid]

	upReaches <- which(ws$connectivity[rid,] == 1)
	downReaches <- which(ws$connectivity[,rid] == 1)
	direction <- rep(0, length(inreach)) ## 1 is downstream, -1 is up

	adj <- ws$adjacency[inreach,inreach]
	rownames(adj) <- colnames(adj) <- 1:nrow(adj)
	dsPixes <- downstreamPixelIds(list(adjacency = adj))
	ds <- Matrix::which(Matrix::colSums(adj) == 0)
	us <- which(inreach == x)
	downPixels <- connectCPP(dsPixes, us, ds)
	direction[downPixels] <- 1
	direction[direction == 0] <- -1
	direction[which(inreach == x)] <- 0

	upPixes <- c(ws$data$id[ws$data$reachID %in% upReaches], inreach[direction == -1])
	downPixes <- c(ws$data$id[ws$data$reachID %in% downReaches], inreach[direction == 1])

	list(upstream = upPixes, downstream = downPixes)
}
