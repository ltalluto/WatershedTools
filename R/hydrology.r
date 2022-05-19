#' Compute discharge from flow measurements
#' 
#' 
#' @param depth a vector of depth measurements
#' @param velocity a vector of velocity measurements
#' @param flowID a vector of IDs (e.g., site numbers) corresponding to depth/velocity measurements
#' @param width a vector of stream widths
#' @param widthID a vector of IDs for width measurements to match to flowID
#' @param na.rm should NAs be removed before computation?
#' @import data.table
#' @keywords internal
#' 
#' @return A data.table giving discharge by ID
q_from_flow = function(depth, velocity, flowID = rep(1, length(depth)), width, 
		widthID = rep(1, length(width)), na.rm = FALSE) {

	if(length(depth) != length(velocity) | length(depth) != length(flowID))
		stop("depth, velocity, and flowID must have the same length")
	if(length(width) != length(widthID))
		stop("width and widthID must have the same length")

	# require(data.table)
	widthByID = data.table(id = widthID, width = width)
	widthByID = widthByID[, .(width = mean(width, na.rm = na.rm)), keyby=id]

	flow = data.table(id = flowID, velocity = velocity, depth = depth, key = 'id')
	if(any(!complete.cases(flow)) & na.rm) {
		flow = flow[complete.cases(flow)]
	}

	## need to divide width by the number of measurements
	n = tapply(flow$id, flow$id, length)
	n = data.table(id =names(n), n = as.vector(n), key = 'id')
	if(is.numeric(flowID))
		n$id = as.integer(n$id)
	setkey(n, id)

	widthByID = n[widthByID]

	flow = widthByID[flow]
	flow[, .(Q = sum((width/n)*velocity*depth)), id]
}
