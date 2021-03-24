#' Runs a transport-reaction model
#' @param ws A [Watershed] object
#' @param initial A vector of initial conditions
#' @param lateral Vector of values of the state variable for lateral input for each pixel
#' @param times The times at which to compute the state variable
#' @param method The integration method to use; default is euler
#' @param dt The size of the time step for integration (euler method only)
#' @param rxn A function giving the time derivative of the reaction component
#' @param rxnParams A list of parameters to pass to `rxn`
#' 
#' @details Because the appropriate time step will vary based on the units and the specific
#'    problem, the recommended approach is to first try a single run using lsoda (perhaps on
#'    a small subset), then run euler at varying time steps to find the dt that produces
#'    acceptable results.
#' 
#' If used, `rxn` must be a function taking `t` (the current time) and `y` 
#'   (the state of the system) as its first few arguments. Other arguments can be added by 
#'   name via the `rxnParams` argument.
#' @return A state variable matrix (with one column per time step)
#' @useDynLib WatershedTools
#' @importFrom Rcpp sourceCpp
#' @export
transport <- function(ws, initial, lateral, times, method = c('euler', 'lsoda'), dt = 1, 
	rxn = NULL, rxnParams = list()) {
	method <- match.arg(method)
	if(method == "lsoda" && !requireNamespace("deSolve"))
		stop("method = 'lsoda' requires the deSolve package; please install it and try again")
	if(times[1] != 0)
		stop('first time step must be 0')
	if(method == 'euler' && sum(times %% dt) != 0)
		stop('dt must divide evenly into all values in times')
	if(!'discharge' %in% names(ws))
		stop('discharge is required for the transport model')
	if(!'csArea' %in% names(ws)) {
		if(!('width' %in% names(ws)) & !('depth' %in% names(ws))) {
			stop("depth and width, or csArea, are required attributes")
		} else
			ws$data$csArea <- ws$data$width * ws$data$depth
	}
	if(!('adjacency_q' %in% ls(ws)))
		ws$adjacency_q <- Matrix::t(Matrix::t(ws$adjacency) * ws$data$discharge)

	parms <- list(qout = Matrix::colSums(ws[["adjacency_q"]]), 
		qin = Matrix::rowSums(ws[["adjacency_q"]]), lateral = lateral, 
		csArea = ws[['data']][['csArea']], dx = ws[['data']][["length"]])
	# add discharge for the outlet back in (outlet drains to nowhere)
	parms$qout[which(parms$qout == 0)] <- ws[['data']]$discharge[which(parms$qout == 0)]
	
	if(method == "euler") {
		parms$adjacencyQ <- ws[["adjacency_q"]]
		state <- transport_euler(initial, parms, times, dt, rxn, rxnParams)
	} else {
		parms$adjacencyQ <- cbind(Matrix::which(ws[["adjacency_q"]] != 0, arr.ind = TRUE), 
			ws[["adjacency_q"]][Matrix::which(ws[["adjacency_q"]] != 0)])
		state <- t(deSolve::ode(initial, times = times, func = dCdt_transport, parms = parms, 
			rxn = rxn, rxnParams = rxnParams))
		state <- state[2:nrow(state),] # dropping the times row that desolve adds on
	}
	return(state)
}


#' Do euler integration for the transport model
#' @param initial A vector of initial conditions
#' @param parms Integration parameters and data
#' @param times The times at which to compute the state variable
#' @param dt The size of the time step for integration
#' @param rxn A function giving the time derivative of the reaction component
#' @param rxnParams A list of parameters to pass to `rxn`
#' @return A state variable matrix (with one column per time step)
#' @keywords internal
transport_euler <- function(initial, parms, times, dt, rxn = NULL, rxnParams = list()) {
	state <- matrix(NA, ncol=length(times), nrow=length(initial))
	state[,1] <- statet <- initial
	t <- 0
	for(i in 2:ncol(state)) {
		while(t < times[i]) {
			statet <- statet + dt * dCdt_transport_r(t, statet, parms)
			if(!is.null(rxn)) {
				rxnList <- c(list(t = t, y = statet), (rxnParams))
				statet <- statet + dt * do.call(rxn, rxnList)
			}
			t <- t + dt
		}
		state[,i] = statet[,1]
	}
	return(state)
}


#' @name dCdt_transport
#' @aliases dCdt_transport_r
#' @title Compute flux due to transport
#' @param t The time step
#' @param y The state variable
#' @param parms Model parameters and data
#' @param rxn A function giving the time derivative of the reaction component
#' @param rxnParams A list of parameters to pass to `rxn`
#' @details At present, only advective transport is considered
#' @return Flux vector, in units of state variable per unit time, for each spatial unit
#' @keywords internal
dCdt_transport <- function(t, y, parms, rxn = NULL, rxnParams = list()) {
	advection <- dCdt_transport_cpp(t, y, parms$adjacencyQ, parms$qout, parms$qin, 
		parms$lateral, parms$csArea, parms$dx)
	if(!is.null(rxn)) {
		rxnList <- c(list(t = t, y = y), (rxnParams))
		reaction <- do.call(rxn, rxnParams)
	} else {
		reaction <- 0
	}

	return(list(advection + reaction))
}

#' @rdname dCdt_transport
#' @keywords internal
dCdt_transport_r <- function(t, y, parms) {
	inputMass <- parms$adjacencyQ %*% y
	totalInputMass <- inputMass + (parms$qout - parms$qin) * parms$lateral;
	totalOutputMass <- parms$qout * y;
	advection <- (-1/parms$csArea) * (totalOutputMass - totalInputMass)/parms$dx;
	return(advection)
}







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
q_from_flow <- function(depth, velocity, flowID = rep(1, length(depth)), width, 
		widthID = rep(1, length(width)), na.rm = FALSE) {

	if(length(depth) != length(velocity) | length(depth) != length(flowID))
		stop("depth, velocity, and flowID must have the same length")
	if(length(width) != length(widthID))
		stop("width and widthID must have the same length")

	# require(data.table)
	widthByID <- data.table(id = widthID, width = width)
	widthByID <- widthByID[, .(width = mean(width, na.rm = na.rm)), keyby=id]

	flow <- data.table(id = flowID, velocity = velocity, depth = depth, key = 'id')
	if(any(!complete.cases(flow)) & na.rm) {
		flow <- flow[complete.cases(flow)]
	}

	## need to divide width by the number of measurements
	n <- tapply(flow$id, flow$id, length)
	n <- data.table(id =as.integer(names(n)), n = n, key = 'id')

	widthByID <- n[widthByID]

	flow <- widthByID[flow]
	flow[, .(Q = sum((width/n)*velocity*depth)), id]
}
