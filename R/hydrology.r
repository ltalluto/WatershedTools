#' Runs a transport-only model
#' @param ws A [Watershed] object
#' @param initial A vector of initial conditions
#' @param lateral Vector of values of the state variable for lateral input for each pixel
#' @param times The times at which to compute the state variable
#' @param method The integration method to use; default is euler
#' @param dt The size of the time step for integration (euler method only)
#' 
#' @details Because the appropriate time step will vary based on the units and the specific
#'    problem, the recommended approach is to first try a single run using lsoda (perhaps on
#'    a small subset), then run euler at varying time steps to find the dt that produces
#'    acceptable results.
#' @return A state variable matrix (with one column per time step)
#' @useDynLib WatershedTools
#' @importFrom Rcpp sourceCpp
#' @export
transport <- function(ws, initial, lateral, times, method = c('euler', 'lsoda'), dt = 1) {
	method <- match.arg(method)
	if(method == "lsoda" && !requireNamespace("deSolve"))
		stop("method = 'lsoda' requires the deSolve package; please install it and try again")
	if(times[1] != 0)
		stop('first time step must be 0')
	if(method == 'euler' && sum(times %% dt) != 0)
		stop('dt must divide evenly into all values in times')
	if(!('adjacency_q' %in% names(ws)))
		ws[["adjacency_q"]] <- Matrix::t(Matrix::t(ws[['adjacency']]) * ws$data$discharge)

	parms <- list(qout = Matrix::colSums(ws[["adjacency_q"]]), 
		qin = Matrix::rowSums(ws[["adjacency_q"]]), lateral = lateral, 
		csArea = ws[['data']][['csArea']], dx = ws[['data']][["dx"]])
	# add discharge for the outlet back in (outlet drains to nowhere)
	parms$qout[which(parms$qout == 0)] <- ws[['data']]$discharge[which(parms$qout == 0)]
	
	if(method == "euler") {
		parms$adjacencyQ <- ws[["adjacency_q"]]
		state <- transport_euler(initial, parms, times, dt)
	} else {
		parms$adjacencyQ <- cbind(Matrix::which(ws[["adjacency_q"]] != 0, arr.ind = TRUE), 
			ws[["adjacency_q"]][Matrix::which(ws[["adjacency_q"]] != 0)])
		state <- t(deSolve::ode(initial, times = times, func = dCdt_transport, parms = parms))
		state <- state[2:nrow(state),] # dropping the times row that desolve adds on
	}
	return(state)
}

#' Do euler integration for the transport model
#' @param initial A vector of initial conditions
#' @param parms Integration parameters and data
#' @param times The times at which to compute the state variable
#' @param dt The size of the time step for integration
#' @return A state variable matrix (with one column per time step)
#' @keywords internal
transport_euler <- function(initial, parms, times, dt) {
	state <- matrix(NA, ncol=length(times), nrow=length(initial))
	state[,1] <- statet <- initial
	t <- 0
	for(i in 2:ncol(state)) {
		while(t < times[i]) {
			statet <- statet + dt * dCdt_transport_r(t, statet, parms)
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
#' @details At present, only advective transport is considered
#' @return Flux vector, in units of state variable per unit time, for each spatial unit
#' @keywords internal
dCdt_transport <- function(t, y, parms) {
	advection <- dCdt_transport_cpp(t, y, parms$adjacencyQ, parms$qout, parms$qin, 
		parms$lateral, parms$csArea, parms$dx)
	return(list(advection))
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

