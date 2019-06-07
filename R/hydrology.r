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
	if(!('adjacency_q' %in% ls(ws)))
		ws$adjacency_q <- Matrix::t(Matrix::t(ws$adjacency) * ws$data$discharge)

	parms <- list(qout = Matrix::colSums(ws[["adjacency_q"]]), 
		qin = Matrix::rowSums(ws[["adjacency_q"]]), lateral = lateral, 
		csArea = ws[['data']][['csArea']], dx = ws[['data']][["dx"]])
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

#' Compute hydraulic geometry from discharge
#' @param discharge A vector of discharge values.
#' @param pars Optional parameters for the scaling relationships; see 'details'
#' @details If given, `pars` should be a named `list` with three elements: `velocity`, `depth`
#' 		and `width`. Each element must be a vector with two elements, the first is the intercept
#' 		and the second the slope. If omitted, values from Raymond et al will be used.
#' 
#' For the Raymond et al parameters, discharge should be in units of \eqn{m^3 s^{-1}}. 
#' 		All other variables are in meters.
#' @references Raymond, PA et al. 2012. Scaling the gas transfer velocity and hydraulic 
#' 		geometry in streams and small rivers. *Limnology and Oceanography: Fluids and
#' 		Environments*. **2**:41-53.
#' @return A data frame with discharge, velocity, depth, and width
#' @export
hydraulic_geometry <- function(discharge, pars) {
	if(missing(pars))
	{
		va <- -1.64
		vb <- 0.285
		da <- -0.895
		db <- 0.294
		wa <- 2.56
		wb <- 0.423
	} else {
		va <- pars$velocity[1]
		vb <- pars$velocity[2]
		da <- pars$depth[1]
		db <- pars$depth[2]
		wa <- pars$width[1]
		wb <- pars$width[2]
	}

	velocity <- exp(va + vb * log(discharge))
	depth <- exp(da + db * log(discharge))
	width <- exp(wa + wb * log(discharge))
	ind <- which(discharge == 0)
	velocity[ind] <- depth[ind] <- width[ind] <- 0
	data.frame(discharge, velocity, depth, width)
}


#' Compute discharge from catchment area
#' @param A vector of catchment area, in \eqn{m^2}
#' @param calib A data.frame with two elements, `A` (catchment area) and `Q` (observed discarge)
#' @details Computes discharge from catchment area as \eqn{\log Q = \log b + m \log A}.
#' 		If a single point in `calib` is included, the relationship will be re-parameterised by 
#'  	adjusting the intercept parameter `b` so that the calibration point falls on the line  
#' 		(while keeping the slope the same).
#' 
#' 		With multiple points, a bayesian regression model is fit with rstanarm. The model uses
#' 		Burgers et al parameters as an informative prior.
#' 
#' 		The default parameters used by this function come from Burgers et al (2014). 
#' 		For these parameters, catchment area units are expected to be in 
#' 		\eqn{m^2}, and discharge will be computed in \eqn{m^3 s^{-1}}.
#' @references Burgers HE et al. 2014. Size relationships of water discharge in rivers: scaling
#' 		of discharge with catchment area, main-stem lengthand precipitation. 
#' 		*Hydrological Processes*. **28**:5769-5775.
#' @export
discharge_scaling <- function(A, calib)
{

	# params from Burgers et al 2014; assuming normal distribution on log scale
	logB_mu <- log(6.3e-6)
	logB_sd <- 0.4137399
	m_mu <- 0.78
	m_sd <- 0.03571429

	# convert units to match pars from Burgers et al
	m2perkm2 <- 1000^2
	m3perkm3 <- 1000^3
	secperday <- 60*60*24
	A <- A / m2perkm2 ## now in square kilometers
	calib$A <- calib$A / m2perkm2
	calib$logA <- log(calib$A)
	calib$Q <- (calib$Q / m3perkm3) * secperday # now in cubic kilometers per day
	calib$logQ <- log(calib$Q)

	if(nrow(calib) == 1) {
		logB_mu <- log(calib$Q) - m_mu * log(calib$A)
		Q <- exp(logB_mu + m_mu * log(A))
	} else {
		if(!requireNamespace("rstanarm"))
			stop("This functionality requires the rstanarm package; 
				please install it and try again")
		# stan doesn't play nice with data tables
		calib <- as.data.frame(calib)
		fit <- rstanarm::stan_glm(logQ ~ logA, data = calib,
				prior_intercept = rstanarm::normal(logB_mu, logB_sd),
				prior = rstanarm::normal(m_mu, m_sd))
		Q <- exp(predict(fit, newdata = data.frame(logA = log(A))))
	}

	# convert to m^3 per second
	Q <- (Q * m3perkm3) / secperday
	return(Q)
}
