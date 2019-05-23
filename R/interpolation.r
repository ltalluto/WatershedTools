#' Provide interpolation in time and space
#' @param ws A watershed
#' 
#' @details The data are expected in 'tall' format, so `sites`, `times`, and `values` must all
#' have the same length if they are specified. It is required to specify `values` plust at least
#' one of `sites` and `times. See below to determine how prediction proceeds based on
#' which optional parameters are specified.
#' 
#' If only `sites` is specified:
#' 		If `sites_out` is missing, prediction will be to all sites in the watershed. Otherwise
#' 		prediction is to the sites in `sites_out`. `times_out` will be ignored.
#'
#' If only `times` is specified:
#' 		Prediction will be made either to the times specified in `times_out` 
#' 		(which is required). Prediction will always be only for the input site,
#' 		`sites_out` is ignored.
#' 
#' If `sites` and `times` are both provided:
#' 		Prediction will be to `sites_out` if not missing, otherwise all sites, and to `times_out`
#' 		if not missing, otherwise `times`. If interpolating in both space and time, interpolation
#' 		is first done in time to the input sites, then in space.
#' 
#' When interpolation in space is desired, we weight by discharge and the inverse square of the 
#' distance.
#' nearest 
#' 
#' @param ws A watershed
#' @param values A vector of measured values to interpolate
#' @param sites An option vector of pixel IDs where measurements have been performed.
#' @param times An optional vector of times (e.g., in minutes) at which measurements were made.
#' @param sites_out Optional vector of target sites. See "details." 
#' @param times_out Optional vector of target times. See "details." 
#' @param distMatrix Optional distance matrix of the sites to all points in the watershed, as
#' created by [downstreamDist()]
interpolate <- function(ws, values, sites, times, sites_out, times_out, distMatrix) {
	if(missing(sites) & missing(times))
		stop("either sites or times must be specified")

	if(!missing(sites) && length(values) != length(sites))
		stop("sites and values must have the same length")

	if(!missing(times) && length(values) != length(times))
		stop("times and values must have the same length")

	if(!missing(sites) & missing(distMatrix))
		distMatrix <- downstreamDist(ws, sites)

	if(!missing(times) & missing(sites)) {
		if(missing(times_out))
			stop("times_out must be provided if only interpolating in time")
		result <- matrix(approx(x = times, y = values, xout = times_out), nrow=1)

	} else if(!missing(sites) & missing(times)) {
		if(missing(sites_out))
			sites_out <- ws[,'id']

		STEPS for each site_out
		1. compute nearest downstream neighbor site
		2. compute all independent upstream neighbor sites -- this requires a complete site by pixel distMatrix
		3. weithted avg, where weights are (Q_upstream/Q_downstream) * (1/d^2)
	}

	return(result)
}

#' @param ws A watershed
#' @param
interpSpace <- function(ws, values, sites, distMatrix) {

}
