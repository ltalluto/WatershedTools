#' Provide interpolation in time and space
#' @param ws A watershed
#' 
#' @details The data are expected in 'tall' format, so `sites`, `times`, and `values` must all
#' have the same length.
#' 
#' It is mandatory to specify at least one of `sites_out` or `times_out`.
#' 
#' @param sites A vector of pixel IDs where measurements have been performed
#' @param times An optional vector of times (e.g., in minutes) at which measurements were made.
#' @param values A vector of measured values to interpolate
#' @param sites_out Either missing, NA, or a vector of pixelIDs. See "details." 
#' @param times_out Either `missing, NA, or a vector of times. See "details" 
interpolate <- function(ws, sites, times, values, sites_out = NA, times_out) {
	if(!(length(sites) == length(times) & length(sites) == length(values)))
		stop("sites, times, and values must all have the same length")
}