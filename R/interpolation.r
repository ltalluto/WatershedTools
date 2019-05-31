#' Provide interpolation in time
#' @param ws A watershed
#' @param values A matrix, with sites in rows and times in columns.
#' @param times_out A vector of target times for interpolation
#' @param sites An optional vector of sites with length = nrow(values).
#' @param times An optional vector of times at which observations have been made.
#' @param ... Additional parameters to pass to approx
interpolate <- function(ws, values, times_out, sites = as.integer(rownames(values)), 
		times = as.integer(colnames(values)), ...) {
	res <- t(apply(values, 1, function(y) approx(x = times, y = y, xout=times_out, ...)))
	rownames(res) <- sites
	colnames(res) <- times_out
	res
}

