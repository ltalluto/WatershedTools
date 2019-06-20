#' @name [.Watershed
#' @aliases names.Watershed
#' @aliases head.Watershed
#' @aliases summary.Watershed
#' @aliases plot.Watershed
#' 
#' @title S3 methods for watersheds
#' @param x A watershed object
#' @param i A vector of pixel IDs
#' @param j Variable names or column numbers to extract
#' @param variable Optional, name of variable to plot
#' @param transform A transformation function (e.g., `log`) to apply to the variable 
#' @param size Size of the plotted points
#' @param xlim Optional x-axis limits
#' @param ylim Optional y-axis limits
#' @param ... Additional parameters to pass to other methods 
#' 		(usually from [data.frame] or [SpatialPointsDataFrame])
#' @return A data.frame of the watershed data
#' @rdname methods
#' @export
`[.Watershed` <- function(x, i, j, ...) {
	if(missing(i)) {
		as.data.frame(x$data[,, ...])[,j]
	} else {
		as.data.frame(x$data)[match(i, x$data$id),j]
	}
}

#' @rdname methods
#' @return A vector of variable names
#' @export
names.Watershed <- function(x) {
	names(x$data)
}

#' @rdname methods
#' @return As [head.data.frame()]
#' @export
head.Watershed <- function(x, ...) head(x$data, ...)

#' @rdname methods
#' @return As [summary.data.frame()]
#' @export
summary.Watershed <- function(x, ...) summary(x$data, ...)

#' @rdname methods
#' @export
plot.Watershed <- function(x, variable, transform, size = 0.5, xlim, ylim, viridis_col = "inferno")
{
	if(!missing(transform))
		x$data[[variable]] <- transform(x$data[[variable]])
	pl <- ggplot2::ggplot(as.data.frame(x$data), ggplot2::aes(x=x, y=y))
	pl <- pl + ggplot2::geom_point(size=size)
	if(!missing(variable)){
		pl <- pl + ggplot2::aes_string(color=variable)
		if(requireNamespace("viridisLite")) {
			if(is.factor(x$data[[variable]])) {
				pl <- pl + ggplot2::scale_colour_viridis_d(option = viridis_col)
			} else {
				pl <- pl + ggplot2::scale_colour_viridis_c(option = viridis_col)
			}
		}
	}
	if(!missing(xlim))
		pl <- pl + ggplot2::xlim(xlim)
	if(!missing(ylim))
		pl <- pl + ggplot2::ylim(ylim)
	pl
}


## Dont forget assignment methods
## e.g., 
## `[<-.Watershed` <- function(x, i, j, value) {
	# x$data[i,j] <- value
	# x
## }
## also `$<-.Watershed`
## also add a subset function