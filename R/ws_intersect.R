#' Intersect a river feature with a polygon feature
#' 
#' `Watershed`-friendly wrapper to [watershed::w_intersect()], see the help file in the 
#' `watershed` package for details on the function.
#' @param ws A Watershed
#' @param ... Additional parameters for [watershed::w_intersect()]
#' @return A data.table summarising the polygons in each layer in areas along the river provided in x.
#' @export
ws_intersect = function(ws, ...) {
	if(!requireNamespace("watershed")) {
		stop("The companion package 'watershed' is required for this analysis. Please install it ",
			"using devtools::install_github('flee-group/watershed'). ",
			"See https://github.com/flee-group/watershed for more information.")
	}
	riv_r = raster::raster(ws$data)
	riv = as.sf.Watershed(ws)

	args = list(...)
	args$riv = riv
	# get the outlet of each reach
	if(!pts %in% names(args)) {
		out = outlets(ws, unique(ws[, 'reachID']), "Spatial")
		ags$pts = sf::st_as_sf(out)
	}

	do.call(watershed::w_intersect, args)
}



