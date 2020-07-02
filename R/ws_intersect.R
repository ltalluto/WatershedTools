#' Intersect a river feature with a polygon feature
#' 
#' ws_intersect handles the use case where we want to summarize areal features, such as land use, soils, etc, across
#' reaches of a river, while taking into account drainage area, flow, etc.
#' 
#' There are three potential outputs, which will be identified in the output by the 'method' column
#' * `riparian`: The layers in `areas` are summarized within a riparian buffer surrounding each reach; 
#'  		provided if `rip_buffer` is not NA.
#' * `riparian_upstream`: The layers in `areas` are summarized within a riparian buffer for each focual reach and all
#'          upstream reaches; provided if `rip_buffer` is not NA and `catchment` is TRUE
#' * `catchment`: The layers in `areas` are summarized within the entire catchment for each reach, provided if
#'          `catchment` is TRUE
#' 
#' @param x A river network; this should be a vector gis feature, either from the `sf` or `sp` packages
#' @param areas A named list of areal features to summarize; should be polygon features from `sf` or `sp`
#' @param reach_id Name of the field in x that identifies stream segments
#' @param rip_buffer The width of the riparian buffer; if NA, then no buffered analysis will be done, see 'details'
#' @param catchment Boolean, if TRUE catchment area-based analyses will be done; see 'details'
#' @param use_sf Boolean, if TRUE, functions from the sf package will be used whenever possible, if FALSE, GRASS will be used
#' @param ... Additional parameters to pass to [GrassSession()]
ws_intersect = function(x, areas, reach_id = "a_cat_", rip_buffer = 100, catchment = FALSE, use_sf = FALSE, ...) {
	if(is(x, "Spatial"))
		x = sf::st_as_sf(x)
	if(!is(x, "sf"))
		stop("x must inherit from either sf or Spatial")
	
	# step 1: intersection
	if(!is.na(rip_buffer)) {
		if(use_sf) {
			x_buff = sf::st_buffer(x, rip_buffer)
		} else {
			gs = GrassSession(as(areas[[1]], "Spatial"), ...)
		}
	}
}
