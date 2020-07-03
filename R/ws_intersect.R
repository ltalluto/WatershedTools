#' Intersect a river feature with a polygon feature
#' 
#' ws_intersect handles the use case where we want to summarize areal features, such as land use, soils, etc, across
#' reaches of a river, while taking into account drainage area, flow, etc.
#' 
#' There are three potential outputs, which will be identified in the output by the 'method' column
#' * `riparian`: The layers in `areas` are summarized within a riparian buffer surrounding each reach; 
#'  		provided if `rip_buffer` is not NA.
#' * `riparian_upstream`: The layers in `areas` are summarized within a riparian buffer for each focual reach and all
#'          upstream reaches; provided if `rip_buffer` is not NA and `drainage` is not NA
#' * `catchment`: The layers in `areas` are summarized within the entire catchment for each reach, provided if
#'          `drainage` is not NA
#' 
#' @param x A river network; this should be a vector gis feature, either from the `sf` or `sp` packages
#' @param areas A named list of areal features to summarize; should be polygon features from `sf` or `sp`
#' @param reach_id Name of the field in x that identifies stream segments
#' @param rip_buffer The width of the riparian buffer; if NA, then no buffered analysis will be done, see 'details'
#' @param drainage Optional drainage direction raster for calculating catchment; see 'details'
#' @param use_sf Boolean, if TRUE, functions from the sf package will be used whenever possible, if FALSE, GRASS will be used
#' @param ... Additional parameters to pass to [GrassSession()]
#' @export
ws_intersect = function(x, areas, reach_id = "a_cat_", rip_buffer = 100, drainage = NA, use_sf = FALSE, ...) {
	if(is(x, "Spatial"))
		x = sf::st_as_sf(x)
	if(!is(x, "sf"))
		stop("x must inherit from either sf or Spatial")
	
	if(use_sf) {
		gs = NULL
	} else {
		if(is.na(drainage)) {
			ras = raster::raster(ext = raster::extent(x), crs = sp::proj4string(areas[[1]]), resolution  = 25)
			layername = NA
		} else {
			ras = drainage
			layername = 'drainage'
		}
		gs = GrassSession(ras, layerName = layername, ...)
	}
	
	# step 1: buffer and intersection
	if(!is.na(rip_buffer)) {
		rb_isect = mapply(.riv_buff_intersect, x = x, y = areas, moreArgs = list(rip_buffer = rip_buffer, gs = gs))
	}

	
	# for(i in seq_along(areas))
	# 	rgrass7::writeVECT(areas[[i]], vname = names(areas)[i])
}

#' Helper function for creating river buffers
#' @param x A river network; this should be a vector gis feature, either from the `sf` or `sp` packages
#' @param y The polygon feature to intersect, sf class
#' @param rip_buffer The width of the riparian buffer
#' @param gs A grass session, if missing or NA, then sf functions will be used
#' @return If gs is null or grass_output is "sf", then a buffer polygon layer of class "sf", if grass_output is "string",
#' then the name of the vector layer in the grass session
#' @keywords internal
.riv_buff_intersect = function(x, y, rip_buffer, gs) {
	if(missing(gs) || is.na(gs)) {
		x_buff = sf::st_buffer(x, rip_buffer)
		x_isect = sf::st_intersection(x_buff, y)
	} else {
		rgrass7::writeVECT(as(x, "Spatial"), vname = "river", ignore.stderr=TRUE, v.in.ogr_flags="overwrite")
		area_name = 'ri_polygon'
		rgrass7::writeVECT(as(y, "Spatial"), vname = area_name, ignore.stderr=TRUE, v.in.ogr_flags="overwrite")
		buff_name = 'riv_buff'
		int_name = 'riv_isect'
		rgrass7::execGRASS("v.buffer", flags = c('t', 'overwrite', 'quiet'), input = 'river', output = buff_name, 
						   distance = rip_buffer, ignore.stderr=TRUE, Sys_ignore.stdout=TRUE)
		rgrass7::execGRASS("v.overlay", flags = c('overwrite', 'quiet'), ainput = buff_name, binput = area_name, output = int_name, 
						   operator='and', ignore.stderr=TRUE, Sys_ignore.stdout=TRUE)
		x_isect = rgrass7::readVECT(int_name, ignore.stderr=TRUE)
		if(is(x_isect, "Spatial"))
			x_isect = sf::st_as_sf(x_isect)
	}
	return(x_isect)
}

