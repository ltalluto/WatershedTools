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
#' @param area_id Character vector; the names of one or more columns in `areas` that should be summarised along
#' the river feature; partial matching with regular expressions is allowed, but each item in `area_id` must identify 
#' a single column in `areas` 
#' @param reach_id Name of the field in x that identifies stream segments
#' @param rip_buffer The width of the riparian buffer; if NA, then no buffered analysis will be done, see 'details'
#' @param drainage Optional drainage direction raster for calculating catchment; see 'details'
#' @param use_sf Boolean, if TRUE, functions from the sf package will be used whenever possible, if FALSE, GRASS will be used
#' @param ... Additional parameters to pass to [GrassSession()]
#' @return A data.table summarising the polygons in each layer in areas along the river provided in x.
#' @export
ws_intersect = function(x, areas, area_id, reach_id = "a_cat_", rip_buffer = 100, drainage = NA, use_sf = FALSE, ...) {
	if(is(x, "Spatial"))
		x = sf::st_as_sf(x)
	if(!is(x, "sf"))
		stop("x must inherit from either sf or Spatial")

	if(use_sf) {
		res = .ws_intersect_sf(x, areas, area_id, reach_id, rip_buffer, drainage)
	} else {
		res = .ws_intersect_grass(x, areas, area_id, reach_id, rip_buffer, drainage)
	}
	
	res
}


.ws_intersect_sf = function(x, areas, area_id, reach_id, rip_buffer, drainage, ws, ...) {
	res = list()
	x_buffer = NA
	if(!is.na(rip_buffer)) {
		x_buffer = sf::st_buffer(x, rip_buffer)
		x_isect = mapply(sf::st_intersection, list(x_buffer), areas, SIMPLIFY=FALSE)
		res = c(res, mapply(.sf_summary, 
				x = x_isect, by = area_id, layer = names(areas), 
				MoreArgs = list(subunit = reach_id, method = 'riparian'), SIMPLIFY = FALSE
		))
		x_buffer = sf::st_union(x_buffer)
	}
	
	if(!is.na(drainage)) {
		dname = 'drainage'
		gs = GrassSession(drainage, layerName = dname, override = TRUE, ...)
		rids = unique(x[[reach_id]])
		reach_bottoms = lapply(rids, function(i) {
			ind = which(ws$data$vReachNumber == rid)
			pt = as.integer(names(which(colSums(ws$adjacency[ind,ind]) == 0)))
			as(ws$data[pt,], "SpatialPoints")
		})
		for(i in length(rids)) {
			cment = catchment(reach_bottoms[[i]], dname, gs, areas = FALSE, vector = TRUE)
			cment_polys = mapply(function(xx, yy) sf::st_intersection(xx, yy), xx = list(cment), yy = areas, SIMPLIFY = FALSE)
			.... summarize
			if(!is.na(x_buffer)) {
				cment_polys_buff = mapply(function(xx, yy) sf::st_intersection(xx, yy), 
					xx = list(x_buffer), yy = cment_polys)
				.... summarize
				.... combine with cment_polys summary
			}
			.... rbind and return result (and convert to function and use lapply/mapply)
		}
		
	}
	
	res = data.table::rbindlist(res)
	return(res)
}



.ws_intersect_grass = function(x, areas, area_id, reach_id, rip_buffer, drainage) {
	stop("not done yet")
	if(is.na(drainage)) {
		ras = raster::raster(ext = raster::extent(x), crs = sf::st_crs(areas[[1]])$proj4string, resolution  = 25)
		layername = NA
	} else {
		ras = drainage
		layername = 'drainage'
	}
	gs = GrassSession(ras, layerName = layername, ...)
	res = list()
	res = c(res, mapply(.riv_buff_intersect, x = list(riv), y = areas, by = area_id, layer = names(areas),
						MoreArgs = list(rip_buffer = rip_buffer, gs = gs, 
										summarise = TRUE, subunit = reach_id, method = 'riparian'), SIMPLIFY = FALSE))
	res = data.table::rbindlist(res)
}




#' Helper function for creating river buffers, choosing either sf or grass
#' @param x A river network; this should be a vector gis feature, either from the `sf` or `sp` packages
#' @param rip_buffer The width of the riparian buffer
#' @param gs A grass session, if missing or NA, then sf functions will be used
#' @return A polygon layer of class 'sf' if gs is missing or NA, otherwise a character string giving the name
#' of the grass layer
#' @keywords internal
.riv_buff_intersect = function(x, y, rip_buffer, gs, summarise = FALSE, ...) {
	if(missing(gs) || is.na(gs)) {
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
	if(summarise) {
		return(.sf_summary(x_isect, ...))
	} else {
		return(x_isect)
	}
}

#' Summarise an sf polygon feature by attributes
#' @details Additional arguments can be supplied to add columns to the output. For example, adding `category = 2` will add a column
#' named category to the output table, with a value of 2 for all entries. Useful for example when calling .sf_summary from mapply,
#' and it is desirable to know the combinations in mapply resulting in a given line of output
#' @param x An sf polygon
#' @param by Character vector, column names giving factors to summarise along. Partial matching with regular expressions
#' is allowed, as long as the result is unique.
#' @param subunit Optional character vector; If included, it will be used as a grouble variable(s) when computing proportions;
#' variables here will also be used as grouping variables as in 'by'
#' @param ... Additional named arguments, see 'details'
#' @import data.table
#' 
#' @return A data.table with the following columns:
#' * ... ID columns, same names as `by`, giving the id levels of each summarised feature
#' @keywords internal
.sf_summary = function(x, by, subunit, ...) {
	args = list(...)
	if(!missing(subunit))
		by = c(by, subunit)
	cols = sapply(by, grep, x = colnames(x), simplify = TRUE)
	names = colnames(x)[cols]
	if(is.list(cols) || length(unlist(cols)) > ncol(x))
		stop("Non-unique column names given for by, please be more specific with column names")
	areas = cbind(as.data.table(x[,cols]), data.table(area = as.vector(sf::st_area(x))))
	areas = areas[, .(area=sum(area)), names]
	if(missing(subunit)) {
		areas[, proportion := area / sum(area)]
	} else {
		su_cols = sapply(subunit, grep, x = colnames(areas), simplify = TRUE)
		su_cols = colnames(areas)[su_cols]
		areas[, proportion := area / sum(area), su_cols]
	}
	for(nm in names(args)) {
		areas[, nm] = args[[nm]]
	}
	return(areas)
}
