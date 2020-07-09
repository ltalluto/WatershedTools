#' Intersect a river feature with a polygon feature
#' 
#' ws_intersect handles the use case where we want to summarize areal features, such as land use, soils, etc, across
#' reaches of a river, while taking into account drainage area, flow, etc.
#' 
#' There are three potential outputs, which will be identified in the output by the 'method' column
#' * `riparian`: The layers in `areas` are summarized within a riparian buffer surrounding each reach; 
#'  		provided if `rip_buffer` is not NA.
#' * `riparian_upstream`: The layers in `areas` are summarized within a riparian buffer for each focal reach and all
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
#' @param ws A watershed
#' @param ... Additional parameters to pass to [GrassSession()]
#' @return A data.table summarising the polygons in each layer in areas along the river provided in x.
#' @export
ws_intersect = function(x, areas, area_id, reach_id = "a_cat_", rip_buffer = 100, drainage = NA, ws, ...) {
	if(is(x, "Spatial"))
		x = sf::st_as_sf(x)
	if(!is(x, "sf"))
		stop("x must inherit from either sf or Spatial")
	res = .ws_intersect_sf(x, areas, area_id, reach_id, rip_buffer, drainage, ws, ...)
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
	
	if(is(drainage, "RasterLayer")) {
		dname = 'drainage'
		gs = GrassSession(drainage, layerName = dname, override = TRUE, ...)
		rids = unique(x[[reach_id]])
		reach_bottoms = lapply(rids, function(i) {
			ind = which(ws$data$vReachNumber == i)
			pt = as.integer(names(which(colSums(ws$adjacency[ind,ind, drop=F]) == 0)))
			as(ws$data[pt,], "SpatialPoints")
		})

		for(i in seq_len(rids)) {
			pt = reach_bottoms[[i]]
			rid = rids[i]
			cment = catchment(pt, dname, gs, output = "sf")
			cment = sf::st_transform(cment, sf::st_crs(areas[[1]]))
			cment_polys = mapply(function(xx, yy) sf::st_intersection(xx, yy), xx = list(cment), yy = areas, 
								 SIMPLIFY = FALSE)
			out = mapply(.sf_summary, x = cment_polys, by = area_id, layer = names(areas), 
						 MoreArgs = list(a_cat_ = rids[i], method = 'catchment'), SIMPLIFY = FALSE)
		
			if(!is.na(x_buffer)) {
				cment_polys_buff = mapply(function(xx, yy) sf::st_intersection(xx, yy),
										  xx = cment_polys, MoreArgs = list(yy=x_buffer), SIMPLIFY=FALSE)
				out = c(out, mapply(.sf_summary, x = cment_polys, by = area_id, layer = names(areas), 
									MoreArgs = list(a_cat_ = rid, method = 'catchment'), SIMPLIFY = FALSE))
			}
			res = c(res, out)
		}
		
		
		res = c(res, mapply(function(pt, rid) {
			cment = catchment(pt, dname, gs, output = "sf")
			## restore CRS, which gets corrupted a bit by grass
			cment = sf::st_transform(cment, sf::st_crs(areas[[1]]))
			cment_polys = mapply(function(xx, yy) sf::st_intersection(xx, yy), xx = list(cment), yy = areas, 
								 SIMPLIFY = FALSE)
			out = mapply(.sf_summary, x = cment_polys, by = area_id, layer = names(areas), 
					MoreArgs = list(a_cat_ = rids[i], method = 'catchment'), SIMPLIFY = FALSE)
			if(!is.na(x_buffer)) {
				cment_polys_buff = mapply(function(xx, yy) sf::st_intersection(xx, yy),
					xx = cment_polys, MoreArgs = list(yy=x_buffer), SIMPLIFY=FALSE)
				out = c(out, mapply(.sf_summary, x = cment_polys, by = area_id, layer = names(areas), 
					MoreArgs = list(a_cat_ = rid, method = 'catchment'), SIMPLIFY = FALSE))
			}
			out
		}, reach_bottoms, rids))
	}

	res = data.table::rbindlist(res)
	return(res)
}



.ws_intersect_grass = function(x, areas, area_id, reach_id, rip_buffer, drainage, ...) {
	stop("not done yet")
	if(is.na(drainage)) {
		ras = raster::raster(ext = raster::extent(x), crs = sf::st_crs(areas[[1]])$proj4string, resolution  = 25)
		layername = NA
	} else {
		ras = drainage
		layername = 'drainage'
	}
	gs = GrassSession(ras, layerName = layername, ...)
	rgrass7::writeVECT(as(x, "Spatial"), vname = "river", ignore.stderr=TRUE, v.in.ogr_flags="overwrite")
	area_name = 'ri_polygon'
	rgrass7::writeVECT(as(y, "Spatial"), vname = area_name, ignore.stderr=TRUE, v.in.ogr_flags="overwrite")
	buff_name = 'riv_buff'
	int_name = 'riv_isect'
	rgrass7::execGRASS("v.buffer", flags = c('t', 'overwrite', 'quiet'), input = 'river', output = buff_name, 
					   distance = rip_buffer, ignore.stderr=TRUE, Sys_ignore.stdout=TRUE)
	rgrass7::execGRASS("v.overlay", flags = c('overwrite', 'quiet'), ainput = buff_name, binput = area_name, output = int_name, 
					   operator='and', ignore.stderr=TRUE, Sys_ignore.stdout=TRUE)

	
}



#' Summarise an sf polygon feature by attributes
#' @details Additional arguments can be supplied to add columns to the output. For example, adding `category = 2` will add a column
#' named category to the output table, with a value of 2 for all entries. Useful for example when calling .sf_summary from mapply,
#' and it is desirable to know the combinations in mapply resulting in a given line of output
#' @param x An sf polygon
#' @param by Character vector, column names giving factors to summarise along. Partial matching with regular expressions
#' is allowed, as long as the result is unique.
#' @param subunit Optional character vector; If included, it will be used as a group variable(s) when computing proportions;
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
	colnames(areas)[grep(by[1], colnames(areas))] = "category"
	areas$category_colname = by[1]
	return(areas)
}
