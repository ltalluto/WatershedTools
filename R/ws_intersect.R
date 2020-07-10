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
#' @param ws A Watershed
#' @param areas A named list of areal features to summarize; should be polygon features from `sf` or `sp`
#' @param area_id Character vector; the names of one or more columns in `areas` that should be summarised along
#' the river feature; partial matching with regular expressions is allowed, but each item in `area_id` must identify 
#' a single column in `areas` 
#' @param rip_buffer The width of the riparian buffer
#' @param drainage Optional drainage direction raster for calculating catchment; see 'details'
#' @param ... Additional parameters to pass to [GrassSession()]
#' @return A data.table summarising the polygons in each layer in areas along the river provided in x.
#' @export
ws_intersect = function(ws, areas, area_id, rip_buffer = 50, drainage = NA, ...) {
	if(!requireNamespace("fasterize", quietly=TRUE)) {
		stop("The fasterize package is required for this function; \n",
			 "see https://github.com/ecohealthalliance/fasterize for instructions")
	}
	
	areas = mapply(function(a, id) {
		a$wst_category = factor(a[[id]])
		a
	}, areas, area_id, SIMPLIFY = FALSE)
	area_tables = mapply(function (a, id) levels(a$wst_category), areas, area_id, SIMPLIFY = FALSE)
	riv = raster::raster(ws$data, layer = which(names(ws$data) == "reachID"))
	riv_buff = raster::buffer(riv, width=rip_buffer)
	
	areas = lapply(areas, fasterize::fasterize, raster = riv, field = "wst_category")
	areas = raster::stack(areas)

	dname = 'drainage'
	gs = GrassSession(drainage, layerName = dname, override = TRUE, ...)
	
	res = lapply(unique(ws$data$reachID), .do_ws_intersect, y = areas, ws = ws, 
				 width = rip_buffer, gs = gs, dname = dname, riv_buff = riv_buff)
	res = rbindlist(res)
	res$category = ""
	for(lyr in names(area_tables))
		res$category[res$layer == lyr] = area_tables[[lyr]][res[layer == lyr, value]]
	res$value = NULL
	res
}

#' Helper function to perform the watershed intersection on a single reach
#' @param x A reachID to operate on
#' @param y A raster stack to intersect
#' @param ws A watershed
#' @param width Buffer width
#' @param gs A grass session
#' @param dname Drainage layer in the grass session
#' @param riv_buff A buffer layer for the entire river
.do_ws_intersect = function(x, y, ws, width, gs, dname, riv_buff) {
	reach = ws$data[ws$data$reachID == x,]
	reach = sf::st_as_sf(reach)
	rb = sf::st_buffer(reach, width)
	ras = raster::raster(ext=extent(as(rb, "Spatial")), res=raster::res(y))
	rb = fasterize::fasterize(rb, ras)
	yy = raster::crop(y, rb)
	rb = raster::resample(rb, yy)
	riparian = rb * yy
	names(riparian) = names(yy)

	outlet = outlets(ws, x, "Spatial")
	cment = catchment(outlet, dname, gs, output = "raster")
	cment_int = y * cment
	names(cment_int) = names(y)
	riparian_upstream = cment_int * riv_buff
	names(riparian_upstream) = names(cment_int)
	.wsi_summary(list(riparian = riparian, riparian_upstream=riparian_upstream, catchment = cment_int), reachID = x)
}


#' Summarise multiple raster layers by area
#' @name wsi_summary
#' @param ras A named list of raster stacks. The names of this list will be used for the 'method' column.
#' Layer names in each raster stack will give the layer name
#' @param ... Additional (single-valued) attributes to add to the output table
#' @import data.table
#' 
#' @return A data.table with the following columns:
#'  * method: The method used, either 'riparian', 'riparian_upstream', or 'catchment'
#'  * layer: The layer of origin for each row
#'  * value: The numeric value from the raster layer of each row
#'  * area: The area occupied by the category
#'  * proportion: The proportion of area (within the method and layer)
#' @keywords internal
#' 
.wsi_summary = function(ras, ...) {
	more_attr = list(...)
	tab = lapply(ras, .wsi_summary_single)
	tab = rbindlist(tab, idcol='method')
	if(!missing(more_attr)) {
		for(nm in names(more_attr)) {
			tab[[nm]] = more_attr[[nm]]
		}
	}
	return(tab)
}	

#' @rdname wsi_summary
#' @param ras A raster stack
#' @import data.table
.wsi_summary_single = function(ras) {
	tab = lapply(1:raster::nlayers(ras), function(i) table(raster::values(ras[[i]])))
	area = lapply(tab, function(x) x * prod(raster::res(ras)))
	prop = lapply(area, function(x) x / sum(x))
	rbindlist(mapply(function(a, p, lyr) data.table(layer = lyr, value = as.numeric(names(a)), 
			area = as.vector(a), proportion = as.vector(p)), area, prop, names(ras), SIMPLIFY=FALSE))
}
