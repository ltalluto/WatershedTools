#' Intersect a river feature with a polygon feature
#' 
#' ws_intersect handles the use case where we want to summarize areal features, such as land use, soils, etc, across
#' reaches of a river, while taking into account drainage area, flow, etc.
#' 
#' There are three output types:
#' * `riparian`: The layers in `areas` are summarized within a riparian buffer surrounding each reach
#' * `riparian_upstream`: The layers in `areas` are summarized within a riparian buffer for each focal reach and all
#'          upstream reaches
#' * `catchment`: The layers in `areas` are summarized within the entire catchment for each reach
#' 
#' @param ws A Watershed
#' @param areas A named list of areal features to summarize; either a [raster::stack()], or a list of
#' [sp::SpatialPolygons()] layers, or a list of [sf::sf()] layers.
#' @param area_id Character vector; the names of one or more columns in `areas` that should be summarised along
#' the river feature; partial matching with regular expressions is allowed, but each item in `area_id` must identify 
#' a single column in `areas`. Not required if areas is a raster stack.
#' @param rip_buffer The width of the riparian buffer
#' @param drainage Drainage direction raster for calculating catchment
#' @param ... Additional parameters to pass to [GrassSession()]
#' @return A data.table summarising the polygons in each layer in areas along the river provided in x.
#' @export
ws_intersect = function(ws, areas, area_id, rip_buffer = 50, drainage, ...) {
	if(!requireNamespace("fasterize", quietly = TRUE))
		stop("package 'fasterize' is required for this functionality")
	
	riv_r = raster::raster(ws$data)
	
	if(methods::is(areas, "RasterStack")) {
		if(!all(is.factor(areas))) {
			for(i in 1:raster::nlayers(areas))
				areas[[i]] = raster::ratify(areas[[i]])
		}
	} else if(methods::is(areas, "list")) {
		areas = mapply(.process_area, x = areas, id = area_id, MoreArgs = list(ras = riv_r), SIMPLIFY=FALSE)
		areas = raster::stack(areas)
	} else {
		areas = .process_area(areas, area_id, riv_r)
	}

	
	riv = as.sf.Watershed(ws)
	riv_buff = sf::st_sf(sf::st_union(sf::st_buffer(riv, rip_buffer)))
	riv_buff_r = fasterize::fasterize(riv_buff, raster = riv_r)

	dname = 'drainage'
	gs = GrassSession(drainage, layerName = dname, override = TRUE, ...)
	
	res = lapply(unique(ws$data$reachID), .do_ws_intersect, y = areas, ws = ws, 
				 width = rip_buffer, gs = gs, dname = dname, riv_buff = riv_buff_r)
	res = rbindlist(res)
	res$category = ""
	for(i in seq_len(raster::nlayers(areas))) {
		tab = raster::levels(areas)[[i]]
		if(missing(area_id)) {
			id = 1
		} else {
			if(area_id[[i]] %in% colnames(tab)) {
				id = which(colnames(tab) == area_id[[i]])
			} else {
				id = 1
			}
		}
		lyr = names(areas)[i]
		j = which(res[["layer"]] == lyr)
		if(length(j) > 0) {
			res[j, "category"] = tab[res[["value"]][j], id]
		}
	}
	res$value = NULL
	res
}


#' Take whatever format of areal layer, process into a raster for ws_intersect
#' @param x The areal layer to process
#' @param id The names of one or more columns in `x` that should be summarised along
#' @param ras A raster layer template to create rasters from polygons
.process_area = function(x, id, ras) {
	if(methods::is(x, "SpatialPolygons")) {
		x = sf::st_as_sf(x)
	}
	
	if(methods::is(x, "sf")) {
		warning("Vector input can produce inconsistent results; the preferred method is to
				provide a classified rasterlayer")
		fieldname = 'wst_category'
		x[[fieldname]] = factor(x[[id]])
		tab = levels(x[[fieldname]])
		x = fasterize::fasterize(x, raster = ras, field = fieldname)
		x = raster::ratify(x)
		lev = raster::levels(x)[[1]]
		lev[[id]] = tab[lev[,1]]
		levels(x) = list(lev)
	}
	
	if(!methods::is(x, "RasterLayer")) {
		stop("Please provide a raster (preferred), sp, or sf polygon object")
	}

	if(!raster::is.factor(x)) {
		x = raster::ratify(x)
	}
	return(x)
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
	i = which(ws$data$reachID == x)
	reach = ws$data[i,]
	adj = ws$adjacency[i,i,drop=FALSE]
	reach = .reach_to_sf(reach, adj)
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
