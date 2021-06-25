#' Compute zonal statistics on a watershed
#' 
#' Compute zonal statistics, where the zones are determined by watershed reaches. 
#' 
#' There are three output types:
#' * `riparian`: The layers in `areas` are summarized within a riparian buffer surrounding each 
#' reach
#' * `riparian_upstream`: The layers in `areas` are summarized within a riparian buffer for each 
#' focal reach and all upstream reaches
#' * `catchment`: The layers in `areas` are summarized within the entire catchment for each reach
#' @param ws A Watershed]
#' @param ras A raster
#' @param statistics A list of function names; each function must accept a vector and return
#' 	a single value, and must have an na.rm argument
#' @param rip_buffer The width of the riparian buffer
#' @param drainage Drainage direction raster for calculating catchment
#' @param mc.cores Number of cores for parallel processing
#' @param catch_dir Where to store reach catchments
#' @param ... Additional parameters to pass to [GrassSession()]
#' @return A data.table summarising the values in each layer in areas along the river provided in x.
#' @export
ws_zonal = function(ws, ras, statistics = list("mean", "sd"), rip_buffer = 50, drainage, 
			mc.cores = parallel::detectCores(), catch_dir = tempdir(), ...) {
	if(!is(ras, "RasterLayer"))
		stop("Zonal works with a single RasterLayer")

	riv_buff_r = .ws_buffer(ws, rip_buffer)

	dname = 'drainage'
	gs = GrassSession(drainage, layerName = dname, override = TRUE, ...)

	rids = unique(ws$data$reachID)

	## slow step; compute the catchment for each reach
	## use a for loop so we can show progress
	catchments = list()
	pts = outlets(ws, rids, output = "Spatial")
	cat("Computing catchments, this may take a while\n")
	for(i in seq_along(rids)) {
		r = pts[i,]
		rn = rids[i]
		catchments[i] = catchment(r, dname, gs, 
			file = file.path(catch_dir, paste0("c_", rn, ".grd")) , output = "raster")
		cat("   ...", i, "/", length(rids), "(", floor(i/length(rids)*100), "%)...\r")
	}

	res = parallel::mclapply(rids, .do_ws_zonal, ws = ws, ras = ras, width = rip_buffer, gs = gs,
		drainage = dname, riv_buff = riv_buff_r, stats = statistics, mc.cores = mc.cores)
	res = lapply(res, function(x) {
		colnames(x) = statistics
		reshape2::melt(x)})
	res = lapply(res, data.table::as.data.table)
	names(res) = rids
	res = data.table::rbindlist(res, idcol = "reachID")
	colnames(res)[2:3] = c("zone", "statistic")
	res
}


#' Compute a zones for a single reach
#' 
#' @param ws A Watershed
#' @param rid A reachID
#' @param drainage Drainage direction raster, or name of raster from gs
#' @param gs A GrassSession()
#' @param width The width of the riparian buffer
#' @param riv_buff Buffer for the entire river
#' @return list of 3 rasters; the buffer for the reach, the catchment for the reach
#' intersected with the buffer for the whole river, and the catchment for the reach
#' @keywords internal
.reach_zone = function(ws, rid, drainage, gs, width, riv_buff) {
	pt = outlets(ws, rid, output = "Spatial")
	catch = catchment(pt, drainage, gs, output = "raster")
	rb = .ws_buffer(ws, width, rid)
	res = list(catchment = catch, reach_buffer = rb, 
		catchment_buffer = catch * riv_buff)
}


#' Compute zonal stats for a single reach
#' 
#' @param rid A reachID
#' @param ws A Watershed
#' @param ras The raster to compute zonal stats on
#' @param width The width of the riparian buffer
#' @param drainage Drainage direction raster, or name of raster from gs
#' @param gs A GrassSession()
#' @param riv_buff Buffer for the entire river
#' @param stats A list of statistics function names
.do_ws_zonal = function(rid, ws, ras, width, drainage, gs, riv_buff, stats) {
	tryCatch({
		zones = .reach_zone(ws, rid, drainage, gs, width, riv_buff)
		isect = lapply(zones, function(x) {
			if(!raster::compareRaster(x, ras, extent = FALSE, stopiffalse = FALSE))
				x = raster::resample(x, ras)
			x * ras})
		isect = raster::stack(isect)
		isect = raster::values(isect)
		isect = isect[rowSums(is.finite(isect)) > 0,]
		r = sapply(stats, function(x) apply(isect, 2, function(y) get(x)(y, na.rm = TRUE)))
		cat(rid, "done\n")
		r
	}, error = function(e) {
		msg = paste("rid", rid, "-", as.character(e))
		cat(msg, "\n")
		NA
	})
}


#' Compute a riparian buffer around a Watershed or a reach
#' @param ws A Watershed
#' @param buff The buffer width
#' @param rid A reach id, if not missing, the buffer will only include this reach
#' @return A raster
#' @keywords internal
.ws_buffer = function(ws, buff, rid) {
	if(!requireNamespace("fasterize", quietly = TRUE))
		stop("package 'fasterize' is required for this functionality")

	riv_r = raster::raster(ws$data)
	if(missing(rid)) {
		riv = as.sf.Watershed(ws)
		riv_buff = sf::st_sf(sf::st_union(sf::st_buffer(riv, buff)))
		res = fasterize::fasterize(riv_buff, raster = riv_r)
	} else {
		i = which(ws$data$reachID == rid)
		reach = ws$data[i,]
		## if the reach is only a single point, buffer the point instead
		# otherwise convert to linestring
		if(nrow(reach) > 1) {
			adj = ws$adjacency[i,i,drop=FALSE]
			reach = .reach_to_sf(reach, adj)
		} else {
			reach = sf::st_as_sf(reach)
		}
		rb = sf::st_buffer(reach, buff)
		ras = raster::raster(ext=extent(as(rb, "Spatial")), res=raster::res(riv_r))
		res = fasterize::fasterize(rb, ras)
	}
	res
}

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

	
	riv_buff_r = .ws_buffer(ws, rip_buffer)

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


	rb = .ws_buffer(ws, width, x)

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
