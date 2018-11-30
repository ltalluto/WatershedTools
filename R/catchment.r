#' Produce a DEM with basins filled in, as well as a flow direction raster
#'
#' @param dem A [sp::SpatialGridDataFrame] or [raster::RasterLayer] object giving the digital
#' 		elevation model of the study area.
#' @param grass_session An existing grass session (optional)
#' @param dem_name Character; the name of the dem in the GRASS mapset to use; see details
#' @param file The file name of the raster to create; see details.
#' @param ... Additional parameters to pass to [GrassSession()]
#'
#' @details If `dem` is not missing, it will **overwrite** any existing DEM named `dem_name` in 
#' the GRASS mapset. If `dem` is missing, then any existing DEM will be used. 
#'
#' At the moment, the grass_session parameter is ignored; a new anonymous session will be created
#'
#' It is highly recommended to specify the `file` parameter (including the extension to specify
#' file format; e.g., .tif, .grd). If not specified, a temp file will be created and will be
#' lost at the end of the R session
#'
#' @return A [raster::raster], optionally written to `file`
fill_dem <- function(dem, grass_session, dem_name = "dem", file, ...)
{
	grass_session <- GrassSession(dem, layer_name = dem_name, ...)
	filled_dem <- "filled_dem"
	flow_direction <- "flow_direction"
	probs <- "problem_areas"
	rgrass7::execGRASS("r.fill.dir", flags=c("overwrite", "quiet"), input=dem_name, 
		output = filled_dem, direction = flow_direction, areas = probs)
	ras <- raster::raster(rgrass7::readRAST(filled_dem))
	if(!missing(file)) 
		ras <- raster::writeRaster(ras, filename = file)
	return(ras)
}

#' Watershed analysis tools
#'
#' @param dem A [sp::SpatialGridDataFrame] or [raster::RasterLayer] object giving the digital
#' 		elevation model of the study area; for best results, use the output from [fill_dem()].
#' @param dem_name Character; the name of the dem in the GRASS mapset to use; see details
#' @param threshold Minimum size of an exterior watershed, in *cells*
#' @param file The file name of the raster stack to create; see details.
#' @param ... Additional parameters to pass to [GrassSession()]
#' @details If `dem` is not missing, it will **overwrite** any existing DEM named `dem_name` in 
#' the GRASS mapset. If `dem` is missing, then any existing DEM will be used. 
#'
#' It is highly recommended to specify the `file` parameter (including the extension to specify
#' file format; e.g., .tif, .grd). If not specified, a temp file will be created and will be
#' lost at the end of the R session
#' @return A [raster::stack()] with two layers, flow accumulation ('accumulation') and drainage
#' 		direction ('drainage')
watershed <- function(dem, dem_name = "dem", threshold = 250, file, ...)
{
	grass_session <- GrassSession(dem, layer_name = dem_name, ...)
	accu_name <- "accumulation"
	drain_name <- "drainage_direction"
	rgrass7::execGRASS("r.watershed", flags=c("overwrite", "quiet"), elevation = dem_name, 
		accumulation = accu_name, drainage = drain_name)
	accu_ras <- raster::raster(rgrass7::readRAST(accu_name))
	drain_ras <- raster::raster(rgrass7::readRAST(drain_name))
	ras <- raster::stack(list(accumulation = accu_ras, drainage_direction = drain_ras))
	if(!missing(file)) 
		ras <- raster::writeRaster(ras, file)
	return(ras)
}

#' Catchment deliniation
#' 
#' @param x An [sp::SpatialPoints] objector a matrix of x-y coordinates
#' @param drainage Drainage direction raster, from e.g., [watershed]
#' @param drain_name Character; the name of the drainage raster in the GRASS mapset to use
#' @param grass_session An existing grass session (optional, recommended)
#' @param areas Logical; if `TRUE`, returns catchment area for each point, otherwise returns
#'  a [raster::RasterStack], one layer per point, with each layer delimiting the catchment for its
#'  respective point.
#' @details If `grass_session` is specified, the existing session will be used. In this case
#'   `drainage` is specified it will be written to the grass session with the name `drain_name`.
#'   If `drainage` is missing, then the grass session will be searched for an existing layer with
#'   a name given by `drain_name` (this is the fastest mode of operation).
#'
#' If no `grass_session` is specified, then a new session will be created, in which case `drainage`
#'   is required and will  be copied to the session.
#' 
#' This tool will perform much better if it is run with an existing grass session with 
#'	drainage direction already added:
#' 		`gs <- GrassSession(layer = drainage, layer_name = "drainage")`
#' 		`cAreas <- catchment(points, drain_name = "drainage", grass_session = gs, areas = TRUE)`
#' 	However note that one raster file is created *per point* in x, thus it is a good idea
#' 	with many points to specify area = TRUE to get area (rather than the raster) for each point
#' @return A vector of catchment areas (if `areas = TRUE`), otherwise a [raster::stack()] of 
#'   delineated catchments
catchment <- function(x, drainage, drain_name = "drainage_direction", grass_session, 
	areas = TRUE, file, ...)
{
	if(missing(grass_session)) {
		grass_session <- GrassSession(drainage, layer_name = drain_name, ...)
	} else {
		if(!missing(drainage))
			GSAddRaster(drain_name, drain_name, grass_session)
	}

	if(!is.matrix(x) & is.numeric(x)) {
		x <- matrix(x, ncol=2)
	} else if(!is.matrix(x)) {
		x <- sp::coordinates(x)
	}

	result <- if(areas) numeric(nrow(x)) else list()
	catch_name <- "catchment"
	for(i in 1:nrow(x))
	{
		rgrass7::execGRASS("r.water.outlet", flags=c("overwrite"), input = drain_name, 
			output = catch_name, coordinates = x[i,])
		if(areas) {
			result[i] <- catchmentArea(catch_name)
		} else {
			result <- c(result, ras)
		}
	}
	if(!areas) {
		result <- stack(result)
		if(!missing(file))
			result <- writeRaster(result, file)
	}
	return(result)
}


#' Compute the area of a delineated catchment
#' @param layer_name Character; the name of the grass catchment from which to compute area
#' @return Numeric; area of the catchment
#' @keywords internal
catchmentArea <- function(layer_name)
{
	vname <- paste0('v_', layer_name)

	rgrass7::execGRASS("r.to.vect", flags=c("overwrite", "quiet", "s"), 
		input = layer_name, output = vname, type='area', column = 'one')
	res <- rgrass7::execGRASS("v.to.db", flags=c("quiet", "p"), map = vname, 
				option = "area", intern=TRUE)
	idres <- as.numeric(sub("^(-?[0-9])+\\|(.+)", "\\1", res))
	res <- as.numeric(sub(".+?([0-9\\.]+)$", "\\1", res))
	val <- rgrass7::execGRASS("v.to.db", flags=c("quiet", "p"), map = vname, 
				option = "query", intern=TRUE, query_column="one")
	idval <- as.numeric(sub("^(-?[0-9])+\\|(.+)", "\\1", val))
	suppressWarnings(val <- as.numeric(sub(".+?([0-9\\.]+?)$", "\\1", val)))
	keep <- merge(data.frame(id = idres, area = res), 
		data.frame(id = idval, value = val), all.x = TRUE)

	## clean up
	rgrass7::execGRASS("g.remove", flags = c("f", "quiet"), type="vector", name=vname)
	sum(keep$area[keep$value == 1])
}