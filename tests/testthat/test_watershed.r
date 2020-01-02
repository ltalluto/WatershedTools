context("Watershed Object")
library("WatershedTools")

getGISBase <- function(appPath = "/Applications", suffix = "Contents/Resources") {
	grassBase <- list.files(appPath, pattern='GRASS')
	gVersion <- as.numeric(sub(".*(7\\.[0-9\\.]+)\\.app", "\\1", grassBase))
	if(length(grassBase) > 1)
		grassBase <- grassBase[which.max(gVersion)]
	file.path(appPath, grassBase, suffix)
}
gisBase <- getGISBase()

test_that("Creation of a basic Watershed proceeds without error", {
	skip_on_cran()
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	gs <- GrassSession(testDEM, layerName = "dem", gisBase = gisBase)
	gs <- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs)
	gs <- drainageAccumulation("filledDEM", accumulation = "accum", drainage = "drain",
		gs = gs)
	gs <- extractStream(dem = "filledDEM", accumulation = "accum", qthresh = 0.95,
		outputName = "streamRas", gs = gs)
	streamRas <- GSGetRaster("streamRas", gs)
	drainage <- GSGetRaster("drain", gs)
	accum <- GSGetRaster("accum", gs)
	coords <- sp::coordinates(accum)[which.max(raster::values(accum)),, drop=FALSE]
	streamCrop <- cropToCatchment(coords, streamRaster = streamRas, drainage = "drain", gs = gs)
	expect_error(testWS <<- Watershed(streamCrop, drainage), regex=NA)
})

test_that("All created subobjects have correct types", {
	skip_on_cran()
	objNames <- c('adjacency', 'data', 'reach_adjacency', 'reach_connectivity')
	objClasses <- list(adjacency = "dgCMatrix", 
		data = "SpatialPixelsDataFrame", 
		reach_adjacency = "ngCMatrix",
		reach_connectivity = "dgCMatrix")
	expect_true(all(objNames %in% ls(testWS)))
	expect_true(objClasses$adjacency %in% class(testWS$adjacency))
	expect_true(objClasses$data %in% class(testWS$data))
	expect_true(objClasses$reach_adjacency %in% class(testWS$reach_adjacency))
	expect_true(objClasses$reach_connectivity %in% class(testWS$reach_connectivity))
})

test_that("Topology functions", {
	skip_on_cran()
	# testWS <- readRDS(system.file("testdata/testWatershed.rds", package="WatershedTools"))
	points <- confluences(testWS)[1:10, 'id']
	expect_error(dm <- siteByPixel(testWS, points), regex=NA)
	expect_equal(sum(dm, na.rm=TRUE), 71.558, tolerance = 0.001)
})

