context("Watershed Object")
library("WatershedTools")

gisBase <- tryCatch(system2("grass76", args=c("--config", "path"), stdout=TRUE), 
	error = function(e) {
		tryCatch(system2("grass74", args=c("--config", "path"), stdout=TRUE), 
			error = function(e) "/Applications/GRASS-7.4.1.app/Contents/Resources/")
})

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
	expect_error(testWS <- Watershed(streamCrop, drainage), regex=NA)
})


test_that("Topology functions", {
	skip_on_cran()
	testWS <- readRDS(system.file("testdata/testWatershed.rds", package="WatershedTools"))
	points <- confluences(testWS)[1:10, 'id']
	expect_error(dm <- siteByPixel(testWS, points), regex=NA)
	expect_equal(sum(dm, na.rm=TRUE), 71.558, tolerance = 0.001)
})

