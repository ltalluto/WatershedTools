library("WatershedTools")

tryCatch(rgrass7::use_sp(), error = function (e) NULL)

test_that("Test the catchment area workflow with rasters", {
	skip_on_cran()
	gisBase <<- WatershedTools:::getGISBase()
	testDEM <<- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	gs <<- GrassSession(testDEM, layerName = "dem", gisBase = gisBase)
	gs <<- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs)
	gs <<- drainageAccumulation("filledDEM", accumulation = "accum", drainage = "drain", gs = gs)
	accum <<- GSGetRaster("accum", gs)
	coords <<- sp::coordinates(accum)[which.max(raster::values(accum)),, drop=FALSE]
	expect_error(catchArea <- catchment(coords, drainage = "drain", gs = gs), regex = NA)
	expect_equal(length(catchArea), nrow(coords))
	expect_gt(catchArea, 95563900)
	expect_lt(catchArea, 95563905)
})

test_that("Crop to catchment", {
	skip_on_cran()
	stream <<- extractStream(dem = "filledDEM", accumulation = accum, qthresh = 0.95, 
							 gs = gs, type='both')
	# sitesSnap <- snapToStream(sites, thurStream$raster, buff= 400)
	# this line produces a text conversion warning
	expect_warning(streamCrop <- cropToCatchment(coords, streamRaster = stream$raster, 
												 streamVector = stream$vector, drainage = "drain", gs = gs))
	vals <- raster::values(streamCrop$raster)
	expect_equal(sum(vals > 0, na.rm = T), 5777)
	expect_equal(sum(vals == 0, na.rm = T), 0)
	expect_equal(sum(vals > 0, na.rm = T) + sum(is.na(vals)), raster::ncell(streamCrop$raster))
})


