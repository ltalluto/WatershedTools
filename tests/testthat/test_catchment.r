context("Catchment Area")
library("WatershedTools")

test_that("Test the catchment area workflow with rasters", {
	skip_on_cran()
	gisBase <<- "/Applications/GRASS-7.4.1.app/Contents/Resources/"
	testDEM <<- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	gs <<- GrassSession(testDEM, layerName = "dem", gisBase = gisBase)
	gs <<- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs)
	gs <<- accumulate("filledDEM", accumulation = "accum", drainage = "drain", gs = gs)
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
	expect_error(streamCrop <- cropToCatchment(coords, streamRaster = stream$raster, 
		streamVector = stream$vector, drainage = "drain", gs = gs), regex=NA)
	vals <- raster::values(streamCrop$raster)
	expect_equal(sum(vals > 0, na.rm = T), 5777)
	expect_equal(sum(vals == 0, na.rm = T), 0)
	expect_equal(sum(vals > 0, na.rm = T) + sum(is.na(vals)), ncell(streamCrop$raster))
	vectPr <- spTransform(streamCrop$vector, sp::CRS("+init=epsg:32632"))
	expect_equal(as.integer(rgeos::gLength(vectPr)), 184657)
})

