context("Catchment Area")
library("WatershedTools")

test_that("Test the catchment area workflow with rasters", {
	skip_on_cran()
	gisBase <- "/Applications/GRASS-7.4.1.app/Contents/Resources/"
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	gs <- GrassSession(testDEM, layerName = "dem", gisBase = gisBase)
	gs <- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs)
	gs <- watershed("filledDEM", accumulation = "accum", drainage = "drain", gs = gs)
	accum <- GSGetRaster("accum", gs)
	coords <- coordinates(accum)[which.max(values(accum)),, drop=FALSE]
	expect_error(catchArea <- catchment(coords, drainage = "drain", gs = gs), regex = NA)
	expect_equal(length(catchArea), nrow(coords))
	expect_gt(catchArea, 95563900)
	expect_lt(catchArea, 95563905)
})