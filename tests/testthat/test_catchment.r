context("Watershed Deliniation Functions")
library("WatershedTools")

test_that("Test the catchment area workflow", {
	skip_on_cran()
	gisBase <- "/Applications/GRASS-7.4.1.app/Contents/Resources/"
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	gs <- GrassSession(testDEM, gisBase = gisBase)
	expect_error(filled <- fill_dem(testDEM, gisBase = gisBase), regex=NA)
	expect_error(flow <- watershed(filled, gisBase = gisBase), regex=NA)
	coords <- coordinates(flow)[which.max(values(flow[[1]])),, drop=FALSE]
	expect_error(catchArea <- catchment(coords, drain_name = "drainage_direction", 
		grass_session = gs), regex = NA)
	expect_equal(length(catchArea), 1)
	expect_gt(catchArea, 95716950)
	expect_lt(catchArea, 95716955)
})