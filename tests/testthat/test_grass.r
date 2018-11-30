context("GRASS GIS Interface")
library("WatershedTools")

test_that("Anonymous Grass session can be started", {
	skip_on_cran()
	gisBase <- "/Applications/GRASS-7.4.1.app/Contents/Resources/"
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	expect_error(gs <- GrassSession(testDEM, gisBase = gisBase), regex = NA)
	err <- "A GRASS location is already in use"
	expect_error(gs2 <- GrassSession(testDEM, override = FALSE), regex = err)
	expect_error(gs <- GSAddRaster(testDEM, "dem", gs), regex = NA)
	expect_match(gs$layers, "dem", all=FALSE)
})
