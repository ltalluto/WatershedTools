context("GRASS GIS Interface")
library("WatershedTools")

tryCatch(rgrass7::use_sp(), error = function (e) NULL)

test_that("Anonymous Grass session can be started", {
	skip_on_cran()
	gisBase <- getGISBase()
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	expect_error(gs <- GrassSession(testDEM, gisBase = gisBase), regex = NA)
	err <- "A GRASS location is already in use"
	expect_error(gs2 <- GrassSession(testDEM, override = FALSE, gisBase = gisBase), regex = err)
	expect_error(gs <- GSAddRaster(testDEM, "dem", gs), regex = NA)
	expect_match(gs$layers, "dem", all=FALSE)
})
