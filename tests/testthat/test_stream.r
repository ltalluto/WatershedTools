context("Watershed Deliniation Functions")
library("WatershedTools")

gisBase <- "/Applications/GRASS-7.4.1.app/Contents/Resources/"

test_that("Fill DEM (raster object)", {
	skip_on_cran()
	testDEM <<- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	expect_error(filled <<- fillDEM(testDEM, gisBase = gisBase), regex=NA)
})

test_that("Flow accumulation (rasters)", {
	skip_on_cran()
	expect_error(flow <<- watershed(filled$filledDEM, gisBase = gisBase), regex=NA)
})

test_that("Test the watershed delineation (raster object) workflow", {
	skip_on_cran()
	expect_error(streamRas <- extractStream(dem = filled$filledDEM,  
		accumulation = flow$accumulation, qthresh = 0.95,gisBase = gisBase), regex=NA)
	expect_equal(sum(!is.na(values(streamRas))), 51915)
})


## test the workflow using GRASS layers instead of rasters