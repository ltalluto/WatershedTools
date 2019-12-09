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
	expect_error(flow <<- drainageAccumulation(filled$filledDEM, gisBase = gisBase), regex=NA)
})

test_that("Test the watershed delineation (raster object) workflow", {
	skip_on_cran()
	expect_error(streamRas <- extractStream(dem = filled$filledDEM,  
		accumulation = flow$accumulation, qthresh = 0.95, gisBase = gisBase), regex=NA)
	expect_equal(sum(!is.na(raster::values(streamRas))), 51915)
})


# test the workflow using GRASS layers instead of rasters
test_that("Test the watershed delineation (GRASS object) workflow", {
	skip_on_cran()
	gisBase <- "/Applications/GRASS-7.4.1.app/Contents/Resources/"
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	expect_error(gs <- GrassSession(testDEM, layerName = "dem", gisBase = gisBase), regex=NA)
	expect_error(gs <- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs),
		regex=NA)
	expect_error(gs <- fillDEM("demNotHere", filledDEM = "filledDEMNotHere", probs = "problems",
		gs = gs))
	expect_error(gs <- drainageAccumulation("filledDEM", accumulation = "accum", drainage = "drain",
		gs = gs), regex=NA)
	expect_error(gs <- extractStream(dem = "filledDEM", accumulation = "accum", qthresh = 0.95,
		outputName = "streamRas", gs = gs), regex=NA)
	expect_error(streamRas <- GSGetRaster("streamRas", gs), regex=NA)
	expect_equal(sum(!is.na(raster::values(streamRas))), 51915)
	ext <- raster::extent(testDEM)
	pt <- data.frame(x=mean(ext[1:2]), y=mean(ext[3:4]), elev=NA)
	sp::coordinates(pt) <- c(1,2)
	sp::proj4string(pt) <- sp::proj4string(testDEM)
	pt$elev <- raster::extract(testDEM, pt)
	expect_warning(ptsnap <- snapToStream(pt, streamRas, buff= 100), 
		regex="geographic coordinates")
	expect_equivalent(sp::coordinates(ptsnap), sp::coordinates(pt))
	expect_equal(ptsnap$elev, pt$elev)
	expect_warning(ptsnap <- snapToStream(pt, streamRas, buff= 500), 
		regex="geographic coordinates")
	expect_false(isTRUE(all.equal(sp::coordinates(ptsnap), sp::coordinates(pt))))
	expect_equal(ptsnap$elev, pt$elev)
})