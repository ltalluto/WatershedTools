context("Watershed Deliniation Functions")
library("WatershedTools")

gisBase <- getGISBase()

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
	testDEM <- raster::raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
	expect_error(gs <- GrassSession(testDEM, layerName = "dem", gisBase = gisBase), regex=NA)
	expect_error(gs <- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs),
				 regex=NA)
	expect_error(gs <- fillDEM("demNotHere", filledDEM = "filledDEMNotHere", probs = "problems",
							   gs = gs))
	expect_error(gs <- drainageAccumulation("filledDEM", accumulation = "accum", 
											drainage = "drain", gs = gs), regex=NA)
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

test_that("Stream order", {
	skip_on_cran()
	ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))
	expect_error(sord <- strahler(ws, parallel = FALSE), regex=NA)
	expect_error(sord2 <- strahler(ws, parallel = TRUE), regex=NA)

	# parallel should produce same results
	expect_identical(sord, sord2)

	# basic functionality
	expect_equal(range(sord), c(1,2))
	expect_equal(length(sord), nrow(ws$data))

	# make sure that breaking topology also breaks the function
	ind <- which(ws$data$reachID == 4)
	ws$adjacency[ind,] <- ws$adjacency[,ind] <- 0
	expect_error(sord <- strahler(ws), regex='topology')

})
