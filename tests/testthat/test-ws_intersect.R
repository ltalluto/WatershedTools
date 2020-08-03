library(sf)
library(raster)
library(WatershedTools)

riv = readRDS(system.file("testdata/ws_intersect/river.rds", package="WatershedTools"))
litho = readRDS(system.file("testdata/ws_intersect/lith.rds", package="WatershedTools"))
landuse = readRDS(system.file("testdata/ws_intersect/landuse.rds", package="WatershedTools"))
drain = raster(system.file("testdata/ws_intersect/drainage.tif", package="WatershedTools"))
ws = readRDS(system.file("testdata/ws_intersect/ws.rds", package="WatershedTools"))
ws = subcatchment(ws, 29000)

test_that("ws_intersect", {
	skip_on_cran()
	expect_warning(sumtab <- ws_intersect(ws, list(litho=litho, landuse=landuse), area_id = list('xx', 'code_18'), 
			rip_buffer=50, drainage = drain, gisBase = "/Applications/GRASS-7.6.app/Contents/Resources"), regex='extent')

	## try again, this time using raster layers
	expect_warning(areas <- mapply(WatershedTools:::.process_area, list(litho=litho, landuse=landuse), 
			list('xx', 'code_18'), MoreArgs = list(ras=drain)), regex="Vector")
	areas = stack(areas)
	expect_warning(sumtab2 <- ws_intersect(ws, areas,  rip_buffer=50, drainage = drain, 
			gisBase = "/Applications/GRASS-7.6.app/Contents/Resources"), regex='extent')
	

	expect_true(all(c('method', 'layer', 'area', 'proportion', 'reachID', 'category') %in% colnames(sumtab)))
	expect_true(all(unique(sumtab$category) %in% c(unique(landuse$code_18), unique(litho$xx))))
	expect_equal(sum(sumtab[, .(sum=sum(proportion) - 1), .(method, layer, reachID)]$sum),0)
})


