library(sf)
library(raster)
riv = readRDS(system.file("testdata/ws_intersect/river.rds", package="WatershedTools"))
litho = readRDS(system.file("testdata/ws_intersect/lith.rds", package="WatershedTools"))
drain = raster(system.file("testdata/ws_intersect/drainage.tif", package="WatershedTools"))
ws = readRDS(system.file("testdata/ws_intersect/ws.rds", package="WatershedTools"))

test_that("ws_intersect", {
	expect_warning(sumtab <- ws_intersect(riv, list(litho = litho, landuse = litho), 
			area_id = list('xx', 'xx'), reach_id = 'a_cat_', drainage = drain, ws=ws,
			gisBase = "/Applications/GRASS-7.6.app/Contents/Resources"), regex='attribute')
	expect_true(all(c('xx', 'a_cat_', 'area', 'proportion', 'layer', 'method') %in% colnames(sumtab)))
	expect_true(all(sumtab[, .(sum=sum(proportion)), .(a_cat_, layer, method)]$sum == 1))
	
})


