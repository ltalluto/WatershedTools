library(sf)
riv = readRDS(system.file("testdata/ws_intersect/river.rds", package="WatershedTools"))
litho = readRDS(system.file("testdata/ws_intersect/lith.rds", package="WatershedTools"))


test_that("ws_intersect", {
	expect_error(sumtab_sf <- ws_intersect(riv, list(litho = litho, landuse = litho), 
		area_id = list('xx', 'xx'), reach_id = 'a_cat_', use_sf=FALSE), regex=NA)

	expect_warning(sumtab_sf <- ws_intersect(riv, list(litho = litho, landuse = litho), 
			area_id = list('xx', 'xx'), reach_id = 'a_cat_', use_sf=TRUE), regex='attribute')
	expect_true(all(c('xx', 'a_cat_', 'area', 'proportion', 'layer', 'method') %in% colnames(sumtab_sf)))
	expect_true(all(sumtab_sf[, .(sum=sum(proportion)), .(a_cat_, layer)]$sum == 1))
	
})