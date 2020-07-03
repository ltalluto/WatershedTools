library(sf)
riv = readRDS(system.file("testdata/ws_intersect/river.rds", package="WatershedTools"))
litho = readRDS(system.file("testdata/ws_intersect/lith.rds", package="WatershedTools"))

test_that("Buffering & Intersection works", {
	skip_on_cran()
	ras = raster::raster(ext = raster::extent(riv), crs = sf::st_crs(riv)$proj4string, resolution  = 25)
	layername = NA
	gs = GrassSession(ras, layerName = layername, gisBase = WatershedTools:::getGISBase())

	expect_warning(isect_sf <- WatershedTools:::.riv_buff_intersect(riv, litho, 100, gs = NA), regex="attribute")
	expect_error(isect_gr <- WatershedTools:::.riv_buff_intersect(riv, litho, 100, gs = gs), regex=NA)
	
	# these methods are really different, but at least expect similar areas for the intersection fields
	# expect less than 1% error in area on average
	s1 = reshape2::melt(tapply(st_area(isect_sf), list(isect_sf$xx, isect_sf$a_cat_), sum), na.rm=TRUE)
	s2 = reshape2::melt(tapply(st_area(isect_gr), list(isect_gr$b_xx, isect_gr$a_a_cat_), sum), na.rm=TRUE)
	s_merge = merge(s1, s2, by = c("Var1", "Var2"))
	expect_lt(mean(abs(s_merge$value.x - s_merge$value.y) / s_merge$value.x), 0.01)
	
})
