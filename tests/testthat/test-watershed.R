ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))


test_that("Topology functions", {
	skip_on_cran()
	points =confluences(ws)[, 'id']
	expect_error(dm <- siteByPixel(ws, points), regex=NA)
	expect_equal(sum(dm, na.rm=TRUE), 3754613, tolerance = 0.001)
})

