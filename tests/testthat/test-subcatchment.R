
ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))
x <- confluences(ws)[1,'id']

test_that("Subcatchment extraction works as expected", {
	expect_error(subWs <- subcatchment(ws, x), regex=NA)
	expect_true("Watershed" %in% class(subWs))
	expect_equal(nrow(subWs$data), length(connect(ws, downstream=x, upstream=Inf)))
	expect_true(all(1:max(subWs[,'reachID']) %in% subWs[,'reachID']))
	expect_true(all(subWs[,'reachID'] %in% 1:max(subWs[,'reachID'])))
	expect_identical(ls(ws), ls(subWs))
})

