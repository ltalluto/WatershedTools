

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
