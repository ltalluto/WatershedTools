test_that("Reach trimming works", {
	skip_on_cran()
	ws = readRDS(system.file("testdata/ws_intersect/ws.rds", package="WatershedTools"))
	expect_error(ws_trim <- trim_reaches(ws, 200), regex=NA)
	expect_lt(length(unique(ws_trim$data$reachID)), length(unique(ws$data$reachID)))
	expect_equal(nrow(ws_trim$reach_adjacency), length(unique(ws_trim$data$reachID)))
	
	# check headwater reach lengths
	rch_len = tapply(ws_trim$data$length, ws_trim$data$reachID, sum)
	hw = headwaters(ws_trim)$reachID
	rch_len = rch_len[names(rch_len) %in% as.character(hw)]
	expect_gte(min(rch_len), 200)
})

