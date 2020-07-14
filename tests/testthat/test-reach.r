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

ws = readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))
test_that("Reach resizing works", {
	skip_on_cran()
	expect_error(ws_rsz <- resize_reaches(ws, 500, 200), regex=NA)
	## parallel or not shouldn't matter
	expect_identical(ws_rsz, resize_reaches(ws, 500, 200, parallel=FALSE))
	
	rch_len = tapply(ws_rsz$data$length, ws_rsz$data$reachID, sum)
	expect_gt(length(unique(ws_rsz$data$reachID)), length(unique(ws$data$reachID)))
	expect_equal(nrow(ws_rsz$reach_adjacency), length(unique(ws_rsz$data$reachID)))
	expect_gte(min(rch_len), 200)
	# max size is the size + min_size + 2sqrt(2)res(ws)
	expect_lte(max(rch_len), 500+200+28)
})

test_that("Reach Adjacency", {
	hw = headwaters(ws)$reachID
	expect_error(adj <- sapply(hw, function(i) WatershedTools:::reachAdj(ws, i)), regex=NA)
	expect_true(all(sapply(adj, is.null)))
	rch = unique(ws$data$reachID)
	rch = rch[!rch %in% hw]
	expect_error(adj <- sapply(rch, function(i) WatershedTools:::reachAdj(ws, i)), regex=NA)
	expect_gte(nrow(adj), length(rch))
})