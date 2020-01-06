context("Watershed Topology")
library("WatershedTools")

ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))
outletID <- outlets(ws)$id
hwID <- headwaters(ws)[1,'id']

test_that("Accumulation works", {
	expect_error(accum_outlet_us <- 
				 	accumulate(ws, upstream = Inf, downstream = outletID, direction = "up"), regex=NA)
	expect_error(accum_outlet_us_parallel <- 
				 	accumulate(ws, upstream = Inf, downstream = outletID, direction = "up", parallel=TRUE),
				 regex=NA)
	expect_error(accum_outlet_ds <- 
				 	accumulate(ws, upstream = Inf, downstream = outletID, direction = "down"), regex=NA)
	expect_error(accum_hw_ds <- 
				 	accumulate(ws, upstream = hwID, downstream = Inf, direction = "down"), regex=NA)
	expect_error(accum_hw_us <- 
				 	accumulate(ws, upstream = hwID, downstream = Inf, direction = "up"), regex=NA)
	expect_error(accum_both_us <- 
				 	accumulate(ws, upstream = Inf, downstream = c(outletID, hwID), direction = "down"), 
				 regex=NA)
	expect_error(accum_out_hw_ds <- 
				 	accumulate(ws, upstream = hwID, downstream = outletID, direction = "down"), regex=NA)
	expect_error(accum_out_hw_us <- 
				 	accumulate(ws, upstream = hwID, downstream = outletID, direction = "up"), regex=NA)
	
	# results shouldn't depend on parallel computing or not
	expect_identical(accum_outlet_us, accum_outlet_us_parallel)
	
	# from outlet to all headwaters should recover all points
	expect_identical(sort(accum_outlet_us[,1]), sort(ws[,'id']))
	expect_identical(sort(accum_outlet_ds[,1]), sort(ws[,'id']))
	
	# from hw to outlet should produce same pixels no matter which direction
	expect_identical(sort(accum_hw_ds[,1]), sort(accum_hw_us[,1]))
	expect_identical(sort(accum_out_hw_ds[,1]), sort(accum_out_hw_us[,1]))
	expect_identical(accum_out_hw_ds, accum_hw_ds)
	expect_identical(accum_out_hw_us, accum_hw_us)
	
	## upstream always negative, downstream always positive
	expect_gt(0, sum(accum_outlet_us[,2]))
	expect_gt(0, sum(accum_hw_us[,2]))
	expect_gt(sum(accum_outlet_ds[,2]), 0)
	expect_gt(sum(accum_hw_ds[,2]), 0)
	
})

test_that("Watershed Distance", {
	expect_error(dmat <- wsDistance(ws, c(outletID, hwID)), regex=NA)
	expect_equal(dmat[1, outletID], 0)
	expect_equal(dmat[2, hwID], 0)
	expect_equal(dmat[1, hwID], -1*dmat[2, outletID])
	expect_true(all(dmat[1,] <= 0))
	expect_equal(ncol(dmat), nrow(ws$data))
})

test_that("Site By Reach", {
	expect_error(sbyr <- siteByReach(ws, c(outletID, hwID)), regex=NA)
	expect_equal(ncol(sbyr), max(ws[,'reachID']))
	# outlet is connected to all reaches, headwater to itself only
	expect_true(all(sbyr[1,] == 1))
	expect_equal(sum(sbyr[2,]), 1)
	expect_equal(unname(which(sbyr[2,] == 1)), ws[hwID, 'reachID'])
})

test_that("Downstream Neighbor", {
	ussite <- which(ws$adjacency[outletID,]==1)
	expect_error(ndsnb <- nearestDownstreamNeighbor(ws, c(outletID, ussite, hwID)), regex=NA)
	expect_true(!outletID %in% ndsnb[,1])
	expect_equal(unname(ndsnb[match(ussite, ndsnb[,1]),2]), outletID)
	expect_equal(unname(ndsnb[match(hwID, ndsnb[,1]),2]), unname(ussite))
})



# test_that("Nearest Neighbors", {
# 	expect_error(nnsnb <- nearestNeighbors(ws, c(outletID, ussite, hwID), 
# 		sites=c(outletID, ussite, hwID)), regex=NA)
# })

