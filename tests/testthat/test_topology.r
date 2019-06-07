context("Watershed Topology")
library("WatershedTools")

ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))

test_that("Accumulation works", {
	outletID <- outlets(ws)$id
	hwID <- headwaters(ws)[1,'id']
	expect_error(accum_outlet_us <- 
		accumulate(ws, upstream = Inf, downstream = outletID, direction = "up"), regex=NA)
	expect_error(accum_outlet_ds <- 
		accumulate(ws, upstream = Inf, downstream = outletID, direction = "down"), regex=NA)
	expect_error(accum_hw_ds <- 
		accumulate(ws, upstream = hwID, downstream = Inf, direction = "down"), regex=NA)
	expect_error(accum_hw_us <- 
		accumulate(ws, upstream = hwID, downstream = Inf, direction = "up"), regex=NA)
	expect_error(accum_out_hw_ds <- 
		accumulate(ws, upstream = hwID, downstream = outletID, direction = "down"), regex=NA)
	expect_error(accum_out_hw_us <- 
		accumulate(ws, upstream = hwID, downstream = outletID, direction = "up"), regex=NA)

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

# setwd("..")
# library(devtools)
# load_all()

# library(ggplot2)

# # ws <- readRDS("~/work/projects/metabolism/catchment_delineations/vjosa/res/vjosaWatershedSpring2018.rds")
# ws <- readRDS("~/work/projects/catchment_delineations/vjosa/res/vjosaWatershedSpring2018.rds")

# xx <- c(1369, 22506, 13792, 8, 44100, 38600)
# pl <- plot(ws)
# for(i in x)
# 	pl <- pl + geom_point(x = ws[i,'x'], y=ws[i,'y'], colour = 'red')
# pl

# distMatrix <- wsDistance(ws, xx)
# x <- c(897, 43207, 6789)

# test <- nearestNeighbors(ws, x, distMatrix)
# test2 <- nearestNeighbors(ws, x, sites=xx)



# system.time(test1 <- accumulateOne(ws, x, variable = 'length')[y])
# system.time(test2 <- {test2 <- accumulate(ws, x, y); test2[test2[,1] == y, 2]})
# c(test1, test2)

# system.time(test3 <- accumulateOne(ws, x, variable = 'length'))
# system.time(test4 <- accumulate(ws, x, Inf))
# c(test3[1], test4[test4[,1] == 1,2])

# system.time(test5 <- accumulateOne(ws, y, variable = 'length'))
# system.time(test6 <- accumulate(ws, Inf, y, parallel = TRUE))
# c(test5[1,2], test6[1,2])

# test[y]


# ws$data$test <- test
# plot(ws, 'test', transform = function(x) log(abs(x)))


# TESTS TO WRITE
# 		accumulate
#		wsDistance
#		siteByReach
#		nearestDownstreamNeighbor
#		downstreamDist
#		nearestNeighbors
