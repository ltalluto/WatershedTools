context("Hydrology")
library("WatershedTools")

ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))

test_that("transport works with lsoda & euler", {
	initial <- lateral <- rep(0, nrow(ws$data))
	startloc <- 10546
	initial[startloc] <- 500

	## Note --> lsoda tests fail due to a bug in deSolve
	## "unlock_solver" not resolved from current namespace (deSolve)
	# expect_error(tra_lsoda <- transport(ws, initial, lateral, seq(0,600, 2), method='lsoda'), 
	# 	regex=NA)
	expect_error(tra_eul <- transport(ws, initial, lateral, seq(0,600, 2), method='euler'), 
		regex=NA)
	expect_equal(tra_eul[startloc,1], initial[startloc])
	expect_lt(tra_eul[startloc,2], tra_eul[startloc,1])
	# expect_equal(tra_lsoda[startloc,1], initial[startloc])
	# propdiff <- abs(tra_eul[startloc,2] - tra_lsoda[startloc,2]) / tra_lsoda[startloc,2]
	# expect_less_than(propdiff, 0.05)
})

test_that("Hydraulic geometry scaling", {
	ca <- 896 * 1000^2 ## 896 square kilometers
	expect_error(Q <- discharge_scaling(ca), regex=NA)
	expect_equal(Q, 14.6429, tolerance=1e-4)
	expect_error(Q2 <- discharge_scaling(ca, calib=list(A=100 * 1000^2, Q = 4.5)), regex=NA)
	expect_equal(Q2, 24.8893, tolerance = 1e-2)

	expect_error(hg <- hydraulic_geometry(Q), regex=NA)
	expect_equal(hg$discharge[1], Q)
	expect_equal(hg$velocity[1], 0.41683, tolerance = 1e-5)
	expect_equal(hg$depth[1], 0.89950, tolerance = 1e-5)
	expect_equal(hg$width[1], 40.25818, tolerance = 1e-5)
})
