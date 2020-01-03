context("Hydrology")
library("WatershedTools")

ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))
initial <- lateral <- rep(0, nrow(ws$data))
startloc <- headwaters(ws)[1, 'id']
initial[startloc] <- 500


test_that("transport works with lsoda & euler", {

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
	expect_error(Q <- discharge_scaling(ca, calib=data.frame(A=100 * 1000^2, Q = 4.5)), regex=NA)
	expect_equal(Q, 24.8893, tolerance = 1e-2)

	expect_error(hg <- hydraulic_geometry(Q), regex=NA)
	expect_equal(hg$discharge[1], Q)
	expect_equal(hg$velocity[1], 0.4848648, tolerance = 1e-5)
	expect_equal(hg$depth[1], 1.051318, tolerance = 1e-5)
	expect_equal(hg$width[1], 50.3857, tolerance = 1e-5)
})

test_that("Simple transport-reaction model", {
	react <- function(t, y, rate) y * rate
	rt1 <- -0.5
	rt3 <- 1
	expect_error(res <- transport(ws, initial, lateral, seq(0,50, 1), method = 'euler', 
		rxn = react, rxnParams = list(rate = rt1)), regex=NA)
	expect_error(res2 <- transport(ws, initial, lateral, seq(0,50, 1), method = 'euler', 
		rxn = react, rxnParams = list(rate = 0)), regex=NA)
	expect_error(res3 <- transport(ws, initial, lateral, seq(0,50, 1), method = 'euler', 
		rxn = react, rxnParams = list(rate = rt3)), regex=NA)
	tr <- res2[startloc,1] - res2[startloc,2]
	expect_equal(res[startloc,1] - tr, res2[startloc,2])
	expect_equal((res[startloc,1] - tr) + (res[startloc,1] - tr) * rt1, res[startloc,2])
	expect_equal((res3[startloc,1] - tr) + (res3[startloc,1] - tr) * rt3, res3[startloc,2])

})

test_that("Q from flow model is sensible", {
	# single measurement
	# must have same length for depth and velocity
	expect_error(WatershedTools:::q_from_flow(depth = rep(1,5), velocity = rep(1,3), width=1), 
		regex = "length")
	expect_equal(WatershedTools:::q_from_flow(depth = rep(1,5), velocity = rep(1,5), width=1)$Q, 1)

	# multiple measurememts
	expect_equal(WatershedTools:::q_from_flow(depth = rep(1,5), velocity = c(1,1,1,2,2), 
		flowID = c(1,1,1,2,2), width=c(1,1), widthID = c(1,2))$Q, c(1,2))
})
