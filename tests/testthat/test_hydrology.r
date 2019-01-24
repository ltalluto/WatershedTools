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
	expect_less_than(tra_eul[startloc,2], tra_eul[startloc,1])
	# expect_equal(tra_lsoda[startloc,1], initial[startloc])
	# propdiff <- abs(tra_eul[startloc,2] - tra_lsoda[startloc,2]) / tra_lsoda[startloc,2]
	# expect_less_than(propdiff, 0.05)
})

