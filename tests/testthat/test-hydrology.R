
test_that("Q from flow model is sensible", {
	# single measurement
	# must have same length for depth and velocity
	expect_error(q_from_flow(depth = rep(1,5), velocity = rep(1,3), width=1), 
				 regex = "length")
	expect_equal(q_from_flow(depth = rep(1,5), velocity = rep(1,5), width=1)$Q, 1)
	
	# multiple measurememts
	expect_equal(q_from_flow(depth = rep(1,5), velocity = c(1,1,1,2,2), 
				flowID = c(1,1,1,2,2), width=c(1,1), widthID = c(1,2))$Q, c(1,2))
})
