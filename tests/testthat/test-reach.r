ws = readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
