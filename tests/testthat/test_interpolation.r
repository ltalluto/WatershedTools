context("Interpolation")
library("WatershedTools")
library(data.table)
library(lubridate)

ws <- readRDS(system.file("testdata/testWS.rds", package="WatershedTools"))
wtDat <- fread(system.file("testdata/water_temp.txt", package="WatershedTools"))
wtDat[, timestamp := as_datetime(timestamp, tz="Europe/Tirane")]
selDay <- dmy("29-04-2018")
wtDat <- wtDat[date(timestamp) == selDay]
wtDat[, minutes := as.integer((timestamp - min(timestamp))/60)]
pixLookup <- cbind(unique(wtDat$deploymentID), c(30, 7, 145))
wtDat$pixID <- pixLookup[match(wtDat$deploymentID, pixLookup[,1]),2]
target_times <- seq(0, max(wtDat$minutes), 35)


test_that("Interpolation in space", {
	## WRONG - this must be for a SINGLE time; so times needs to be NULL/NA, and need to subset
	## wtDat to have exactly 3 observations
	expect_error(spaceOnly <- interpolate(ws, sites = wtDat$pixID, times =wtDat$minutes, 
		values = wtDat$value, sites_out = 'all', times_out = NA), regex = NA)
	expect_equal(dims(spaceOnly), c(nrow(ws$data), 1))

})

test_that("Interpolation in time", {
	expect_error(timeOnly <- interpolate(ws, sites = wtDat$pixID, times =wtDat$minutes, 
		values = wtDat$value, sites_out = NA, times_out = target_times), regex = NA)
	expect_equal(dims(timeOnly), c(length(unique(wtDat$pixID)), length(target_times)))
	timeMelted <- melt(timeOnly, varnames=c("pixID", "minutes"), value.name = 'value')
	commonPoints <- merge(timeMelted, wtDat, by = c('pixID', 'minutes'), all=FALSE)
	expect_equal(commonPoints$value.x, commonPoints$value.y)

})

test_that("Interpolation in space and time", {
	expect_error(spaceTime <- interpolate(ws, sites = wtDat$pixID, times =wtDat$minutes, 
		values = wtDat$value, sites_out = 'all', times_out = target_times), regex = NA)
	expect_equal(dims(spaceOnly), c(nrow(ws$data), length(target_times)))

})

