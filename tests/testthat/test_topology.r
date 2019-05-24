setwd("..")
library(devtools)
load_all()

ws <- readRDS("~/work/projects/metabolism/catchment_delineations/vjosa/res/vjosaWatershedSpring2018.rds")
x <- 18456

test <- accumulateOne(ws, x)
ws$data$test <- test
plot(ws, 'test', transform = function(x) log(abs(x)))


# TODO; write proper tests for connectAll
# TODO; reconcile this with existing connect() function
# TODO; reconcile this with  dmat; this might be faster and more flexible
# 		just need to grab all upstream and downstream pixels like this
#		then sum the lengths by subsetting the ws object
