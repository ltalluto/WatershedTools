setwd("..")
library(devtools)
load_all()

ws <- readRDS("~/work/projects/catchment_delineations/vjosa/res/vjosaWatershedSpring2018.rds")
ws$data$reachID <- renumberReaches(ws$data$reachID)
ws$connectivity <- reachByReachConn(ws, self = FALSE)
x <- 18456
system.time(conn1 <- WatershedTools:::connectAll(ws, x))
ws$data$test1 <- 0
ws$data$test1[conn1$upstream] <- -1
ws$data$test1[conn1$downstream] <- 1
plot(ws, variable = 'test1')
quartz()


system.time(conn2 <- WatershedTools:::connectAll2(ws, x))
ws$data$test2 <- 0
ws$data$test2[conn2$upstream] <- -1
ws$data$test2[conn2$downstream] <- 1
plot(ws, variable = 'test2')
identical(conn1, conn2)
all(sort(conn1$downstream) == sort(conn2$downstream))

# TODO; write proper tests for connectAll
# TODO; reconcile this with existing connect() function
# TODO; reconcile this with  dmat; this might be faster and more flexible
# 		just need to grab all upstream and downstream pixels like this
#		then sum the lengths by subsetting the ws object
