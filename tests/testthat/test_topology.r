setwd("..")
library(devtools)
load_all()

library(ggplot2)

ws <- readRDS("~/work/projects/metabolism/catchment_delineations/vjosa/res/vjosaWatershedSpring2018.rds")
x <- c(1369, 22506, 13792, 8, 44100, 38600)
pl <- plot(ws)
for(i in x)
	pl <- pl + geom_point(x = ws[i,'x'], y=ws[i,'y'], colour = 'red')
pl


## HUH. 44100 produces a NULL when going upstream. WHY>? doiwns is fine
ws$test <- NA
ws$data$test[downs[[5]][,1]] <- downs[[5]][,2]
plot(ws, variable='test')


system.time(test1 <- accumulateOne(ws, x, variable = 'length')[y])
system.time(test2 <- {test2 <- accumulate(ws, x, y); test2[test2[,1] == y, 2]})
c(test1, test2)

system.time(test3 <- accumulateOne(ws, x, variable = 'length'))
system.time(test4 <- accumulate(ws, x, Inf))
c(test3[1], test4[test4[,1] == 1,2])

system.time(test5 <- accumulateOne(ws, y, variable = 'length'))
system.time(test6 <- accumulate(ws, Inf, y, parallel = TRUE))
c(test5[1,2], test6[1,2])

test[y]


ws$data$test <- test
plot(ws, 'test', transform = function(x) log(abs(x)))


# TODO; write proper tests for connectAll
# TODO; reconcile this with existing connect() function
# TODO; reconcile this with  dmat; this might be faster and more flexible
# 		just need to grab all upstream and downstream pixels like this
#		then sum the lengths by subsetting the ws object
