library("WatershedTools")
library(raster)
library(sp)
gisBase <- getGISBase()

testDEM <- raster(system.file("testdata/testDEM.grd", package="WatershedTools"))
testDEM <- projectRaster(testDEM, crs=CRS("+init=epsg:3035"))
gs <- GrassSession(testDEM, layerName = "dem", gisBase = gisBase)
gs <- fillDEM("dem", filledDEM = "filledDEM", probs = "problems", gs = gs)
gs <- drainageAccumulation("filledDEM", accumulation = "accum", drainage = "drain",
	gs = gs)
gs <- extractStream(dem = "filledDEM", accumulation = "accum", qthresh = 0.998,
	outputName = "streamRas", gs = gs)
streamRas <- GSGetRaster("streamRas", gs)
drainage <- GSGetRaster("drain", gs)
accum <- GSGetRaster("accum", gs)
coords <- coordinates(accum)[which.max(values(accum)),, drop=FALSE]
streamCrop <- cropToCatchment(coords, streamRaster = streamRas, drainage = "drain", gs = gs)

vals <- values(streamCrop)
pts <- coordinates(streamCrop)[!is.na(vals),]
caFile <- "inst/testdata/testCA.rds"
if(file.exists(caFile)) {
	catchArea <- readRDS(caFile)
} else {
	catchArea <- catchment(pts, "drain", gs)
	catchArea <- data.frame(pts, A = catchArea)
	saveRDS(catchArea, caFile)
}


# discharge for northernmost point
q <- 5
A <- catchArea[which.max(catchArea$y),'A']
qMat <- discharge_scaling(catchArea, data.frame(A = A, Q = q))
geom <- hydraulic_geometry(qMat$A)
elevCrop <- crop(testDEM, streamCrop)
elevCrop <- resample(elevCrop, streamCrop, 'bilinear')
accumCrop <- crop(accum, streamCrop)
drainCrop <- crop(drainage, streamCrop)
caRaster <- rasterFromXYZ(catchArea, crs=proj4string(elevCrop))
coordinates(geom) <- pts
proj4string(geom) <- proj4string(elevCrop)
gridded(geom) <- TRUE
geom <- stack(geom)

testWS <- Watershed(streamCrop, drainCrop, elevCrop, accumCrop, caRaster, geom)
saveRDS(testWS, "inst/testdata/testWS.rds")
