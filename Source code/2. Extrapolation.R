library(rgdal)
library(raster)
paths = ""
###############################################################
setwd(paste(paths,"Input_ Data",sep=""))
study_area = readOGR(".",layer="RetunPeriodStudyAreaProjected")

setwd(paste(paths,"Preliminary a u analysis",sep=""))
u = raster("u1.tif")
a = raster("a1.tif")
dem = raster("dem1.tif")

depth_5k  =  u - a* log(-log(1-(1/5000)))
depth_10k  =  u - a* log(-log(1-(1/10000))) 
depth_15k  =  u - a* log(-log(1-(1/15000))) 
depth_20k  =  u - a* log(-log(1-(1/20000))) 

elev_5k = depth_5k + dem
elev_10k = depth_10k + dem
elev_15k = depth_15k + dem
elev_20k = depth_20k + dem
########Moving average##########################
window_size = 31

ma_elev_5k = focal(elev_5k, w=matrix(1,window_size,window_size), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA, NAonly=TRUE)
ma_elev_10k = focal(elev_10k, w=matrix(1,window_size,window_size), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA, NAonly=TRUE)
ma_elev_15k = focal(elev_15k, w=matrix(1,window_size,window_size), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA, NAonly=TRUE)
ma_elev_20k = focal(elev_20k, w=matrix(1,window_size,window_size), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA, NAonly=TRUE)

smooth_elev_5k = focal(ma_elev_5k,w=matrix(1,3,3), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA)
smooth_elev_10k = focal(ma_elev_10k,w=matrix(1,3,3), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA)
smooth_elev_15k = focal(ma_elev_15k,w=matrix(1,3,3), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA)
smooth_elev_20k = focal(ma_elev_20k,w=matrix(1,3,3), fun=mean, na.rm =TRUE, pad= TRUE, padValue = NA)

###
masked_elev_5k = overlay(elev_5k, smooth_elev_5k, fun = function(x,y) { y[!is.na(x)]= NA; y})
masked_depth_5k = overlay(masked_elev_5k,dem,fun = function(x,y) { is.na(y)= NA; x-y} )

masked_elev_10k = overlay(elev_10k, smooth_elev_10k, fun = function(x,y) { y[!is.na(x)]= NA; y})
masked_depth_10k = overlay(masked_elev_10k,dem,fun = function(x,y) { is.na(y)= NA; x-y} )

masked_elev_15k = overlay(elev_15k, smooth_elev_15k, fun = function(x,y) { y[!is.na(x)]= NA; y})
masked_depth_15k = overlay(masked_elev_15k,dem,fun = function(x,y) { is.na(y)= NA; x-y} )

masked_elev_20k = overlay(elev_20k, smooth_elev_20k, fun = function(x,y) { y[!is.na(x)]= NA; y})
masked_depth_20k = overlay(masked_elev_20k,dem,fun = function(x,y) { is.na(y)= NA; x-y} )

##########################
dir.create(paste(path,"Extrapolation",sep=""), showWarnings = FALSE)
setwd(paste(path,"Extrapolation",sep=""))

writeRaster(masked_depth_5k, "depth_5k.tif",overwrite=TRUE)
writeRaster(masked_depth_10k, "depth_10k.tif",overwrite=TRUE)
writeRaster(masked_depth_15k, "depth_15k.tif",overwrite=TRUE)
writeRaster(masked_depth_20k, "depth_20k.tif",overwrite=TRUE)

