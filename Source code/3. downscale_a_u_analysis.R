library(rgdal)
library(raster)
paths = ""
#########################
setwd(paste(paths,"Preliminary a u analysis",sep=""))
u = raster("u1.tif")
depth_100 = raster("qc_depth_100.tif")
depth_500 = raster("qc_depth_500.tif")

#####################################################################
extents = extent(u)
depth_100 = crop(depth_100,extents)
depth_500 = crop(depth_500,extents)

#remove the cells of the first analysis to work only the remaining cells
au_depth_100 = overlay(u, depth_100, fun = function(x,y) { y[!is.na(x)]= NA; y})
au_depth_500 = overlay(u, depth_500, fun = function(x,y) { y[!is.na(x)]= NA; y})

au_depth_100[is.na(au_depth_100)] = threshold_value
au_depth_100[is.na(au_depth_500)] = NA
au_depth_500[is.na(au_depth_500)] = threshold_value

setwd(paste(path,'Extrapolation',sep=""))
au_depth_5k = raster("depth_5k.tif")
au_depth_10k = raster("depth_10k.tif")
au_depth_15k = raster("depth_15k.tif")
au_depth_20k = raster("depth_20k.tif")
#####################################################################
#a and u estimation
rp = c(log(-log(1-(1/100))),log(-log(1-(1/500))), log(-log(1-(1/5000))),log(-log(1-(1/10000))),log(-log(1-(1/15000))),log(-log(1-(1/20000))))
threshold_value = -0.05

df_depth_100 = as.data.frame(au_depth_100, xy=T, na.rm=F)
df_depth_500 = as.data.frame(au_depth_500, xy=F, na.rm=F)
df_depth_5k = as.data.frame(au_depth_5k, xy=F, na.rm=F)
df_depth_10k = as.data.frame(au_depth_10k, xy=F, na.rm=F)
df_depth_15k = as.data.frame(au_depth_15k, xy=F, na.rm=F)
df_depth_20k = as.data.frame(au_depth_20k, xy=F, na.rm=F)

df_all = cbind(df_depth_100,df_depth_500,df_depth_5k,df_depth_10k,df_depth_15k,df_depth_20k)
colnames(df_all) = c('x','y','depth_100','depth_500','depth_5k','depth_10k','depth_15k','depth_20k')
df_all1 = df_all[!is.na(df_all$depth_5k),]

fun2 = function(x,y){
  if(y[3]<0){
    if(y[2]>0){
      u_a = lm(y[1:2]~x[1:2])$coefficients
      depth_100year = u_a[1] + u_a[2] * log(-log(1-(1/100)))
      value = -0.05
      while (depth_100year>0){
        value = value+threshold_value
        y[1] = value
        u_a = lm(y[1:2]~x[1:2])$coefficients
        depth_100year = u_a[1] + u_a[2] * log(-log(1-(1/100)))
      }
    }
    if(y[2]<0){
      u_a = lm(y[3:6]~x[3:6])$coefficients
    }
  }
  if (y[3]>0){
    if(y[2]<0){
      u_a = lm(y[2:6]~x[2:6])$coefficients
      depth_500year = u_a[1] + u_a[2] * log(-log(1-(1/500)))
      value = threshold_value
      while (depth_500year>0){
        value = value+threshold_value
        y[2] = value
        u_a = lm(y[2:6]~x[2:6])$coefficients
        depth_500year = u_a[1] + u_a[2] * log(-log(1-(1/500)))
      }
    }
    if(y[2]>y[3]){
      u_a = lm(y[1:2]~x[1:2])$coefficients
      depth_100year = u_a[1] + u_a[2] * log(-log(1-(1/100)))
      value = threshold_value
      while (depth_100year>0){
        value = value+threshold_value
        y[1] = value
        u_a = lm(y[1:2]~x[1:2])$coefficients
        depth_100year = u_a[1] + u_a[2] * log(-log(1-(1/100)))
      }
    }
    if(y[2]<y[3]){
      u_a = lm(y~x)$coefficients
      depth_100year = u_a[1] + u_a[2] * log(-log(1-(1/100)))
      value = threshold_value
      while (depth_100year>0){
        value = value+threshold_value
        y[1] = value
        u_a = lm(y~x)$coefficients
        depth_100year = u_a[1] + u_a[2] * log(-log(1-(1/100)))
      }
    }
  }
  return(u_a)
}

reg = mapply(fun2,as.data.frame(rp), as.data.frame(t(df_all1[3:8])))
reg = as.data.frame(t(reg))

df_all1 = cbind(df_all1,reg[1],-1*reg[2])
colnames(df_all1) = c('x','y','depth_100','depth_500','depth_5k','depth_10k','depth_15k','depth_20k','u','a')

u_df = cbind(df_all1[1:2], df_all1[9])
a_df = cbind(df_all1[1:2], df_all1[10])
colnames(u_df) = c('x','y','u')
colnames(a_df) = c('x','y','a')

a = rasterFromXYZ (a_df, crs="EPSG:26916" )
u = rasterFromXYZ (u_df, crs="EPSG:26916" )
#####################################################################
dir.create(paste(path,"Downscale a u analysis",sep=""), showWarnings = FALSE)
setwd(paste(path,"Downscale a u analysis",sep=""))
writeRaster(a,'a2.tif',overwrite=TRUE)
writeRaster(u,'u2.tif',overwrite=TRUE)

setwd(path)
a1 = raster("a1.tif")
u1 = raster("u1.tif")

final_a = overlay(a1,a,fun=function(x,y){x[is.na(x)]=y[is.na(x)];x})
final_u = overlay(u1,u,fun=function(x,y){x[is.na(x)]=y[is.na(x)];x})
writeRaster(final_a,'final_a.tif')
writeRaster(final_u,'final_u.tif')
