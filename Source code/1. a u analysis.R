library(rgdal)
library(raster)
paths = ""

###############################################################
#Import rasters and quality control
setwd(paste(paths,"Input_Data",sep=""))
depth_10 = raster("10YearFloodDepth.tif")
depth_50 = raster("50YearFloodDepth.tif")
depth_100 = raster("100YearFloodDepth.tif")
depth_500 = raster("500YearFloodDepth.tif")
dem = raster("DEM.tif")
study_area = readOGR(".",layer="RetunPeriodStudyAreaProjected")

###quality control
qc1 = depth_10>=depth_50
df= na.omit(as.data.frame(qc1))
if(length(df[df$layer == 'TRUE',])>0){
  depth_10[qc1==1] = NA
}
qc2 = depth_10>=depth_100
df= na.omit(as.data.frame(qc2))
if (length(df[df$layer == 'TRUE',])>0){
  depth_10[qc2==1] = NA
}
qc3 = depth_10>=depth_500
df= na.omit(as.data.frame(qc3))
if (length(df[df$layer == 'TRUE',])>0){
  depth_10[qc3==1] = NA
}
qc4 = depth_10>0 & is.na(depth_50)
df= na.omit(as.data.frame(qc4))
if (length(df[df$layer == 'TRUE',])>0){
  depth_10[qc4==1] = NA
}
qc5 = depth_10>0 & is.na(depth_100)
df= na.omit(as.data.frame(qc5))
if (length(df[df$layer == 'TRUE',])>0){
  depth_10[qc5==1] = NA
}
qc6 = depth_10>0 & is.na(depth_500)
df= na.omit(as.data.frame(qc6))
if (length(df[df$layer == 'TRUE',])>0){
  depth_10[qc6==1] = NA
}
qc7 = depth_10<=0
df= na.omit(as.data.frame(qc7))
if (length(df[df$layer == 'TRUE',])>0){
  depth_10[qc7==1] = NA
}
qc8 = depth_50>=depth_100
df= na.omit(as.data.frame(qc8))
if (length(df[df$layer == 'TRUE',])>0){
  depth_50[qc8==1] = NA
}
qc9 = depth_50>=depth_500
df= na.omit(as.data.frame(qc9))
if (length(df[df$layer == 'TRUE',])>0){
  depth_50[qc9==1] = NA
}
qc10 = depth_50>0 & is.na(depth_100)
df= na.omit(as.data.frame(qc10))
if (length(df[df$layer == 'TRUE',])>0){
  depth_50[qc10==1] = NA
}
qc11 = depth_50>0 & is.na(depth_500)
df= na.omit(as.data.frame(qc11))
if (length(df[df$layer == 'TRUE',])>0){
  depth_50[qc11==1] = NA
}
qc12 = depth_50<=0
df= na.omit(as.data.frame(qc12))
if (length(df[df$layer == 'TRUE',])>0){
  depth_50[qc12==1] = NA
}
qc13 = depth_100>=depth_500
df= na.omit(as.data.frame(qc13))
if (length(df[df$layer == 'TRUE',])>0){
  depth_100[qc13==1] = NA
}
qc14 = depth_100>0 & is.na(depth_500)
df= na.omit(as.data.frame(qc14))
if (length(df[df$layer == 'TRUE',])>0){
  depth_100[qc14==1] = NA
}
qc15 = depth_100<=0
df= na.omit(as.data.frame(qc15))
if (length(df[df$layer == 'TRUE',])>0){
  depth_100[qc15==1] = NA
}
qc16 = depth_500<=0
df= na.omit(as.data.frame(qc16))
if (length(df[df$layer == 'TRUE',])>0){
  depth_500[qc16==1] = NA
}
qc_list = list(qc1,qc2,qc3,qc4,qc5,qc6,qc7,qc8,qc9,qc10,qc11,qc12,qc13,qc14,qc15,qc16)
qc_list$fun = sum
qc_list$na.rm = TRUE
qc_all = do.call(mosaic,qc_list)

dir.create(paste(paths,"Preliminary a u analysis",sep=""), showWarnings = FALSE)
setwd(paste(paths,"Preliminary a u analysis",sep=""))
writeRaster(depth_10,'qc_depth_10.tif',overwrite=T)
writeRaster(depth_50,'qc_depth_50.tif',overwrite=T)
writeRaster(depth_100,'qc_depth_100.tif',overwrite=T)
writeRaster(depth_500,'qc_depth_500.tif',overwrite=T)
###############################################################
#Gumbel fit
threshold_value = -0.05
au_depth_10 = depth_10
au_depth_50 = depth_50
au_depth_100 = depth_100
au_depth_500 = depth_500

au_depth_10[is.na(au_depth_10)] = threshold_value
au_depth_10[is.na(au_depth_50)] = NA
au_depth_50[is.na(au_depth_50)] = threshold_value
au_depth_50[is.na(au_depth_100)] = NA

rp = c(log(-log(1-(1/2))),log(-log(1-(1/10))),log(-log(1-(1/50))), log(-log(1-(1/100))),log(-log(1-(1/500))))
df_depth_10 = as.data.frame(au_depth_10, xy=T, na.rm=F)
df_depth_50 = as.data.frame(au_depth_50, xy=F, na.rm=F)
df_depth_100 = as.data.frame(au_depth_100, xy=F, na.rm=F)
df_depth_500 = as.data.frame(au_depth_500, xy=F, na.rm=F)
df_all = cbind(df_depth_10,df_depth_50,df_depth_100,df_depth_500)
df_all1 = df_all[!is.na(df_all$depth_500),]
df_all1$depth_2 = ifelse(df_all1$depth_10>0, -0.05, NA)
df_all1 = df_all1[,c(1,2,7,3,4,5,6)]

fun1 = function(x,y){
  points= 0 
  for (i in y){
   if (!(is.na(i))){
     points = points+1
   } 
  }
  if (points==5){
    u_a = lm(y[2:5]~x[2:5])$coefficients
    if (u_a[1]>0){
      u_a = lm(y~x)$coefficients
      value = threshold_value
      while (u_a[1]>0){
        value = value+threshold_value
        y[1] = value 
        u_a = lm(y~x)$coefficients
      }
    }
  }
  if (points==4){
    u_a = lm(y[3:5]~x[3:5])$coefficients
    depth_10year = u_a[1] + u_a[2] * log(-log(1-(1/10)))
    if(depth_10year>0){
      u_a = lm(y[2:5]~x[2:5])$coefficients
      value = threshold_value
      depth_10year = u_a[1] + u_a[2] * log(-log(1-(1/10)))
      while (depth_10year>0){
        value = value+threshold_value
        y[2] = value
        u_a = lm(y[2:5]~x[2:5])$coefficients
        depth_10year = u_a[1] + u_a[2] * log(-log(1-(1/10)))
      }
    }
  }
  if (points==3){
    u_a = lm(y[4:5]~x[4:5])$coefficients
    depth_50year = u_a[1] + u_a[2] * log(-log(1-(1/50)))
    if(depth_50year>0){
      u_a = lm(y[3:5]~x[3:5])$coefficients
      depth_50year = u_a[1] + u_a[2] * log(-log(1-(1/50)))
      value = threshold_value
      while (depth_50year>0){
        value = value+threshold_value
        y[3] = value
        u_a = lm(y[3:5]~x[3:5])$coefficients
        depth_50year = u_a[1] + u_a[2] * log(-log(1-(1/50)))
      }
    }
  }
  if (points==2||points==1){
    u_a = c(NA,NA)
  }
  return(u_a)
}

reg = mapply(fun1,as.data.frame(rp), as.data.frame(t(df_all1[3:7])))
reg = as.data.frame(t(reg))

u_df = cbind(df_all1[1:2], reg[1])
a_df = cbind(df_all1[1:2], -1*reg[2])
colnames(u_df) = c('x','y','u')
colnames(a_df) = c('x','y','a')

a = rasterFromXYZ (a_df, crs="EPSG:26916" )
u = rasterFromXYZ (u_df, crs="EPSG:26916" )

extents = extent(a)
dem1 = crop(dem,extents)
############################################################################
writeRaster(a,'a1.tif',overwrite=T)
writeRaster(u,'u1.tif',overwrite=T)
writeRaster(u,'dem1.tif',overwrite=T)
