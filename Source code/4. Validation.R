library(rgdal)
library(raster)
paths = ""

setwd(paste(paths,"Input_Data",sep=""))
study_area = readOGR(".",layer="RetunPeriodStudyAreaProjected")

setwd(path)
a = raster("final_a.tif")
u = raster("final_u.tif")

a = mask(a, study_area)
u = mask(u, study_area)

model_depth_10 = u - a* log(-log(1-(1/10)))
model_depth_50 = u - a* log(-log(1-(1/50)))
model_depth_100 = u - a* log(-log(1-(1/100)))
model_depth_500 = u - a* log(-log(1-(1/500)))

depth_10 = crop(depth_10,extent(study_area))
model_depth_10 = crop(model_depth_10,extent(study_area))
depth_50 = crop(depth_50,extent(study_area))
model_depth_50 = crop(model_depth_50,extent(study_area))
depth_100 = crop(depth_100,extent(study_area))
model_depth_100 = crop(model_depth_100,extent(study_area))
depth_500 = crop(depth_500,extent(study_area))
model_depth_500 = crop(model_depth_500,extent(study_area))

masked_depth_10 = overlay(depth_10, model_depth_10, fun = function(x,y) { y[is.na(x)]= NA; y})
masked_depth_50 = overlay(depth_50, model_depth_50, fun = function(x,y) { y[is.na(x)]= NA; y})
masked_depth_100 = overlay(depth_100, model_depth_100, fun = function(x,y) { y[is.na(x)]= NA; y})
masked_depth_500 = overlay(depth_500, model_depth_500, fun = function(x,y) { y[is.na(x)]= NA; y})

difference_10 = (depth_10 - masked_depth_10)
difference_50 = (depth_50 - masked_depth_50)
difference_100 = (depth_100 - masked_depth_100)
difference_500 = (depth_500 - masked_depth_500)

delta_modeled = na.omit(as.data.frame(difference_10))
ds_depth10 = c(mean(delta_modeled$layer),sd(delta_modeled$layer),min(delta_modeled$layer),max(delta_modeled$layer),sqrt(mean(delta_modeled$layer^2)))
delta_modeled = na.omit(as.data.frame(difference_50))
ds_depth50 = c(mean(delta_modeled$layer),sd(delta_modeled$layer),min(delta_modeled$layer),max(delta_modeled$layer),sqrt(mean(delta_modeled$layer^2)))
delta_modeled = na.omit(as.data.frame(difference_100))
ds_depth100 = c(mean(delta_modeled$layer),sd(delta_modeled$layer),min(delta_modeled$layer),max(delta_modeled$layer),sqrt(mean(delta_modeled$layer^2)))
delta_modeled = na.omit(as.data.frame(difference_500))
ds_depth500 = c(mean(delta_modeled$layer),sd(delta_modeled$layer),min(delta_modeled$layer),max(delta_modeled$layer),sqrt(mean(delta_modeled$layer^2)))

Metrics = c("Average","Std_Dev","Min","Max","RMSE")
df_u1 = cbind(Metrics,ds_depth10,ds_depth50,ds_depth100,ds_depth500)
df_u1 = t(df_u1)
########################################
Metrics = c("Average","Std_Dev","Min","Max","RMSE")
df_u1 = cbind(Metrics,ds_depth500)

df_u1 = cbind(df_u1,ds_depth500)
colnames(df_u1) = c('Metrics','MA','OK','NN','IDW')

setwd(path)
write.csv(df_u1, "Validation_table_interpolations_500.csv")

