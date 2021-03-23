library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)

setwd("E:/currentcoral")
coraloriginal <- raster("D:/coast/coral_rocky/coral_raster_Shift3.tif") 
coral <- reclassify(coraloriginal,c(-Inf,0,0,1,Inf,1),filename="currentcoral.tif")
reef <- st_read("D:/distribution/all_reef_pt.shp")
r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


# cores<- detectCores()-1
# cl <- makeCluster(cores, output="") 
# registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef,ras=coral,funct=function(x){mean(x,na.rm=T)},parallel=T)
st_write(reefext, dsn = "currentcoral", driver = "ESRI Shapefile")

dfext <- as.data.frame(reefext)
names(dfext) <- c("pointid","grid_code", "geometry","coral")
coralpt <- inner_join(reef,dfext,by="pointid")
coraldf <- as.data.frame(coralpt)
write.csv(coraldf,"coral.csv", row.names=FALSE)
rasterize(coralpt,r008333,field="coral",fun=max,na.rm=T,filename="coral.tif",overwrite=T)

# stopCluster(cl)



#get count of suitability
count(dfext, vars=coral)
