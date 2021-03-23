library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)

setwd("E:/bathy")
bathy <- raster("E:/bathy/gebco_2020_all.tif") 
reefs <- st_read("E:/currentcoral/reefs.shp")


cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=bathy,funct=function(x){mean(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
  names(dfext) <- c("pointid","grid_code", "iscoral","ID","bathy","geometry")
  dfext$suit <- with(dfext, ifelse(bathy >= -40,1, 0)) #suitable depth is shallower than 40m (>40)
  #bathypt <- inner_join(reef,dfext,by="pointid")
  #bathydf <- as.data.frame(bathypt)
  write.csv(dfext,"bathy.csv", row.names=FALSE)
  #rasterize(bathypt,r008333,field="bathy",fun=max,na.rm=T,filename="bathy1.tif",overwrite=T)
  #rasterize(bathypt,r008333,field="suit",fun=max,na.rm=T,filename="bathy_suit1.tif",overwrite=T)

stopCluster(cl)



#get count of suitability
count(dfext, vars=suit)

  a <- read.csv(file="E:/bathy/bathy.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",5),rep("NULL",2),"numeric"))
  names(a) <- c("pointid","grid_code","iscoral","ID","bathy","suit")
  b <- a[!(a$iscoral==0),]
  count(b, vars="suit")




#make 0,360 for suit raster 
bathy1 <- crop(bathy, extent(-180, 0, -90,90))
bathy2 <- crop(bathy, extent(0, 180, -90,90))   
extent(bathy1) <- c(180,360,-90,90)
bathy360 <- merge(bathy1,bathy2)
writeRaster(bathy360,filename="E:/bathy/bathy360.tif")


