library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)
library(beepr)

###empirical

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/empirical")

reefs <- st_read("E:/currentcoral/reefs.shp")
oa <- raster("E:/oa/empirical/OAempiricalfnew.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="")
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)

  dfi <- subset(dfext,select=c("pointid","grid_code","iscoral","ID","oa_OAempiricalfnew"))
  names(dfi) <- c("pointid","grid_code","iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  #oapt <- inner_join(reefs,dfi,by="ID")
  #oadf <- as.data.frame(oapt)
  write.csv(dfi,"emp_2005new_suit.csv", row.names=FALSE)

stopCluster(cl)
finish <- Sys.time()
finish-start
beep()


#get count of suitability
  b <- dfi[!(dfi$iscoral==0),]
  c2 <- count(b, vars=suit)

beep()






###historic model

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/hist")

reefs <- st_read("E:/currentcoral/reefs.shp")
oa <- brick("E:/oa/co2sys_input/historical/OAhist_modelmedianf.tif")
#histrout <- paste("hist_",1845+(1:16)*10,"_oa.tif",sep="")
#histsuitout <- paste("hist_",1845+(1:16)*10,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1845+(1:16)*10,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:16){
  x <- paste0("oa_OAhist_modelmedianf.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  #oapt <- inner_join(reef,dfi,by="pointid")
  #oadf <- as.data.frame(oapt)
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  #rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  #rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start
beep()



#get count of suitability
csvfiles <- list.files("E:/oa/hist",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars="suit")
  y <- add_column(c,year=paste0(i*10+1845),.before=T)
  count <- rbind(y,count)
} 





###rcp26

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/rcp26")

reefs <- st_read("E:/currentcoral/reefs.shp")
oa <- brick("E:/oa/co2sys_input/rcp26/OArcp26_modelmedianf.tif")
#rcp26rout <- paste("rcp26_",2000+(1:20)*5,"_oa.tif",sep="")
#rcp26suitout <- paste("rcp26_",2000+(1:20)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("oa_OArcp26_modelmedianf.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  #oapt <- inner_join(reef,dfi,by="pointid")
  #oadf <- as.data.frame(oapt)
  write.csv(dfi,rcp26suitcsvout[i], row.names=FALSE)
  #rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
  #rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start
beep()


#get count of suitability
csvfiles <- list.files("E:/oa/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars="suit")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 




###rcp45

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/rcp45")

reefs <- st_read("E:/currentcoral/reefs.shp")
oa <- brick("E:/oa/co2sys_input/rcp45/OArcp45_modelmedianf.tif")
#rcp45rout <- paste("rcp45_",2000+(1:20)*5,"_oa.tif",sep="")
#rcp45suitout <- paste("rcp45_",2000+(1:20)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("oa_OArcp45_modelmedianf.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  #oapt <- inner_join(reef,dfi,by="pointid")
  #oadf <- as.data.frame(oapt)
  write.csv(dfi,rcp45suitcsvout[i], row.names=FALSE)
  #rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  #rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/oa/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars="suit")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 





###rcp85

####
#extract reef suit vals
#find omega arag for each point
#classify as suitable or not - save as .csv
###

setwd("E:/oa/rcp85")

reefs <- st_read("E:/currentcoral/reefs.shp")
oa <- brick("E:/oa/co2sys_input/rcp85/OArcp85_modelmedianf.tif")
#rcp85rout <- paste("rcp85_",2000+(1:20)*5,"_oa.tif",sep="")
#rcp85suitout <- paste("rcp85_",2000+(1:20)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=oa,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("oa_OArcp85_modelmedianf.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","oa")
  dfi$suit <- with(dfi, ifelse(oa >= 3.3, 1, 0)) #if omega arag is > 3.3, it is suitable(1). otherwise not (0)
  #oapt <- inner_join(reef,dfi,by="pointid")
  #oadf <- as.data.frame(oapt)
  write.csv(dfi,rcp85suitcsvout[i], row.names=FALSE)
  #rasterize(oapt,r008333,field="oa",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  #rasterize(oapt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/oa/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars="suit")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 
