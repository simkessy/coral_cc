library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)



###
#historic model - hyde 3.1
###

setwd("E:/pop_hr/hist_hyde")

popd1850 <- raster("popd_1850AD.asc")


#calculate 5 year averages
histhydefiles <- list.files("E:/pop_hr/hist_hyde", pattern="popd_", full.names=T)
histhyde5yrnames <- paste("histhyde_",1845+10*(1:15),"_dens.tif",sep="")

for (i in 1:16){
  x <- raster(histhydefiles[i])
  n <- i+1
  y <- raster(histhydefiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=histhyde5yrnames[i])
}




####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("E:/currentcoral/reefs.shp")
IDmatch <- read.csv("E:/land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
hydefiles <- list.files("E:/pop_hr/hist_hyde",pattern="5_dens_bc", full.names=T)
nhydefiles <- length(hydefiles)
hydesuitcsvout <- paste("histhyde_",1845+10*(1:16),"_suit_bc.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nhydefiles){
  pop <- raster(hydefiles[i])
  crs(pop) <- crs(reefs)
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,hydesuitcsvout[i], row.names=FALSE)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("E:/pop_hr/hist_hyde",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(1845+10*i),.before=T)
  count <- rbind(y,count)
} 




###
#historic empirical
###

setwd("E:/pop_hr/hist_gpw")

histgpw <- raster("E:/pop_hr/hist_gpw/gpw_v4_population_density_rev11_2005_30_sec.tif")
#histgpw3 <- raster("E:/pop_hr/hist_gpw/gpw_v3.tif")
#histgpw1 <- raster("E:/pop_hr/hist_gpw/popdynamics-global-pop-density-time-series-estimates-2000-geotiff/popdynamics-global-pop-density-time-series-estimates_2000.tif")

####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("E:/currentcoral/reefs.shp")
IDmatch <- read.csv("E:/land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

  pop <- histgpw
  crs(pop) <- crs(reefs)
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,"histgpw_v1_2005_suit.csv", row.names=FALSE)
  gc()

stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
  suitcor <- d[!(d$iscoral==0),]
  count(suitcor, vars="suit")

  #make shapefile
  #suitpt <- merge(reefs,d,by="ID")
  #st_write(suitpt, dsn = "emp_popgpw3_2005_suit", driver = "ESRI Shapefile")



  
  
#make coarse empirical dataset for comparison against historic model
hyde2005 <- raster("E:/pop_hr/hist_hyde/popd_2005_dens.asc")
gpw4v2005 <- raster("E:/pop_hr/hist_gpw/gpw_v4_population_density_rev11_2005_30_sec.tif")
gpw4v2005a <- aggregate(gpw4v2005,fact=10,fun=mean,na.rm=T,filename="E:/pop_hr/hist_gpw/gpw_v4_population_density_083.tif")





###
#SSP1
###

setwd("E:/pop_hr/ssp1")

#make area raster
r1.2010 <- raster("D:/human_population/population_projection/ssp1/SSP1_1km/ssp1_total_2010.nc4")
area <- area(reefbuff, filename="area_1km.tif", na.rm=FALSE) #raster layer where Cell values represent the size of the cell in km2

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp1 <- "D:/human_population/population_projection/ssp1/SSP1_1km/"
ssp1files <- list.files(path.ssp1, pattern="_total_", full.names=T)
nfiles <- length(ssp1files)
ssp1_list  <- paste("ssp1.20",1:nfiles,"0",sep="")
ssp1_dens_list <- paste("ssp1.20",1:nfiles,"0d",sep="")
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp1densfiles <- paste("ssp1_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp1files[i]) #open each raster
  assign(ssp1_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp1densfiles[i]) #calculate density
  assign(ssp1_dens_list[i],d)
}

#calculate 5 year averages
ssp1dfiles <- list.files("E:/pop_hr/ssp1", pattern="_dens.tif", full.names=T)
ssp15yrnames <- paste("ssp1_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp1dfiles[i])
  n <- i+1
  y <- raster(ssp1dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp15yrnames[i])
}


#calculate 2005 pop density using the 2000 base year
#calculate 2000 pop density, then calculate average with 2010
base2000 <- raster("E:/pop_hr/baseYr_total_2000.tif")
area <- raster("E:/pop_hr/area_1km.tif")
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
y <- extend(base2000,extent(reefbuff),value=NA) #change extent so matches reef buffer 
extent(y) <- extent(reefbuff)
z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
dens2000 <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename="E:/pop_hr/base_2000_dens.tif") #calculate density
#calculate 2005 density
ssp1dens2010 <- raster("ssp1_2010_dens.tif")
ssp1dens2005 <- overlay(ssp1dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp1_2005_dens.tif") #calculate 2005 pop dens



####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("E:/currentcoral/reefs.shp")
IDmatch <- read.csv("E:/land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp1files <- list.files("E:/pop_hr/ssp1",pattern="_dens", full.names=T)
nssp1files <- length(ssp1files)
#ssp1rout <- paste("ssp1_",2000+(1:nssp1files+1)*5,"_logpop.tif",sep="")
#ssp1suitout <- paste("ssp1_",2000+(1:nssp1files+1)*5,"_suit.tif",sep="")
ssp1suitcsvout <- paste("ssp1_",2000+(1:nssp1files)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nssp1files){
  pop <- raster(ssp1files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp1suitcsvout[i], row.names=FALSE)
  #rasterize(poppt,pop,field="logpop",fun=max,na.rm=T,filename=ssp1rout[i],overwrite=T)
  #rasterize(poppt,pop,field="suit",fun=max,na.rm=T,filename=ssp1suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("E:/pop_hr/ssp1",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 
 


#####
#####
####
#SSP2
####

setwd("E:/pop_hr/ssp2")

area <- raster("E:/pop_hr/area_1km.tif")

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp2 <- "D:/human_population/population_projection/ssp2/SSP2_1km/"
ssp2files <- list.files(path.ssp2, pattern="_total_", full.names=T)
nfiles <- length(ssp2files)
ssp2_list  <- paste("ssp2.20",1:nfiles,"0",sep="")
ssp2_dens_list <- paste("ssp2.20",1:nfiles,"0d",sep="")
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp2densfiles <- paste("ssp2_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp2files[i]) #open each raster
  assign(ssp2_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp2densfiles[i]) #calculate density
  assign(ssp2_dens_list[i],d)
}

#calculate 5 year averages
ssp2dfiles <- list.files("E:/pop_hr/ssp2", pattern="_dens.tif", full.names=T)
ssp25yrnames <- paste("ssp2_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp2dfiles[i])
  n <- i+1
  y <- raster(ssp2dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp25yrnames[i])
}

#calculate 2005 density
dens2000 <- raster("E:/pop_hr/base_2000_dens.tif")
ssp2dens2010 <- raster("ssp2_2010_dens.tif")
ssp2dens2005 <- overlay(ssp2dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp2_2005_dens.tif") #calculate 2005 pop dens


####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("E:/currentcoral/reefs.shp")
IDmatch <- read.csv("E:/land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp2files <- list.files("E:/pop_hr/ssp2",pattern="_dens", full.names=T)
nssp2files <- length(ssp2files)
#ssp2rout <- paste("ssp2_",2000+(1:nssp2files+1)*5,"_logpop.tif",sep="")
#ssp2suitout <- paste("ssp2_",2000+(1:nssp2files+1)*5,"_suit.tif",sep="")
ssp2suitcsvout <- paste("ssp2_",2000+(1:nssp2files)*5,"_suit.csv",sep="")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nssp2files){
  pop <- raster(ssp2files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp2suitcsvout[i], row.names=FALSE)
  #rasterize(poppt,pop,field="logpop",fun=max,na.rm=T,filename=ssp2rout[i],overwrite=T)
  #rasterize(poppt,pop,field="suit",fun=max,na.rm=T,filename=ssp2suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("E:/pop_hr/ssp2",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 






#####
#####
####
#SSP3
####

setwd("E:/pop_hr/ssp3")

area <- raster("E:/pop_hr/area_1km.tif")

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp3 <- "D:/human_population/population_projection/ssp3/SSP3_1km/"
ssp3files <- list.files(path.ssp3, pattern="_total_", full.names=T)
nfiles <- length(ssp3files)
ssp3_list  <- paste("ssp3.20",1:nfiles,"0",sep="")
ssp3_dens_list <- paste("ssp3.20",1:nfiles,"0d",sep="")
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp3densfiles <- paste("ssp3_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp3files[i]) #open each raster

  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp3densfiles[i]) #calculate density
  assign(ssp3_dens_list[i],d)
}

#calculate 5 year averages
ssp3dfiles <- list.files("E:/pop_hr/ssp3", pattern="_dens.tif", full.names=T)
ssp35yrnames <- paste("ssp3_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:9){
  x <- raster(ssp3dfiles[i])
  n <- i+1
  y <- raster(ssp3dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp35yrnames[i])
}

#calculate 2005 density
ssp3dens2010 <- raster("ssp3_2010_dens.tif")
ssp3dens2005 <- overlay(ssp3dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp3_2005_dens.tif") #calculate 2005 pop dens



####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("E:/currentcoral/reefs.shp")
IDmatch <- read.csv("E:/land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp3files <- list.files("E:/pop_hr/ssp3",pattern="_dens", full.names=T)
nssp3files <- length(ssp3files)
#ssp3rout <- paste("ssp3_",2000+(1:nssp3files+1)*5,"_logpop.tif",sep="")
#ssp3suitout <- paste("ssp3_",2000+(1:nssp3files+1)*5,"_suit.tif",sep="")
ssp3suitcsvout <- paste("ssp3_",2000+(1:nssp3files)*5,"_suit.csv",sep="")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:nssp3files){
  pop <- raster(ssp3files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp3suitcsvout[i], row.names=FALSE)
  #rasterize(poppt,pop,field="logpop",fun=max,na.rm=T,filename=ssp3rout[i],overwrite=T)
  #rasterize(poppt,pop,field="suit",fun=max,na.rm=T,filename=ssp3suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("E:/pop_hr/ssp3",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars="suit")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 





#####
#####
####
#SSP5
####

setwd("E:/pop_hr/ssp5")

area <- raster("E:/pop_hr/area_1km.tif")

#crop to only the coastal areas near reef sites. then calculate pop density for all years
path.ssp5 <- "D:/human_population/population_projection/ssp5/SSP5_1km/"
ssp5files <- list.files(path.ssp5, pattern="_total_", full.names=T)
nfiles <- length(ssp5files)
ssp5_list  <- paste("ssp5.20",1:nfiles,"0",sep="")
ssp5_dens_list <- paste("ssp5.20",1:nfiles,"0d",sep="")
reefbuff <- raster("E:/pop_hr/reeflandbuffer.tif")
ssp5densfiles <- paste("ssp5_20",1:nfiles,"0_dens.tif",sep="")

for (i in 1:nfiles){
  
  x <- raster(ssp5files[i]) #open each raster
  assign(ssp5_list[i],x)
  
  y <- extend(x,extent(reefbuff),value=NA) #change extent so matches reef buffer 
  extent(y) <- extent(reefbuff)
  z <- overlay(y,reefbuff,fun=function(x,y){(x*y)}) #get rid of inland areas we aren't interested in
  
  d <- overlay(z,area,fun=function(x,y,na.rm=T){return(x/y)},filename=ssp5densfiles[i]) #calculate density
  assign(ssp5_dens_list[i],d)
}

#calculate 5 year averages
ssp5dfiles <- list.files("E:/pop_hr/ssp5", pattern="_dens.tif", full.names=T)
ssp55yrnames <- paste("ssp5_20",1:nfiles,"5_dens.tif",sep="")

for (i in 1:nfiles){
  x <- raster(ssp5dfiles[i])
  n <- i+1
  y <- raster(ssp5dfiles[n])
  xy5 <- overlay(x,y,fun=function(x){ mean(x,na.rm=T)},filename=ssp55yrnames[i])
}


#calculate 2005 density
ssp5dens2010 <- raster("ssp5_2010_dens.tif")
ssp5dens2005 <- overlay(ssp5dens2010,dens2000,fun=function(x){ mean(x,na.rm=T)},filename="ssp5_2005_dens.tif") #calculate 2005 pop dens



####
#extract reef suit vals
#find maximum population density within 50km radius of reef site - save max pop as raster
#classify as suitable or not - save as .csv
###

polyfiles <- list.files("E:/reef_polygon",pattern="\\.shp$",full.names=T)
nfiles <- length(polyfiles)
reefs <- st_read("E:/currentcoral/reefs.shp")
IDmatch <- read.csv("E:/land_hurtt/polygonIDmatch.csv",header=T)
names(IDmatch) <- c("ObjectID","iscoral","ID","pointid")
ssp5files <- list.files("E:/pop_hr/ssp5",pattern="_dens", full.names=T)
nssp5files <- length(ssp5files)
#ssp5rout <- paste("ssp5_",2000+(1:nssp5files+1)*5,"_logpop.tif",sep="")
#ssp5suitout <- paste("ssp5_",2000+(1:nssp5files+1)*5,"_suit.tif",sep="")
ssp5suitcsvout <- paste("ssp5_",2000+(1:nssp5files)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 2:17){
  pop <- raster(ssp5files[i])
  df.extract1 <- foreach::foreach(i=1:5, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract2 <- foreach::foreach(i=6:10, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract3 <- foreach::foreach(i=11:15, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract4 <- foreach::foreach(i=16:20, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract5 <- foreach::foreach(i=20:25, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  df.extract6 <- foreach::foreach(i=26:30, .packages=c("sf","rgis","velox"), .combine=rbind) %dopar% {
    poly <- st_read(polyfiles[i])
    reefext <- rgis::fast_extract(sf=poly,ras=pop,funct=function(x){max(x,na.rm=T)},parallel=T)
  }
  gc()
  dff <- rbind(df.extract1,df.extract2)
  dff1 <- rbind(dff,df.extract3)
  dff2 <- rbind(dff1,df.extract4)
  dff3 <- rbind(dff2,df.extract5)
  dff4 <- rbind(dff3,df.extract6)
  dfext <- as.data.frame(dff4)
  dfext <- subset(dfext,select=-c(grid_code, BUFF_DIST, ORIG_FID,Shape_Leng,Shape_Area,geometry))
  names(dfext) <- c("pointid","zone","maxpop")
  result <- transform(dfext,pointid=as.numeric(pointid),maxpop=as.numeric(maxpop))
  result$logpop <- with(result, ifelse(maxpop==0 | is.na(maxpop),0,log10(maxpop)))
  result$suit <- with(result, ifelse(logpop <2, 1, 0))
  result$suit <- with(result, ifelse(is.na(logpop),1,suit))
  b<-merge(IDmatch,result,by="pointid",all.x=T,all.y=F)
  c <- b[!duplicated(b$ID),]
  d <- c[base::order(c$ID),]
  write.csv(d,ssp5suitcsvout[i], row.names=FALSE)
  #rasterize(poppt,pop,field="logpop",fun=max,na.rm=T,filename=ssp5rout[i],overwrite=T)
  #rasterize(poppt,pop,field="suit",fun=max,na.rm=T,filename=ssp5suitout[i],overwrite=T)
  gc()
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/pop_hr/ssp5",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  b <- a[!(a$iscoral==0),]
  c <- count(b, vars=suit)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}  




