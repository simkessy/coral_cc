library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)
library(rgis)
library(dplyr)
library(tibble)
library(plyr)


#####
#####
###RCP85-ssp3
#####
#####

setwd("F:/suit/rcp85-ssp3")

#calculate overall suit

#load coral points shapefile & sample raster
reefs <- st_read("E:/currentcoral/reefs.shp")
#r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 1's
#recalculate individual variable suits given filter


#load sst
sstfiles <- list.files("E:/sst/rcp85",pattern="_suit.csv$",full.names=T)
#load dhw
dhwfiles <- list.files("E:/dhw/rcp85",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hurtt/rcp85",pattern="_suit.csv$",full.names=T)
#load oa
oafiles <- list.files("E:/oa/rcp85",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp3",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp85",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")
pointsout <- paste("rcp85_",2000+(1:20)*5,"_suit",sep="")
#rasterccout <- paste("rcp85_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
#rasternewout <- paste("rcp85_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  #sst
  sst <- read.csv(file=sstfiles[i],row.names=NULL,header=T,
                  colClasses=c("numeric","NULL","numeric","numeric","NULL","numeric"))
  names(sst) <- c("pointid","iscoral","ID","sst_suit")
  sst[is.na(sst)] <- 1
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(dhw) <- c("ID","dhw_suit")
  dhw[is.na(dhw)] <- 1
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=T,
                   colClasses=c(rep("NULL",3),"numeric",rep("NULL",2),"numeric"))
  names(land) <- c("ID","land_suit")
  land[is.na(land)] <- 1
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=T,
                 colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(oa) <- c("ID","oa_suit")
  oa[is.na(oa)] <- 1
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric"))
  names(pop) <- c("ID","pop_suit")
  pop[is.na(pop)] <- 1
  #storms
    storms <- read.csv(file=stormfiles[i],row.names=NULL,header=T,
                       colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
    names(storms) <- c("ID","storm_suit")
    storms[is.na(storms)] <- 1
    #join together overall suit
    suitdf <- join_all(list(sst,dhw,land,oa,pop,storms),by="ID")
    suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
    suitdf$suitall <- suitdf$sst_suit *suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit

  #save csv & shapefile
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reefs,suitdf,by="ID")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("F:/suit/rcp85-ssp3",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2000),.before=T)
  countnew <- rbind(z,countnew)
} 






######
######
#RCP85-ssp5
######
######


setwd("F:/suit/rcp85-ssp5")

#calculate overall suit

#load coral points shapefile & sample raster
reefs <- st_read("E:/currentcoral/reefs.shp")
#r008333 <- raster("E:/bathy/gebco_2020_all.tif")



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 1's
#recalculate individual variable suits given filter


#load sst
sstfiles <- list.files("E:/sst/rcp85",pattern="_suit.csv$",full.names=T)
#load dhw
dhwfiles <- list.files("E:/dhw/rcp85",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hurtt/rcp85",pattern="_suit.csv$",full.names=T)
#load oa
oafiles <- list.files("E:/oa/rcp85",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp5",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp85",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp85_ssp5",2000+(1:20)*5,"_suit.csv",sep="")
pointsout <- paste("rcp85_ssp5",2000+(1:20)*5,"_suit",sep="")
#rasterccout <- paste("rcp85_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
#rasternewout <- paste("rcp85_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save
#rasterize. rotate so Pacific Ocean is centered
#print raster 
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  #sst
  sst <- read.csv(file=sstfiles[i],row.names=NULL,header=T,
                  colClasses=c("numeric","NULL","numeric","numeric","NULL","numeric"))
  names(sst) <- c("pointid","iscoral","ID","sst_suit")
  sst[is.na(sst)] <- 1
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(dhw) <- c("ID","dhw_suit")
  dhw[is.na(dhw)] <- 1
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=T,
                   colClasses=c(rep("NULL",3),"numeric",rep("NULL",2),"numeric"))
  names(land) <- c("ID","land_suit")
  land[is.na(land)] <- 1
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=T,
                 colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(oa) <- c("ID","oa_suit")
  oa[is.na(oa)] <- 1
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric"))
  names(pop) <- c("ID","pop_suit")
  pop[is.na(pop)] <- 1
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=T,
                     colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(storms) <- c("ID","storm_suit")
  storms[is.na(storms)] <- 1
  #join together overall suit
  suitdf <- join_all(list(sst,dhw,land,oa,pop,storms),by="ID")
  suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$sst_suit *suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit
  
  #save csv & shapefile
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reefs,suitdf,by="ID")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("F:/suit/rcp85-ssp5",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2000),.before=T)
  countnew <- rbind(z,countnew)
} 





#####
#####
###RCP45-ssp2
#####
#####

setwd("F:/suit/rcp45-ssp2")

#calculate overall suit



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter



#load sst
sstfiles <- list.files("E:/sst/rcp45",pattern="_suit.csv$",full.names=T)
#load dhw
dhwfiles <- list.files("E:/dhw/rcp45",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hurtt/rcp45",pattern="_suit.csv$",full.names=T)
#load oa
oafiles <- list.files("E:/oa/rcp45",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp2",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp45",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")
pointsout <- paste("rcp45_",2000+(1:20)*5,"_suit",sep="")
#rasterccout <- paste("rcp85_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
#rasternewout <- paste("rcp85_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  #sst
  sst <- read.csv(file=sstfiles[i],row.names=NULL,header=T,
                  colClasses=c("numeric","NULL","numeric","numeric","NULL","numeric"))
  names(sst) <- c("pointid","iscoral","ID","sst_suit")
  sst[is.na(sst)] <- 1
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(dhw) <- c("ID","dhw_suit")
  dhw[is.na(dhw)] <- 1
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=T,
                   colClasses=c(rep("NULL",3),"numeric",rep("NULL",2),"numeric"))
  names(land) <- c("ID","land_suit")
  land[is.na(land)] <- 1
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=T,
                 colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(oa) <- c("ID","oa_suit")
  oa[is.na(oa)] <- 1
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric"))
  names(pop) <- c("ID","pop_suit")
  pop[is.na(pop)] <- 1
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=T,
                     colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(storms) <- c("ID","storm_suit")
  storms[is.na(storms)] <- 1
  #join together overall suit
  suitdf <- join_all(list(sst,dhw,land,oa,pop,storms),by="ID")
  suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$sst_suit *suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit
  
  #save csv & shapefile
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reefs,suitdf,by="ID")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  
}

stopCluster(cl)
finish <- Sys.time()
finish-start



#get count of suitability
csvfiles <- list.files("F:/suit/rcp45-ssp2",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2000),.before=T)
  countnew <- rbind(z,countnew)
} 






#####
#####
###RCP26-ssp1
#####
#####

setwd("F:/suit/rcp26-ssp1")


#calculate overall suit



#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter



#load sst
sstfiles <- list.files("E:/sst/rcp26",pattern="_suit.csv$",full.names=T)
#load dhw
dhwfiles <- list.files("E:/dhw/rcp26",pattern="_suit.csv$",full.names=T)
#load land
landfiles <- list.files("E:/land_hurtt/rcp26",pattern="_suit.csv$",full.names=T)
#load oa
oafiles <- list.files("E:/oa/rcp26",pattern="_suit.csv$",full.names=T)
#load pop
popfiles <- list.files("E:/pop_hr/ssp1",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/rcp26",pattern="_suit.csv$",full.names=T)

#name out files
csvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")
pointsout <- paste("rcp26_",2000+(1:20)*5,"_suit",sep="")
#rasterccout <- paste("rcp85_",2000+(1:19+1)*5,"_coralsuit.tif",sep="")
#rasternewout <- paste("rcp85_",2000+(1:19+1)*5,"_newsuit.tif",sep="")


#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:20){
  #sst
  sst <- read.csv(file=sstfiles[i],row.names=NULL,header=T,
                  colClasses=c("numeric","NULL","numeric","numeric",rep("NULL",3),"numeric"))
  names(sst) <- c("pointid","iscoral","ID","sst_suit")
  sst[is.na(sst)] <- 1
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(dhw) <- c("ID","dhw_suit")
  dhw[is.na(dhw)] <- 1
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=T,
                   colClasses=c(rep("NULL",3),"numeric",rep("NULL",2),"numeric"))
  names(land) <- c("ID","land_suit")
  land[is.na(land)] <- 1
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=T,
                 colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(oa) <- c("ID","oa_suit")
  oa[is.na(oa)] <- 1
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric"))
  names(pop) <- c("ID","pop_suit")
  pop[is.na(pop)] <- 1
  #storms
  storms <- read.csv(file=stormfiles[i],row.names=NULL,header=T,
                     colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(storms) <- c("ID","storm_suit")
  storms[is.na(storms)] <- 1
  #join together overall suit
  suitdf <- join_all(list(sst,dhw,land,oa,pop,storms),by="ID")
  suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$sst_suit *suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit
  
  #save csv & shapefile
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reefs,suitdf,by="ID")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
  
}

stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("F:/suit/rcp26-ssp1",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*5+2000),.before=T)
  countnew <- rbind(z,countnew)
} 






#####
#####
###historic model
#####
#####

setwd("F:/suit/hist")

#calculate overall suit

#load coral points shapefile & sample raster
reefs <- st_read("E:/currentcoral/reefs.shp")


#load all rasters
#dhw, oa, storms, pop, land, sst, bathy
#check they are all same resolution and extent
#make sure all NAs -> 0's
#recalculate individual variable suits given filter


#use historic model years: 1855-2005
#load sst
sstfiles <- list.files("E:/sst/historic",pattern="_suit.csv$",full.names=T) 
#load dhw
dhwfiles <- list.files("F:/dhwhist",pattern="_suit.csv$",full.names=T) 
#load land
landfiles <- list.files("E:/land_hurtt/historic",pattern="_suit.csv$",full.names=T) 
#load oa
oafiles <- list.files("E:/oa/hist",pattern="_suit.csv$",full.names=T) 
#load pop
popfiles <- list.files("E:/pop_hr/hist_hyde",pattern="_suit.csv$",full.names=T)
#load storms
stormfiles <- list.files("E:/storms/historic",pattern="_suit.csv$",full.names=T)
stormfiles <- stormfiles[-c(1,3,5,7,9,11)] #don't need these years

#name out files
csvout <- paste("hist_",1845+(1:16)*10,"_suit.csv",sep="")
pointsout <- paste("hist",1845+(1:16)*10,"_suit",sep="")



#load each file per year
#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for (i in 1:length(sstfiles)){
  #sst
  sst <- read.csv(file=sstfiles[i],row.names=NULL,header=T,
                  colClasses=c("numeric","NULL","numeric","numeric","NULL","numeric"))
  names(sst) <- c("pointid","iscoral","ID","sst_suit")
  sst[is.na(sst)] <- 1
  #dhw
  dhw <- read.csv(file=dhwfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(dhw) <- c("ID","dhw_suit")
  dhw[is.na(dhw)] <- 1
  #land
  land <- read.csv(file=landfiles[i],row.names=NULL,header=T,
                   colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric","NULL","NULL"))
  names(land) <- c("ID","land_suit")
  land[is.na(land)] <- 1
  #oa
  oa <- read.csv(file=oafiles[i],row.names=NULL,header=T,
                 colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(oa) <- c("ID","oa_suit")
  oa[is.na(oa)] <- 1
  #pop
  pop <- read.csv(file=popfiles[i],row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric"))
  names(pop) <- c("ID","pop_suit")
  pop[is.na(pop)] <- 1
  #if year is >1955, include storms. else don't (storms dataset starts 1955)
  if (i >= 11) {
    storms <- read.csv(file=stormfiles[i-10],row.names=NULL,header=T,
                       colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
    names(storms) <- c("ID","storm_suit")
    storms[is.na(storms)] <- 1
    #join together overall suit
    suitdf <- join_all(list(sst,dhw,land,oa,pop,storms),by="ID")
    suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
    suitdf$suitall <- suitdf$sst_suit *suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit
  } else {
    #join together overall suit
    suitdf <- join_all(list(sst,dhw,land,oa,pop),by="ID")
    suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$iscoral
    suitdf$suitall <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit 
  }
  #save csv & shapefile
  write.csv(suitdf,csvout[i], row.names=FALSE)
  suitpt <- merge(reefs,suitdf,by="ID")
  st_write(suitpt, dsn = pointsout[i], driver = "ESRI Shapefile")
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("F:/suit/hist",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
countcoral<-data.frame()
countnew<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  c <- count(a, vars="suitcoral")
  y <- add_column(c,year=paste0(i*10+1845),.before=T)
  countcoral <- rbind(y,countcoral)
  d <- count(a,vars="suitall")
  z <- add_column(d,year=paste0(i*10+1845),.before=T)
  countnew <- rbind(z,countnew)
} 





#####
#####
###empirical
#####
#####

setwd("F:/suit/emp")

#calculate overall suit

#load coral points shapefile & sample raster
reefs <- st_read("E:/currentcoral/reefs.shp")


#load all dataframes
#dhw, oa, storms, pop, land, sst, bathy
#make sure all NAs -> 0's


#use year 2005

#count overall suit sites for percentage
#csv method - multiply columns together. save csv
#join with shapefile. save

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
  #sst
  sst <- read.csv(file="E:/sst/empirical/hist_2005_suit.csv",row.names=NULL,header=T,
                  colClasses=c("numeric","NULL","numeric","numeric","NULL","numeric"))
  names(sst) <- c("pointid","iscoral","ID","sst_suit")
  sst[is.na(sst)] <- 1
  #dhw
  dhw <- read.csv(file="E:/dhw/historic/empirical/hist_2005_suit.csv",row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(dhw) <- c("ID","dhw_suit")
  dhw[is.na(dhw)] <- 1
  #land
  land <- read.csv(file="E:/land_hurtt/empirical/emp_2005_suit.csv",row.names=NULL,header=T,
                   colClasses=c(rep("NULL",3),"numeric",rep("NULL",2),"numeric"))
  names(land) <- c("ID","land_suit")
  land[is.na(land)] <- 1
  #oa
  oa <- read.csv(file="E:/oa/empirical/emp_2005new_suit.csv",row.names=NULL,header=T,
                 colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(oa) <- c("ID","oa_suit")
  oa[is.na(oa)] <- 1
  #pop
  pop <- read.csv(file="E:/pop_hr/hist_gpw/histgpw_v3_2005_suit.csv",row.names=NULL,header=T,
                  colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"numeric"))
  names(pop) <- c("ID","pop_suit")
  pop[is.na(pop)] <- 1
  #storms
  storms <- read.csv(file="E:/storms/empirical/emp_2005_suit.csv",row.names=NULL,header=T,
                     colClasses=c(rep("NULL",3),"numeric","NULL","numeric"))
  names(storms) <- c("ID","storm_suit")
  storms[is.na(storms)] <- 1
  #join together overall suit
  suitdf <- join_all(list(sst,dhw,land,oa,pop,storms),by="ID")
  suitdf$suitcoral <- suitdf$sst_suit * suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit * suitdf$iscoral
  suitdf$suitall <- suitdf$sst_suit *suitdf$dhw_suit * suitdf$land_suit * suitdf$oa_suit * suitdf$pop_suit * suitdf$storm_suit

  #save csv & shapefile
  write.csv(suitdf,"emp_2005_suit.csv", row.names=FALSE)
  suitpt <- merge(reefs,suitdf,by="ID")
  st_write(suitpt, dsn = "emp_2005_suit", driver = "ESRI Shapefile")
stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
  #suitdf <- read.csv(file="F:/suit/emp/emp_2005_suit.csv",row.names=NULL,header=T)
  countcoral <- count(suitdf, vars="suitcoral")
  countnew <- count(suitdf,vars="suitall")

