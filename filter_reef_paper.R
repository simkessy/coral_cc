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

#make sst filter - add together suit for rcp85 2010 & 2100
csvfiles <- list.files("E:/sst/historic",pattern="suit.csv$",full.names=T)
  sst1955 <- read.csv(file=csvfiles[11],row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",4),"numeric"))
  names(sst1955) <- c("pointid","sst1955suit")
  sst1955$ID <- seq.int(nrow(sst1955))
  

csvfiles <- list.files("E:/sst/rcp85",pattern="suit.csv$",full.names=T)
  sst2100 <- read.csv(file=csvfiles[20],row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",4),"numeric"))
  names(sst2100) <- c("pointid","sst2100suit")
  sst2100$ID <- seq.int(nrow(sst2100))
  

sstfilter <- merge(sst1955,sst2100,by="ID")
sstfilter$sst1955suit[is.na(sstfilter$sst1955suit)] <- 0
sstfilter$sst2100suit[is.na(sstfilter$sst2100suit)] <- 0
sstfilter$sstsuit <- sstfilter$sst1955suit * sstfilter$sst2100suit
count(sstfilter, vars=sstsuit)


#load up all the tables we want to filter (bathy, sst, coral)
#make sure all NA -> 0
  bathy <- read.csv(file="E:/bathy/bathy.csv",row.names=NULL,header=T,colClasses=c("numeric",rep("NULL",4),"numeric"))
  names(bathy) <- c("pointid","bathysuit")
  count(bathy, vars=bathysuit)
  bathy$ID <- seq.int(nrow(bathy))
  
#merge the tables & multiply them all together into new column filter
  filter <- merge(sstfilter,bathy, by="ID")
  filter$sstbathy <- filter$sstsuit * filter$bathysuit
  count(filter, vars=sstbathy)

#check filter - make sure all current coral sites are included. which percent are current coral sites?
  currentcoral <- read.csv(file="E:/currentcoral/coral.csv",row.names=NULL,header=F,skip=1,colClasses=c("numeric",rep("NULL",2),"numeric",rep("NULL",4)))
  names(currentcoral)<-c("pointid","iscoral")
  currentcoral$iscoral[is.na(currentcoral$iscoral)] <- 0
  count(currentcoral,vars=iscoral)

  coralfilter <- merge(filter,currentcoral,by="pointid")  
  coralfilter$iscoraltest <- coralfilter$iscoral * 10
  coralfilter$test <- coralfilter$sstbathy + coralfilter$iscoraltest
  count(coralfilter, vars=test)
  coralfilter$sstbathycoral <- ifelse(coralfilter$test>=1,1,0)
  count(coralfilter, vars=sstbathycoral)
  #approx 23% of potential sites are current coral sites
  
#make the csv
  write.csv(coralfilter,"filter.csv", row.names=FALSE)
  