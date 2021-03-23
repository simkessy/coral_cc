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
library(plotrix)
library(caret)

#define threshold
#load all noaa climatology files
noaa.1 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_january")
noaa.2 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_february")
noaa.3 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_march")
noaa.4 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_april")
noaa.5 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_may")
noaa.6 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_june")
noaa.7 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_july")
noaa.8 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_august")
noaa.9 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_september")
noaa.10 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_october")
noaa.11 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_november")
noaa.11[noaa.11 > 300] <- NA #get rid of 3 cells that have values 326 and 327
noaa.12 <- raster("E:/ct5km_climatology_v3.1.nc",varname="sst_clim_december")

#get mean of climatology files
noaabrick <- brick(noaa.1,noaa.2,noaa.3,noaa.4,noaa.5,noaa.6,noaa.7,noaa.8,noaa.9,noaa.10,noaa.11,noaa.12)
crs(noaabrick) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#find sst values at coral sites
reefs <- st_read("E:/currentcoral/reefs.shp")
reefext <- rgis::fast_extract(sf=reefs,ras=noaabrick,funct=function(x){max(x,na.rm=T)},parallel=T)


#get 95% confidence interval. these are our lower and upper threshold limits for sst


###RCP26

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/rcp26")


#reef <- st_read("D:/distribution/all_reef_pt.shp")
reefs <- st_read("E:/currentcoral/reefs.shp")
sst <- brick("E:/sst/sstRCP26med_mean5.tif")
sstr <- rotate(sst)
#rcp26rout <- paste("rcp26_",2000+(1:20)*5,"_sst.tif",sep="")
#rcp26suitout <- paste("rcp26_",2000+(1:20)*5,"_suit.tif",sep="")
rcp26suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

  reefext <- rgis::fast_extract(sf=reefs,ras=sstr,funct=function(x){max(x,na.rm=T)},parallel=T)
  dfext <- as.data.frame(reefext)
  for (i in 1:20){
    x <- paste0("sstr_sstRCP26med_mean5.",i)
    dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID", "geometry",x))
    names(dfi) <- c("pointid","grid_code", "iscoral","ID", "geometry","sst")
    dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
    #sstpt <- inner_join(reef,dfi,by="pointid")
    #sstdf <- as.data.frame(sstpt)
    write.csv(dfi,rcp26suitcsvout[i], row.names=FALSE)
    #rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=rcp26rout[i],overwrite=T)
    #rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp26suitout[i],overwrite=T)
  }

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/sst/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=F,skip=1,colClasses=c(rep("numeric",4),rep("NULL",2),"numeric","numeric"))
  names(a) <- c("pointid","grid_code","iscoral","ID","sst","suit")
  a$sst_coral <- a$suit * a$iscoral
  c <- count(a, vars=sst_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 




###RCP45

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/rcp45")


#reef <- st_read("D:/distribution/all_reef_pt.shp")
reefs <- st_read("E:/currentcoral/reefs.shp")
sst <- brick("E:/sst/sstRCP45med_mean5.tif")
sstr <- rotate(sst)
#rcp45rout <- paste("rcp45_",2000+(1:20)*5,"_sst.tif",sep="")
#rcp45suitout <- paste("rcp45_",2000+(1:20)*5,"_suit.tif",sep="")
rcp45suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")
#extent(r008333) <- c(0,360,-90,90)
#res(r008333)<- 0.008333333

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=sstr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("sstr_sstRCP45med_mean5.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
  #sstpt <- inner_join(reef,dfi,by="pointid")
  #sstdf <- as.data.frame(sstpt)
  write.csv(dfi,rcp45suitcsvout[i], row.names=FALSE)
  #rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=rcp45rout[i],overwrite=T)
  #rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp45suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/sst/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$sst_coral <- a$suit * a$iscoral
  c <- count(a, vars=sst_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 






###RCP85

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/rcp85")


reefs <- st_read("E:/currentcoral/reefs.shp")
sst <- brick("E:/sst/sstRCP85med_mean5.tif")
sstr <- rotate(sst)
#rcp85rout <- paste("rcp85_",2000+(1:20)*5,"_sst.tif",sep="")
#rcp85suitout <- paste("rcp85_",2000+(1:20)*5,"_suit.tif",sep="")
rcp85suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")
#extent(r008333) <- c(0,360,-90,90)
#res(r008333)<- 0.008333333

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=sstr,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:20){
  x <- paste0("sstr_sstRCP85med_mean5.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
  #sstpt <- inner_join(reef,dfi,by="pointid")
  #sstdf <- as.data.frame(sstpt)
  write.csv(dfi,rcp85suitcsvout[i], row.names=FALSE)
  #rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=rcp85rout[i],overwrite=T)
  #rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=rcp85suitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/sst/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$sst_coral <- a$suit * a$iscoral
  c <- count(a, vars=sst_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 




###empirical ### 

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/empirical")

sst <- brick("D:/SST/OISST/oisstmasked.nc")
sstr <- rotate(sst)


#make brick with 5 year average. 
#raw file has monthly dates 1981/12 to 2019/6 (450 layers). make averages for 1985, 1990, 1995, 2000, 2005, 2010, 2015

sstd <- dropLayer(sstr, c(1:13,433:450)) #drop layers. we don't need the year 1982 or 2018, 2019
mvlist <- unstack(sstd) # now a list of rasters 
grpsize <- 60 # desired layers per stack - 60 months or five years
# uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
names(mvlist) <- rep(1:(ceiling(length(mvlist) / grpsize )), 
                     each = grpsize, length.out = length(mvlist))
# make list of rasters into a list of stacks. basically each year is now one separate stack
stacks <- lapply(unique(names(mvlist)), function(y) {
  b <- mvlist[names(mvlist) == y]
  stack(b)
})
# find the mean temp for 5 years - w is each stack/5year. 
stacksout <- lapply(stacks, function(w) {
  g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
})
sst_five <- stack(stacksout)
writeRaster(sst_five, filename="oisst_fiveyr.tif", overwrite=TRUE)



reefs <- st_read("E:/currentcoral/reefs.shp")
sst <- brick("oisst_fiveyr.tif")
#histrout <- paste("hist_",1980+(1:7)*5,"_sst.tif",sep="")
#histsuitout <- paste("hist_",1980+(1:7)*5,"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1980+(1:7)*5,"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=sst,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:7){
  x <- paste0("sst_oisst_fiveyr.",i)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 24, ifelse(sst <= 32,1,0), 0)) #in celsius. 24 to 32C. 297 to 305K
  #sstpt <- inner_join(reef,dfi,by="pointid")
  #sstdf <- as.data.frame(sstpt)
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  #rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  #rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/sst/empirical",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$sst_coral <- a$suit * a$iscoral
  c <- count(a, vars=sst_coral)
  y <- add_column(c,year=paste0(i*5+1945),.before=T)
  count <- rbind(y,count)
} 





#####
#####calculate historical scenario model mean/median, make taylor diagram using mean of overlapping period 1981-2005
parent.folder <- "D:/SST/CMIP5/historical"


m <- list.files(parent.folder, pattern=".nc",full.names=T)
for(i in 1:length(m)){
  a <- paste("m", i, sep = "")
  r <- brick(m[i])
  assign(a,r)
}
#add xx empty rasters to mx so it matches the rest with 1855-2005
e <- raster(nrow=180,ncol=360,ext=extent(m1),crs=crs(m1))
e[]<-NA
e1332<-stack(replicate(1332,e))
m1e<-stack(e1332,m1)#add 1332 empty rasters because this one is only 1961-2005
m8d<-dropLayer(m8,c(1873:1932)) #too many layers. remove the last 60 layers
m10d<-dropLayer(m10,c(1873:1932)) #too many layers. remove the last 60 layers
e119<-stack(replicate(119,e))
m12e<-stack(e119,m12) #add 119 layers
e120<-stack(replicate(120,e))
m13e<-stack(e120,m13)#add 120
m14e<-stack(e120,m14) #add120
m15e<-stack(e119,m15)
#m19[m19<270] <- NA #bad dataset. don't use it. has some bizarre interpolation going on with an unmasked land (MIROC5)
#m19d<-dropLayer(m19,c(1873:1956))
m22[m22<270] <- NA
m23[m23<270] <- NA
e12<-stack(replicate(12,e))
m23e<-stack(e12,m23)

#mean&median calculation
modelm <- overlay(m1e,m2,m3,m4,m5,m6,m7,m8d,m9,m10d,m11,m12e,m13e,m14e,m15e,m16,m17,m18,m20,m21,m22,m23e,m24,fun=function(x){ mean(x,na.rm=T)},
                  filename="E:/sst/SSThist_modelmean.tif")
modelmed <-overlay(m1e,m2,m3,m4,m5,m6,m7,m8d,m9,m10d,m11,m12e,m13e,m14e,m15e,m16,m17,m18,m20,m21,m22,m23e,m24,fun=function(x){ median(x,na.rm=T)},
                   filename="E:/sst/SSThist_modelmedian.tif")




#####
#TEST: DONT USE HIST MODELS THAT DONT HAVE A RCP MODEL
#####calculate historical scenario model mean/median, make taylor diagram using mean of overlapping period 1981-2005
parent.folder <- "D:/SST/CMIP5/historical"


m <- list.files(parent.folder, pattern=".nc",full.names=T)
for(i in 1:length(m)){
  a <- paste("m", i, sep = "")
  r <- brick(m[i])
  assign(a,r)
}

#don't use CanCM4, CNRM-CM5, HadCM3, MPI-ESM-LR, MPI-ESM-MR, NorESM1-ME
#m1,m6,m12,m20,m21,m24 ... also m19 bc weird interpolation at coastlines

#add xx empty rasters to mx so it matches the rest with 1855-2005
e <- raster(nrow=180,ncol=360,ext=extent(m1),crs=crs(m1))
e[]<-NA
m8d<-dropLayer(m8,c(1873:1932)) #too many layers. remove the last 60 layers
m10d<-dropLayer(m10,c(1873:1932)) #too many layers. remove the last 60 layers
e119<-stack(replicate(119,e))
e120<-stack(replicate(120,e))
m13e<-stack(e120,m13)#add 120
m14e<-stack(e120,m14) #add120
m15e<-stack(e119,m15)
m22[m22<270] <- NA
m23[m23<270] <- NA
e12<-stack(replicate(12,e))
m23e<-stack(e12,m23)

#mean&median calculation
modelm <- overlay(m2,m3,m4,m5,m7,m8d,m9,m10d,m11,m13e,m14e,m15e,m16,m17,m18,m22,m23e,fun=function(x){ mean(x,na.rm=T)},
                  filename="E:/sst/SSThist_modelmean_n.tif")
modelmed <-overlay(m2,m3,m4,m5,m7,m8d,m9,m10d,m11,m13e,m14e,m15e,m16,m17,m18,m22,m23e,fun=function(x){ median(x,na.rm=T)},
                   filename="E:/sst/SSThist_modelmedian_n.tif")







###
#empirical & hist taylor diagrams

#calculate accuracy: Taylor diagrams

#use overlapping years with historical scenarios (1981/12 to 2005/12)
empirical <- brick("D:/SST/OISST/oisstmasked.nc")
empd <- dropLayer(empirical,c(290:450))
empdm<- calc(empd,fun=mean,na.rm=T)
empext <- rgis::fast_extract(sf=reefs,ras=empdm,funct=function(x){max(x,na.rm=T)},parallel=T)
empdfext <- as.data.frame(empext)
empdfext$oisst <- empdfext$empdm_layer+273.15
empdf <- empdfext[,c(3,4,7)]
names(empdf) <- c("iscoral","pointid","oisst")



##get mean values for each model, model mean, model median for same test period - 1981-2005
##extract only values at coral points
rasterList <- list(m1e,m2,m3,m4,m5,m6,m7,m8d,m9,m10d,m11,m12e,m13e,m14e,m15e,m16,m17,m18,m20,m21,m22,m23e,m24)
rasterListnew <- list(m2,m3,m4,m5,m7,m8d,m9,m10d,m11,m13e,m14e,m15e,m16,m17,m18,m22,m23e)
for(i in 1:length(rasterListnew)){
  r<-rasterListnew[[i]]
  rd <- dropLayer(r,c(1:1583))
  rm<- calc(rd,fun=mean,na.rm=T)
  reefext <- rgis::fast_extract(sf=reefs,ras=rm,funct=function(x){max(x,na.rm=T)},parallel=T)
  dfext <- as.data.frame(reefext)
  dfext$SSTC <- dfext$rm_layer-273.15
  rdf <- dfext[,c(4,5)]
  names(rdf) <- c("pointid",paste("model",i,sep=""))
  a <- paste("mod", i, sep = "")
  assign(a,rdf)
}

rd <- dropLayer(modelm,c(1:1583))
rm<- stackApply(rd,indices=rep(1,nlayers(rd)),fun="mean",na.rm=T) #raster::calc(meand,fun=mean,na.rm=T) #this is resulting in an error for some reason
reefext <- rgis::fast_extract(sf=reefs,ras=rm,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
dfext$SSTC <- dfext$rm_index_1-273.15
mndf <- dfext[,c(4,5)]
names(mndf) <- c("pointid","meanmod")

rd <- dropLayer(modelmed,c(1:1583))
rm<- stackApply(rd,indices=rep(1,nlayers(rd)),fun="mean",na.rm=T) #raster::calc(meand,fun=mean,na.rm=T) #this is resulting in an error for some reason
reefext <- rgis::fast_extract(sf=reefs,ras=rm,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
dfext$SSTC <- dfext$rm_index_1-273.15
medf <- dfext[,c(4,5)]
names(medf) <- c("pointid","medianmod")



##join all (by pointid) into one dataframe for comparison. filter to just coral sites
alldata <- plyr::join_all(list(empdf,mndf,medf,mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,
                         mod15,mod16,mod17),by="pointid")
names(alldata) <- c("iscoral","pointid","oisst","meanmod","medianmod","model1","model2","model3","model4","model5","model6","model7","model8","model9","model10",
                    "model11","model12","model13","model14","model15","model16","model17")
alldf<-alldata[!(alldata$iscoral==0),]
#noNAdata <- na.omit(alldata)





oldpar<-taylor.diagram(alldf$oisst,alldf$oisst,ref.sd=T,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="blue",gamma.col="blue",
                       main="Taylor Diagram - SST historical, coral sites")
oldpar<-taylor.diagram(alldf$oisst,alldf$model1,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model2,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model3,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model4,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model5,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model6,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model7,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model8,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model9,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model10,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model11,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model12,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model13,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model14,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model15,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model16,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model17,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model18,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model20,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model21,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model22,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model23,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$model24,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$oisst,alldf$meanmod,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,col="red",pch=19)
oldpar<-taylor.diagram(alldf$oisst,alldf$medianmod,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,col="green",pch=19)



###historical model ### 

####
#extract reef suit vals
#find sst for each point
#classify as suitable or not - save as .csv
###

setwd("E:/sst/historic")

sst <- brick("E:/sst/SSThist_modelmedian_n.tif")
filen <- "ssthistmed_mean10_n.tif"

#make brick with 5 year average. 
#raw file has monthly dates 1981/12 to 2019/6 (450 layers). make averages for 1985, 1990, 1995, 2000, 2005, 2010, 2015

#calculate ten yr mean of historic model for suitability analysis
tenyrmeanMt <- function(rasterxyz,filen) {
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  tenyr <- 120 # desired layers per stack -  10 years (or 120 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  names(meanz_years) <- rep(1:(ceiling(length(meanz_years) / tenyr )), 
                            each = tenyr, length.out = length(meanz_years))
  # make list of rasters into a list of stacks. basically each decade is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  # find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- stackApply(w,indices=rep(1,nlayers(w)),fun="mean",na.rm=T)
    #g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_tenyr <- stack(stacksout)
  names(meanz_tenyr) <- c(paste(1845+10*(1:16)))
  writeRaster(meanz_tenyr, filename=filen, format="GTiff", overwrite=TRUE)
  
}

tenyrmeanMt(sst,filen)


###calculate historical model suitability
sst <- brick("E:/sst/historic/ssthistmed_mean10_n.tif")
reefs <- st_read("E:/currentcoral/reefs.shp")
reef360 <- st_shift_longitude(reefs) #make longitude 0-360
#histrout <- paste("hist_",1845+10*(1:16),"_sst.tif",sep="")
#histsuitout <- paste("hist_",1845+10*(1:16),"_suit.tif",sep="")
histsuitcsvout <- paste("hist_",1845+10*(1:16),"_suit.csv",sep="")
#r008333 <- raster("E:/pop_hr/ssp1/ssp1_2010_dens.tif")
#extent(r008333) <- c(0,360,-90,90)
#res(r008333)<- 0.008333333


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reef360,ras=sst,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
for (i in 1:16){
  n <- i*10+1845
  x <- paste0("sst_X",n)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","sst")
  dfi$suit <- with(dfi, ifelse(sst >= 297, ifelse(sst <= 305,1,0), 0)) #in kelvin. 24 to 32C. 297 to 305K
  #sstpt <- inner_join(reef360,dfi,by="pointid")
  #sstdf <- as.data.frame(sstpt)
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  #rasterize(sstpt,r008333,field="sst",fun=max,na.rm=T,filename=histrout[i],overwrite=T)
  #rasterize(sstpt,r008333,field="suit",fun=max,na.rm=T,filename=histsuitout[i],overwrite=T)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/sst/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$sst_coral <- a$suit * a$iscoral
  c <- count(a, vars="sst_coral")
  y <- add_column(c,year=paste0(i*10+1845),.before=T)
  count <- rbind(y,count)
} 





