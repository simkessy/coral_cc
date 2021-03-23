library(raster)
library(ncdf4)
library(doParallel)
library(velox)
library(sf)
library(foreach)


#### SST BIAS CALCULATION



#NOAA MMM file
noaa.mmm <- brick("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

#open oisst. only keep 1982-2005
#use oisst monthly
oisst <- brick("D:/SST/OISST/oisstmasked.nc")
oisst <- rotate(oisst)
oisst<- dropLayer(oisst,c(1,290:450))

#break it up by months. get mean monthly values
jans <- c(1+12*(0:23))

for(j in 1:12){
  months <- stack(oisst[[jans[1:24]+j-1]])
  calc(months,fun=mean,na.rm=T,filename=paste0("F:/sst/bias/oisst_mean_",j,".tif",sep=""))
}




#open historic models. 

setwd("F:/sst/hist")

parent.folder <- "D:/SST/CMIP5/historical"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(1,6,7,12,19,20,21,23,24)]
for(i in 1:length(n)){
  a <- paste("mh.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}

#adjust rasters so have same temporal range
#clean up datasets
e <- raster(nrow=180,ncol=360,ext=extent(mh.1),crs=crs(mh.1))
e[]<-NA
mh.5<-dropLayer(mh.5,c(1873:1932)) #too many layers. remove the last 60 layers
mh.7<-dropLayer(mh.7,c(1873:1932)) #too many layers. remove the last 60 layers
e120<-stack(replicate(120,e))
mh.9<-stack(e120,mh.9)#add 120
mh.10<-stack(e120,mh.10) #add120
e119<-stack(replicate(119,e))
mh.11<-stack(e119,mh.11)
mh.15[mh.15<270] <- NA
mh.5[mh.5<270] <- NA


models <- lapply(paste0('mh.',1:length(n)),get)
modelnames <- c("CanESM2","CMCC-CESM","CMCC-CM","CMCC-CMS","GISS-E2-H-CC","GISS-E2-H","GISS-E2-R-CC","GISS-E2-R",
                "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","inmcm4","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")


#calculate mean monthly val for 1982-2005
#calculate bias (difference between monthly val and OISST)

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for(i in 1:length(models)){
  mod <- models[[i]]
  for(j in 0:11){
    months <- stack(mod[[jans[1:24]+j+1584]])
    meanmonth <- calc(months,fun=mean,na.rm=T)
    meanmonth[meanmonth <270] <- NA
    i_brickC <- calc(meanmonth, fun=function(x){x - 273.15}) #change from K to C
    i_brickCd <- disaggregate(i_brickC,fact=2)
    i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
    i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
    extent(i_brick1) <- c(-180,-0.5,-90,90)
    i_brickCc <- merge(i_brick1,i_brick2)
    i_brickf <- focal(i_brickCc,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
    oisstmonth<- raster(paste0("F:/sst/bias/oisst_mean_",j+1,".tif",sep=""))
    oisstmonthf <- focal(oisstmonth,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
    oisstd <- disaggregate(oisstmonthf,fact=2)
    bias <- raster::overlay(oisstd,i_brickf,fun=function(x,y,na.rm=T){return((x-y))},filename=paste0("F:/sst/bias/",modelnames[i],"_bias",j+1,".tif",sep=""))
  }
}
stopCluster(cl)
finish <- Sys.time()
finish-start


#adjust bias to models. 
#calculate DHW
#calculate mean year SST. save
bias.folder <- "F:/sst/bias/"
monthlabel <- c(rep(1:12,156))
memory.limit(size=50000000) #memory free = 113524476

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for(i in 1:length(models)){
  #get all 12 bias rasters for each month
  for(j in 1:12){
    a <- paste("bias.", j, sep = "")
    r <- raster(paste0(bias.folder,modelnames[i],"_bias",j,".tif",sep=""))
    assign(a,r)
  }
  biasmonths <- lapply(paste0('bias.',1:12),get)
  #get model into comparable raster format
  i_brickC <- calc(models[[i]], fun=function(x){x - 273.15}) #change from K to C
  i_brickCd <- disaggregate(i_brickC,fact=2)
  i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
  i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
  extent(i_brick1) <- c(-180,-0.5,-90,90)
  i_brickCc <- merge(i_brick1,i_brick2)
  foreach::foreach(l=1:156, .packages=c("raster","foreach")) %dopar% {
  #for(l in 1:156){
    l1 <- 12*(l-1)+1
    l2 <- 12*l
    SSTbc <- foreach::foreach(k=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
      i_brickf <- focal(i_brickCc[[k]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
      raster::overlay(i_brickf,biasmonths[[monthlabel[k]]],fun=function(x,y,na.rm=T){return((x+y))})
    }
    DHWx <- foreach::foreach(h=1:12, .packages=c("raster"), .combine=raster::stack) %dopar% {
      d <- disaggregate(SSTbc[[h]],fact=10)
      DHW <- raster::overlay(d,noaa.mmm,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
      reclassify(DHW, c(-Inf,0,0))
    }
    calc(DHWx,fun=sum,filename=paste0("F:/dhw/hist/",modelnames[i],"/",modelnames[i],"_",1849+l,".tif",sep=""))
    calc(SSTbc,fun=mean,filename=paste0("F:/sst/hist/",modelnames[i],"/",modelnames[i],"_bc_",1849+l,".tif",sep=""))
    gc()
  }
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#30 min per model
#total time = 6.78 hours



#calculate model mean and model median
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:156, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("F:/dhw/hist/",modelnames[j],"/",modelnames[j],"_",1849+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("F:/sst/hist/",modelnames[j],"/",modelnames[j],"_bc_",1849+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("F:/dhw/hist/modelmean/modelmean_bc_",1849+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("F:/dhw/hist/modelmedian/modelmedian_bc_",1849+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("F:/sst/hist/modelmean/modelmean_bc_",1849+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("F:/sst/hist/modelmedian/modelmedian_bc_",1849+i,".tif",sep=""),overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 1.79 hrs

#calculate 10 year mean
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:16, .packages=c("raster")) %dopar% {
  one = i*10-5
  ten = i*10+5
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:ten){
    dhwmr <- raster(paste0("F:/dhw/hist/modelmean/modelmean_bc_",1845+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("F:/dhw/hist/modelmedian/modelmedian_bc_",1845+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("F:/sst/hist/modelmean/modelmean_bc_",1845+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("F:/sst/hist/modelmedian/modelmedian_bc_",1845+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("F:/dhw/hist/modelmean_bc_10y",1845+i*10,".tif",sep=""))
  calc(dhwdb,fun=mean,filename=paste0("F:/dhw/hist/modelmedian_bc_10y",1845+i*10,".tif",sep=""))
  calc(sstmb,fun=mean,filename=paste0("F:/sst/hist/modelmean_bc_10y",1845+i*10,".tif",sep=""))
  calc(sstdb,fun=mean,filename=paste0("F:/sst/hist/modelmedian_bc_10y",1845+i*10,".tif",sep=""))
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 3 min





##rcp85

#open rcp85 models. 

setwd("F:/sst/rcp85")

parent.folder <- "D:/SST/CMIP5/RCP85/grid"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(1,2,7,8,10,11,12,20,21,22,27)]
for(i in 1:length(n)){
  a <- paste("m85.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}

#adjust rasters so have same temporal range
#clean up datasets
m85.11 <- dropLayer(m85.11, c(1141))

models <- lapply(paste0('m85.',1:length(n)),get)
modelnames <- c("CanESM2","CMCC-CESM","CMCC-CM","CMCC-CMS","GISS-E2-H-CC","GISS-E2-H","GISS-E2-R-CC","GISS-E2-R",
                "HadGEM2-AO","HadGEM2-CC","inmcm4","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")





#adjust bias to models. 
#calculate DHW
#calculate mean year SST. save
bias.folder <- "F:/sst/bias/"
monthlabel <- c(rep(1:12,95))
memory.limit(size=50000000) #memory free = 113524476

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for(i in 1:length(models)){
  #get all 12 bias rasters for each month
  for(j in 1:12){
    a <- paste("bias.", j, sep = "")
    r <- raster(paste0(bias.folder,modelnames[i],"_bias",j,".tif",sep=""))
    assign(a,r)
  }
  biasmonths <- lapply(paste0('bias.',1:12),get)
  #get model into comparable raster format
  i_brickC <- calc(models[[i]], fun=function(x){x - 273.15}) #change from K to C
  i_brickCd <- disaggregate(i_brickC,fact=2)
  i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
  i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
  extent(i_brick1) <- c(-180,-0.5,-90,90)
  i_brickCc <- merge(i_brick1,i_brick2)
  foreach::foreach(l=1:95, .packages=c("raster","foreach")) %dopar% {
    l1 <- 12*(l-1)+1
    l2 <- 12*l
    SSTbc <- foreach::foreach(k=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
      i_brickf <- focal(i_brickCc[[k]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
      raster::overlay(i_brickf,biasmonths[[monthlabel[k]]],fun=function(x,y,na.rm=T){return((x+y))})
    }
    DHWx <- foreach::foreach(h=1:12, .packages=c("raster"), .combine=raster::stack) %dopar% {
      d <- disaggregate(SSTbc[[h]],fact=10)
      DHW <- raster::overlay(d,noaa.mmm,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
      reclassify(DHW, c(-Inf,0,0))
    }
    calc(DHWx,fun=sum,filename=paste0("F:/dhw/rcp85/",modelnames[i],"/",modelnames[i],"_",2005+l,".tif",sep=""))
    calc(SSTbc,fun=mean,filename=paste0("F:/sst/rcp85/",modelnames[i],"/",modelnames[i],"_bc_",2005+l,".tif",sep=""))
    gc()
  }
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#46 min per model
#total time = 3.7 hr


#calculate model mean and model median
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:95, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("F:/dhw/rcp85/",modelnames[j],"/",modelnames[j],"_",2005+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("F:/sst/rcp85/",modelnames[j],"/",modelnames[j],"_bc_",2005+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("F:/dhw/rcp85/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("F:/dhw/rcp85/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("F:/sst/rcp85/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("F:/sst/rcp85/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 1.9 hr

#calculate 10 year mean
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  one = i*5-2
  five = i*5+2
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("F:/dhw/rcp85/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("F:/dhw/rcp85/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("F:/sst/rcp85/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("F:/sst/rcp85/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("F:/dhw/rcp85/modelmean_bc_5y",2000+i*5,".tif",sep=""))
  calc(dhwdb,fun=mean,filename=paste0("F:/dhw/rcp85/modelmedian_bc_5y",2000+i*5,".tif",sep=""))
  calc(sstmb,fun=mean,filename=paste0("F:/sst/rcp85/modelmean_bc_5y",2000+i*5,".tif",sep=""))
  calc(sstdb,fun=mean,filename=paste0("F:/sst/rcp85/modelmedian_bc_5y",2000+i*5,".tif",sep=""))
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 1.98 min


##rcp45

#open rcp45 models. 

setwd("F:/sst/rcp45")

parent.folder <- "D:/SST/CMIP5/RCP45/grid"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(1,2,6:12,21:23,26)]
for(i in 1:length(n)){
  a <- paste("m45.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}

#adjust rasters so have same temporal range
#clean up datasets

models <- lapply(paste0('m45.',1:length(n)),get)
modelnames <- c("CanESM2","CMCC-CM","CMCC-CMS","GISS-E2-H-CC","GISS-E2-H","GISS-E2-R-CC","GISS-E2-R",
                "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","inmcm4","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")






#adjust bias to models. 
#calculate DHW
#calculate mean year SST. save
bias.folder <- "F:/sst/bias/"
monthlabel <- c(rep(1:12,95))
memory.limit(size=50000000) #memory free = 113524476

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for(i in 1:length(models)){
  #get all 12 bias rasters for each month
  for(j in 1:12){
    a <- paste("bias.", j, sep = "")
    r <- raster(paste0(bias.folder,modelnames[i],"_bias",j,".tif",sep=""))
    assign(a,r)
  }
  biasmonths <- lapply(paste0('bias.',1:12),get)
  #get model into comparable raster format
  i_brickC <- calc(models[[i]], fun=function(x){x - 273.15}) #change from K to C
  i_brickCd <- disaggregate(i_brickC,fact=2)
  i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
  i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
  extent(i_brick1) <- c(-180,-0.5,-90,90)
  i_brickCc <- merge(i_brick1,i_brick2)
  foreach::foreach(l=1:95, .packages=c("raster","foreach")) %dopar% {
    l1 <- 12*(l-1)+1
    l2 <- 12*l
    SSTbc <- foreach::foreach(k=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
      i_brickf <- focal(i_brickCc[[k]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
      raster::overlay(i_brickf,biasmonths[[monthlabel[k]]],fun=function(x,y,na.rm=T){return((x+y))})
    }
    DHWx <- foreach::foreach(h=1:12, .packages=c("raster"), .combine=raster::stack) %dopar% {
      d <- disaggregate(SSTbc[[h]],fact=10)
      DHW <- raster::overlay(d,noaa.mmm,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
      reclassify(DHW, c(-Inf,0,0))
    }
    calc(DHWx,fun=sum,filename=paste0("F:/dhw/rcp45/",modelnames[i],"/",modelnames[i],"_",2005+l,".tif",sep=""))
    calc(SSTbc,fun=mean,filename=paste0("F:/sst/rcp45/",modelnames[i],"/",modelnames[i],"_bc_",2005+l,".tif",sep=""))
    gc()
  }
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#46 min per model
#total time = 


#calculate model mean and model median
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("F:/dhw/rcp45/",modelnames[j],"/",modelnames[j],"_",2005+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("F:/sst/rcp45/",modelnames[j],"/",modelnames[j],"_bc_",2005+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("F:/dhw/rcp45/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("F:/dhw/rcp45/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("F:/sst/rcp45/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("F:/sst/rcp45/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 1.895363 hours

#calculate 5 year mean
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  one <- i*5-2
  five <- i*5+2
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("F:/dhw/rcp45/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("F:/dhw/rcp45/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("F:/sst/rcp45/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("F:/sst/rcp45/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("F:/dhw/rcp45/modelmean_bc_5y",2000+i*5,".tif",sep=""))
  calc(dhwdb,fun=mean,filename=paste0("F:/dhw/rcp45/modelmedian_bc_5y",2000+i*5,".tif",sep=""))
  calc(sstmb,fun=mean,filename=paste0("F:/sst/rcp45/modelmean_bc_5y",2000+i*5,".tif",sep=""))
  calc(sstdb,fun=mean,filename=paste0("F:/sst/rcp45/modelmedian_bc_5y",2000+i*5,".tif",sep=""))
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 2 min




##rcp26

#open rcp26 models. 

setwd("F:/sst/rcp26")

parent.folder <- "D:/SST/CMIP5/RCP26/grid"
m <- list.files(parent.folder, pattern=".nc",full.names=T)
n <- m[-c(2:6,10,11,12)]
for(i in 1:length(n)){
  a <- paste("m26.", i, sep = "")
  r <- brick(n[i])
  assign(a,r)
}

#use NOAA MMM file
noaa.mmm <- raster("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

models <- lapply(paste0('m26.',1:length(n)),get)
modelnames <- c("CanESM2","GISS-E2-H","GISS-E2-R","HadGEM2-AO","MIROC-ESM-CHEM","MIROC-ESM","MRI-CGCM3")





#adjust bias to models. 
#calculate DHW
#calculate mean year SST. save
bias.folder <- "F:/sst/bias/"
monthlabel <- c(rep(1:12,95))
memory.limit(size=50000000) #memory free = 113524476

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)
for(i in 1:length(models)){
  #get all 12 bias rasters for each month
  for(j in 1:12){
    a <- paste("bias.", j, sep = "")
    r <- raster(paste0(bias.folder,modelnames[i],"_bias",j,".tif",sep=""))
    assign(a,r)
  }
  biasmonths <- lapply(paste0('bias.',1:12),get)
  #get model into comparable raster format
  i_brickC <- calc(models[[i]], fun=function(x){x - 273.15}) #change from K to C
  i_brickCd <- disaggregate(i_brickC,fact=2)
  i_brick1 <- crop(i_brickCd,extent(180,359.5,-90,90))
  i_brick2 <- crop(i_brickCd,extent(-0.5,180,-90,90))
  extent(i_brick1) <- c(-180,-0.5,-90,90)
  i_brickCc <- merge(i_brick1,i_brick2)
  foreach::foreach(l=1:95, .packages=c("raster","foreach")) %dopar% {
    l1 <- 12*(l-1)+1
    l2 <- 12*l
    SSTbc <- foreach::foreach(k=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
      i_brickf <- focal(i_brickCc[[k]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T)
      raster::overlay(i_brickf,biasmonths[[monthlabel[k]]],fun=function(x,y,na.rm=T){return((x+y))})
    }
    DHWx <- foreach::foreach(h=1:12, .packages=c("raster"), .combine=raster::stack) %dopar% {
      d <- disaggregate(SSTbc[[h]],fact=10)
      DHW <- raster::overlay(d,noaa.mmm,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
      reclassify(DHW, c(-Inf,0,0))
    }
    calc(DHWx,fun=sum,filename=paste0("F:/dhw/rcp26/",modelnames[i],"/",modelnames[i],"_",2005+l,".tif",sep=""))
    calc(SSTbc,fun=mean,filename=paste0("F:/sst/rcp26/",modelnames[i],"/",modelnames[i],"_bc_",2005+l,".tif",sep=""))
    gc()
  }
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#46 min per model
#total time = 


#calculate model mean and model median
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:95, .packages=c("raster")) %dopar% {
  dhwb <- raster()
  sstb <- raster()
  for(j in 1:length(modelnames)){
    dhwr <- raster(paste0("F:/dhw/rcp26/",modelnames[j],"/",modelnames[j],"_",2005+i,".tif",sep=""))
    dhwb <- stack(dhwb,dhwr)
    sstr <- raster(paste0("F:/sst/rcp26/",modelnames[j],"/",modelnames[j],"_bc_",2005+i,".tif",sep=""))
    sstb <- stack(sstb,sstr)
  }
  calc(dhwb,fun=mean,na.rm=T,filename=paste0("F:/dhw/rcp26/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(dhwb,fun=median,na.rm=T,filename=paste0("F:/dhw/rcp26/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=mean,na.rm=T,filename=paste0("F:/sst/rcp26/modelmean/modelmean_bc_",2005+i,".tif",sep=""),overwrite=T)
  calc(sstb,fun=median,na.rm=T,filename=paste0("F:/sst/rcp26/modelmedian/modelmedian_bc_",2005+i,".tif",sep=""),overwrite=T)
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 2 hr

#calculate 5 year mean
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
  one = i*5-2
  five = i*5+2
  dhwmb <- raster()
  dhwdb <- raster()
  sstmb <- raster()
  sstdb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("F:/dhw/rcp26/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
    dhwdr <- raster(paste0("F:/dhw/rcp26/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    dhwdb <- stack(dhwdr,dhwdb)
    sstmr <- raster(paste0("F:/sst/rcp26/modelmean/modelmean_bc_",2000+j,".tif",sep=""))
    sstmb <- stack(sstmr,sstmb)
    sstdr <- raster(paste0("F:/sst/rcp26/modelmedian/modelmedian_bc_",2000+j,".tif",sep=""))
    sstdb <- stack(sstdr,sstdb)
  }
  calc(dhwmb,fun=mean,filename=paste0("F:/dhw/rcp26/modelmean_bc_5y",2000+i*5,".tif",sep=""))
  calc(dhwdb,fun=mean,filename=paste0("F:/dhw/rcp26/modelmedian_bc_5y",2000+i*5,".tif",sep=""))
  calc(sstmb,fun=mean,filename=paste0("F:/sst/rcp26/modelmean_bc_5y",2000+i*5,".tif",sep=""))
  calc(sstdb,fun=mean,filename=paste0("F:/sst/rcp26/modelmedian_bc_5y",2000+i*5,".tif",sep=""))
}
stopCluster(cl)
finish <- Sys.time()
finish-start
#total time = 2 min




####
###empirical

setwd("F:/dhw/emp")

#use NOAA MMM file
noaa.mmm <- brick("D:/SST/CMIP5/RCP26/ct5km_climatology_v3.1.nc", varname = "sst_clim_mmm", lvar=4)
crs(noaa.mmm) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "

#use oisst monthly
oisst <- brick("D:/SST/OISST/oisstmasked.nc")
oisst <- rotate(oisst)
oisst<- dropLayer(oisst,c(1,446:450))

#calculate mean year SST. save
monthlabel <- c(rep(1:12,37))
memory.limit(size=50000000) #memory free = 113524476

start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(l=1:37, .packages=c("raster","foreach")) %dopar% {
  l1 <- 12*(l-1)+1
  l2 <- 12*l
  DHWx <- foreach::foreach(h=l1:l2, .packages=c("raster"), .combine=raster::stack) %dopar% {
    d <- disaggregate(oisst[[h]],fact=20)
    DHW <- raster::overlay(d,noaa.mmm,fun=function(x,y,na.rm=T){return((x-y)*4.34)})
    reclassify(DHW, c(-Inf,0,0))
  }
  calc(DHWx,fun=sum,filename=paste0("F:/dhw/emp/annual/emp_dhw_",1981+l,".tif",sep=""))
  gc()
}

stopCluster(cl)
finish <- Sys.time()
finish-start
#6.4 min



#calculate 5 year mean
start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

foreach::foreach(i=1:7, .packages=c("raster")) %dopar% {
  one = i*5-2
  five = i*5+2
  dhwmb <- raster()
  for(j in one:five){
    dhwmr <- raster(paste0("F:/dhw/emp/annual/emp_dhw_",1980+j,".tif",sep=""))
    dhwmb <- stack(dhwmb,dhwmr)
  }
  calc(dhwmb,fun=mean,filename=paste0("F:/dhw/emp/emp_oisst_5y",1980+i*5,".tif",sep=""))
}
stopCluster(cl)
finish <- Sys.time()
finish-start



