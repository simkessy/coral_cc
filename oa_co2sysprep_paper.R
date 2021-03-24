#prep cmip5 data for use in co2sys
library(raster)
library(ncdf4)
library(rgdal)
library(plyr)
library(plotrix)

#prep each variable: make each into 5 year averages
#for each experiment (rcp26,45,85,hist), make a dataframe that has all needed columns to copy into CO2SYS
#columns: x, y, salinity, SST (c), pressure, totalP, totalSi, SST, Pressure, TA, TCO2

##function to calculate 5 year average from yearly data + units adjustment for dissic##
fiveyrmeanYd <- function(rasterxyz) {
  #unstack to groups of 5, average, then restack 
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  fiveyr <- 5 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(0,0,rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                    each = fiveyr, length.out = length(meanz_years)))
  list <- list[-c(95:96)]
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each 5 year period is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  #find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_fiveyr <- brick(stacksout)
  names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                           '2065','2070','2075','2080','2085','2090','2095','2100')
  dissic_mean_5yr <- calc(meanz_fiveyr, fun=function(x){x*1000/1.025}) #change from mol/m3 to umol/kgSW. x*1000/1.025
  dissic_mean_5yr <- rotate(dissic_mean_5yr)
  extent(dissic_mean_5yr) <- ext
  writeRaster(dissic_mean_5yr, filename=rasternewout[1], format="GTiff", overwrite=TRUE)
  
}

##so
fiveyrmeanMs <- function(rasterxyz) {
  #unstack to groups of 5 (e.g. 2006-2010,2011-2015, 2016-2020, 2020-2025), average, then restack 
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  fiveyr <- 60 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(rep(0,24),rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                    each = fiveyr, length.out = length(meanz_years)))
  list <- list[-c(1141:1164)]
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each 5 year period is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  #find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_fiveyr <- brick(stacksout)
  names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                           '2065','2070','2075','2080','2085','2090','2095','2100')
  meanz_fiveyr <- rotate(meanz_fiveyr)
  extent(meanz_fiveyr) <- ext
  writeRaster(meanz_fiveyr, filename=rasternewout[2], format="GTiff", overwrite=TRUE)
  
}

##talk
fiveyrmeanYt <- function(rasterxyz) {
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  fiveyr <- 5 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(0,0,rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                    each = fiveyr, length.out = length(meanz_years)))
  list <- list[-c(95:96)]
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each 5 year period is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  #find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_fiveyr <- brick(stacksout)
  names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                           '2065','2070','2075','2080','2085','2090','2095','2100')
  talk_mean_5yr <- calc(meanz_fiveyr, fun=function(x){x*1000/1.025}) #change from mol/m3 to umol/kgSW. x*1000/1.025
  talk_mean_5yr <- rotate(talk_mean_5yr)
  extent(talk_mean_5yr) <- ext
  writeRaster(talk_mean_5yr, filename=rasternewout[3], format="GTiff", overwrite=TRUE)
  
}

##function to calculate 5 year average from monthly data + units adjust K to C for tos ##
fiveyrmeanMt <- function(rasterxyz) {
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  fiveyr <- 60 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(rep(0,24),rep(1:(ceiling((length(meanz_years)-2) / fiveyr )), 
                          each = fiveyr, length.out = length(meanz_years)))
  list <- list[-c(1141:1164)]
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each decade is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  # find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_fiveyr <- stack(stacksout)
  names(meanz_fiveyr) <- c('2005','2010','2015','2020','2025','2030','2035','2040','2045','2050','2055','2060',
                           '2065','2070','2075','2080','2085','2090','2095','2100')
  sst_mean_5yr <- calc(meanz_fiveyr, fun=function(x){x-273.15})
  sst_mean_5yr <- rotate(sst_mean_5yr)
  extent(sst_mean_5yr) <- ext
  writeRaster(sst_mean_5yr, filename=rasternewout[4], format="GTiff", overwrite=TRUE)
  
}

####################
#silicate (WOA2018)
####################
sili <- raster("D:/omega_arag/co2sys_input/silicate_woa18_all_i00_01.nc", varname = "i_an", lvar=3,level=1)
ext <- extent(sili)


######
#phosphate (WOA2018)
phos <- raster("D:/omega_arag/co2sys_input/phosphate_woa18_all_p00_01.nc",varname="p_an",lvar=3,level=1)






###RCP85
#
parent.folder <- "E:/oa/co2sys_input/rcp85"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,full.names=T)
  r1 <- brick(f[1],varname="dissic",level=1, lvar=3)
  r2 <- brick(f[2],varname="so",level=1, lvar=3)  
  r3 <- brick(f[3],varname="talk",level=1, lvar=3)
  r4 <- brick(f[4], varname = "tos", lvar=4)
  #files are always in alphabetical order. run fiveyrmeanYd on dissic[1], fiveyrmeanMs on so[2], fiveyrmeanYt on talk[3], fiveyrmeanMt on tos[4]
  rasternewout <- paste((substr(f[1:4],1,nchar(f[1:4])-3)),".tif",sep="")
  fiveyrmeanYd(r1)
  fiveyrmeanMs(r2)
  fiveyrmeanYt(r3)
  fiveyrmeanMt(r4)
}


#create csv files
#columns: x, y, salinity, SST (c), pressure, totalP, totalSi, SST, Pressure, TA, TCO2
parent.folder <- "E:/oa/co2sys_input/rcp85"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
csvfnames <-paste("rcp85_",2000+(1:20)*5,"_co2sysin.csv",sep="")
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,pattern=".tif$",full.names=T)
  r1 <- brick(f[1],varname="dissic",level=1, lvar=3)
  r2 <- brick(f[2],varname="so",level=1, lvar=3)  
  r3 <- brick(f[3],varname="talk",level=1, lvar=3)
  r4 <- brick(f[4], varname = "tos", lvar=4)
  for (i in 1:20){
    dissic <- r1[[i]]
    so <- r2[[i]]
    tos <- r4[[i]]
    talk <- r3[[i]]
    r <- brick(so,tos,phos,sili,talk,dissic)
    names(r)=c("salt","sst","phos","sili","talk","inorgC")
    sysdf <- as.data.frame(r,xy=TRUE)
    sysdf$pressure <- 0
    sysdf$Pressure <- 0
    sysdf$SST <- sysdf$sst
    cosysdf <- sysdf[,c(1,2,3,4,9,5,6,11,10,7,8)]
    write.csv(cosysdf,csvfnames[i],na="")
  }
}





###RCP45
#
parent.folder <- "E:/oa/co2sys_input/rcp45"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,full.names=T)
  r1 <- brick(f[1],varname="dissic",level=1, lvar=3)
  r2 <- brick(f[2],varname="so",level=1, lvar=3)  
  r3 <- brick(f[3],varname="talk",level=1, lvar=3)
  r4 <- brick(f[4], lvar=4)
  #files are always in alphabetical order. run fiveyrmeanYd on dissic[1], fiveyrmeanMs on so[2], fiveyrmeanYt on talk[3], fiveyrmeanMt on tos[4]
  rasternewout <- paste((substr(f[1:4],1,nchar(f[1:4])-3)),".tif",sep="")
  fiveyrmeanYd(r1)
  fiveyrmeanMs(r2)
  fiveyrmeanYt(r3)
  fiveyrmeanMt(r4)
}


#create csv files
#columns: x, y, salinity, SST (c), pressure, totalP, totalSi, SST, Pressure, TA, TCO2
parent.folder <- "E:/oa/co2sys_input/rcp45"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
csvfnames <-paste("rcp45_",2000+(1:20)*5,"_co2sysin.csv",sep="")
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,pattern=".tif$",full.names=T)
  r1 <- brick(f[1],varname="dissic",level=1, lvar=3)
  r2 <- brick(f[2],varname="so",level=1, lvar=3)  
  r3 <- brick(f[3],varname="talk",level=1, lvar=3)
  r4 <- brick(f[4], varname = "tos", lvar=4)
  for (i in 1:20){
    dissic <- r1[[i]]
    so <- r2[[i]]
    tos <- r4[[i]]
    talk <- r3[[i]]
    r <- brick(so,tos,phos,sili,talk,dissic)
    names(r)=c("salt","sst","phos","sili","talk","inorgC")
    sysdf <- as.data.frame(r,xy=TRUE)
    sysdf$pressure <- 0
    sysdf$Pressure <- 0
    sysdf$SST <- sysdf$sst
    cosysdf <- sysdf[,c(1,2,3,4,9,5,6,11,10,7,8)]
    write.csv(cosysdf,csvfnames[i],na="")
  }
}



###RCP26
#
parent.folder <- "E:/oa/co2sys_input/rcp26"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,full.names=T)
  r1 <- brick(f[1],varname="dissic",level=1, lvar=3)
  r2 <- brick(f[2],varname="so",level=1, lvar=3)  
  r3 <- brick(f[3],varname="talk",level=1, lvar=3)
  r4 <- brick(f[4], varname = "tos", lvar=4)
  #files are always in alphabetical order. run fiveyrmeanYd on dissic[1], fiveyrmeanMs on so[2], fiveyrmeanYt on talk[3], fiveyrmeanMt on tos[4]
  rasternewout <- paste((substr(f[1:4],1,nchar(f[1:4])-3)),".tif",sep="")
  fiveyrmeanYd(r1)
  fiveyrmeanMs(r2)
  fiveyrmeanYt(r3)
  fiveyrmeanMt(r4)
}


#create csv files
#columns: x, y, salinity, SST (c), pressure, totalP, totalSi, SST, Pressure, TA, TCO2
parent.folder <- "E:/oa/co2sys_input/rcp26"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
csvfnames <-paste("rcp26_",2000+(1:20)*5,"_co2sysin.csv",sep="")
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,pattern=".tif$",full.names=T)
  r1 <- brick(f[1],varname="dissic",level=1, lvar=3)
  r2 <- brick(f[2],varname="so",level=1, lvar=3)  
  r3 <- brick(f[3],varname="talk",level=1, lvar=3)
  r4 <- brick(f[4], varname = "tos", lvar=4)
  for (i in 1:20){
    dissic <- r1[[i]]
    so <- r2[[i]]
    tos <- r4[[i]]
    talk <- r3[[i]]
    r <- brick(so,tos,phos,sili,talk,dissic)
    names(r)=c("salt","sst","phos","sili","talk","inorgC")
    sysdf <- as.data.frame(r,xy=TRUE)
    sysdf$pressure <- 0
    sysdf$Pressure <- 0
    sysdf$SST <- sysdf$sst
    cosysdf <- sysdf[,c(1,2,3,4,9,5,6,11,10,7,8)]
    write.csv(cosysdf,csvfnames[i],na="")
  }
}



###historical
#

##function to calculate 5 year average from yearly data + units adjustment for dissic##
fiveyrmeanYd <- function(rasterxyz) {
  #unstack to groups of 5, average, then restack 
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  tenyr <- 10 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(rep(1:(ceiling((length(meanz_years)) / tenyr )), 
                    each = tenyr, length.out = length(meanz_years)))
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each 5 year period is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  #find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_tenyr <- brick(stacksout)
  names(meanz_tenyr) <- c(paste(1845+10*(1:16)))
  dissic_mean_10yr <- calc(meanz_tenyr, fun=function(x){x*1000/1.025}) #change from mol/m3 to umol/kgSW. x*1000/1.025
  dissic_mean_10yr <- rotate(dissic_mean_10yr)
  extent(dissic_mean_10yr) <- ext
  writeRaster(dissic_mean_10yr, filename=rasternewout[1], format="GTiff", overwrite=TRUE)
  
}

##so
fiveyrmeanMs <- function(rasterxyz) {
  #unstack to groups of 5 (e.g. 2006-2010,2011-2015, 2016-2020, 2020-2025), average, then restack 
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  tenyr <- 10*12 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(rep(1:(ceiling((length(meanz_years)) / tenyr )), 
                          each = tenyr, length.out = length(meanz_years)))
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each 5 year period is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  #find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_tenyr <- brick(stacksout)
  names(meanz_tenyr) <- c(paste(1845+10*(1:16)))
  meanz_tenyr <- rotate(meanz_tenyr)
  extent(meanz_tenyr) <- ext
  writeRaster(meanz_tenyr, filename=rasternewout[2], format="GTiff", overwrite=TRUE)
  
}

##talk
fiveyrmeanYt <- function(rasterxyz) {
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  tenyr <- 10 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(rep(1:(ceiling((length(meanz_years)) / tenyr )), 
                      each = tenyr, length.out = length(meanz_years)))
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each 5 year period is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  #find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_tenyr <- brick(stacksout)
  names(meanz_tenyr) <- c(paste(1845+10*(1:16)))
  talk_mean_10yr <- calc(meanz_tenyr, fun=function(x){x*1000/1.025}) #change from mol/m3 to umol/kgSW. x*1000/1.025
  talk_mean_10yr <- rotate(talk_mean_10yr)
  extent(talk_mean_10yr) <- ext
  writeRaster(talk_mean_10yr, filename=rasternewout[3], format="GTiff", overwrite=TRUE)
  
}

##function to calculate 5 year average from monthly data + units adjust K to C for tos ##
fiveyrmeanMt <- function(rasterxyz) {
  meanz_years <- unstack(rasterxyz) # now a list of rasters 
  tenyr <- 10*12 # desired layers per stack -  5 years (or 60 months)
  # uses names to group - here, names(mvlist) becomes c('1','1','1','2','2','2', etc) so that each year is in a group
  list <- c(rep(1:(ceiling((length(meanz_years)) / tenyr )), 
                            each = tenyr, length.out = length(meanz_years)))
  names(meanz_years) <- list
  # make list of rasters into a list of stacks. basically each decade is now one separate stack
  stacks <- lapply(unique(names(meanz_years)), function(y) {
    mean_years <- meanz_years[names(meanz_years) == y]
    stack(mean_years)
  })
  # find the mean for each 5 year period
  stacksout <- lapply(stacks, function(w) {
    g <- calc(w, function(x){mean(x)}, forceapply = TRUE)
  })
  meanz_tenyr <- stack(stacksout)
  names(meanz_tenyr) <- c(paste(1845+10*(1:16)))
  sst_mean_10yr <- calc(meanz_tenyr, fun=function(x){x-273.15})
  sst_mean_10yr <- rotate(sst_mean_10yr)
  extent(sst_mean_10yr) <- ext
  writeRaster(sst_mean_10yr, filename=rasternewout[4], format="GTiff", overwrite=TRUE)
  
}

#fix models so all go from 1850-2005
parent.folder <- "E:/oa/co2sys_input/historical/round2"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
e <- raster(nrow=180,ncol=360,ext=extent(m1),crs=crs(m1))
e[]<-NA
e10<-stack(replicate(10,e))
rd1 <- brick("E:/oa/co2sys_input/historical/round2/HadGEM2-CC/dissic_Oyr_HadGEM2-CC_historical_r1i1p1_1860-2005_grid.nc")
rt1 <-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-CC/talk_Oyr_HadGEM2-CC_historical_r1i1p1_1860-2005_grid.nc")
rd2<-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-ES/dissic_Oyr_HadGEM2-ES_historical_r1i1p1_1860-2005_grid.nc")
rt2<-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-ES/talk_Oyr_HadGEM2-ES_historical_r1i1p1_1860-2005_grid.nc")
rd1s<-stack(e10,rd1)
rt1s<-stack(e10,rt1)
rd2s<-stack(e10,rd2)
rt2s<-stack(e10,rt2)
writeRaster(rd1s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-CC/dissic_Oyr_HadGEM2-CC_historical_r1i1p1_1860-2005_grid_e.nc")
writeRaster(rt1s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-CC/talk_Oyr_HadGEM2-CC_historical_r1i1p1_1860-2005_grid_e.nc")
writeRaster(rd2s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-ES/dissic_Oyr_HadGEM2-ES_historical_r1i1p1_1860-2005_grid_e.nc")
writeRaster(rt2s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-ES/talk_Oyr_HadGEM2-ES_historical_r1i1p1_1860-2005_grid_e.nc")
e120<-stack(replicate(120,e))
rs1<-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-CC/so_Omon_HadGEM2-CC_historical_r1i1p1_185912-200511_grid.nc")
rto1<-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-CC/tos_Omon_HadGEM2-CC_historical_r1i1p1_185912-200511_grid.nc")
rs2<-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-ES/so_Omon_HadGEM2-ES_historical_r1i1p1_185912-200512_grid.nc")
rto2<-brick("E:/oa/co2sys_input/historical/round2/HadGEM2-ES/tos_Omon_HadGEM2-ES_historical_r1i1p1_185912-200511_grid.nc")
rs1s<-stack(e120,rs1)
rto1s<-stack(e120,rto1)
rs2s<-stack(e120,rs2)
rto2s<-stack(e120,rto2)
writeRaster(rs1s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-CC/so_Omon_HadGEM2-CC_historical_r1i1p1_185912-200511_grid_e.nc")
writeRaster(rto1s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-CC/tos_Omon_HadGEM2-CC_historical_r1i1p1_185912-200511_grid_e.nc")
writeRaster(rs2s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-ES/so_Omon_HadGEM2-ES_historical_r1i1p1_185912-200512_grid_e.nc")
writeRaster(rto2s,filename="E:/oa/co2sys_input/historical/round2/HadGEM2-ES/tos_Omon_HadGEM2-ES_historical_r1i1p1_185912-200511_grid_e.nc")


parent.folder <- "E:/oa/co2sys_input/historical/round2"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,pattern=".nc$",full.names=T)
  r1 <- brick(f[1],level=1, lvar=3) #dissic
  r2 <- brick(f[2],level=1, lvar=3)  #so
  r3 <- brick(f[3],level=1, lvar=3) #talk
  r4 <- brick(f[4],lvar=4) #tos
  #files are always in alphabetical order. run fiveyrmeanYd on dissic[1], fiveyrmeanMs on so[2], fiveyrmeanYt on talk[3], fiveyrmeanMt on tos[4]
  rasternewout <- paste((substr(f[1:4],1,nchar(f[1:4])-3)),".tif",sep="")
  fiveyrmeanYd(r1)
  fiveyrmeanMs(r2)
  fiveyrmeanYt(r3)
  fiveyrmeanMt(r4)
}


#create csv files
#columns: x, y, salinity, SST (c), pressure, totalP, totalSi, SST, Pressure, TA, TCO2
parent.folder <- "E:/oa/co2sys_input/historical/round2"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
csvfnames <-paste("hist_",1845+10*(1:16),"_co2sysin.csv",sep="")
for(i in sub.folders) {
  setwd(i)
  f<-list.files(i,pattern=".tif$",full.names=T)
  r1 <- brick(f[1],level=1, lvar=3) #dissic
  r2 <- brick(f[2],level=1, lvar=3)  #so
  r3 <- brick(f[3],level=1, lvar=3) #talk
  r4 <- brick(f[4],lvar=4) #tos
  for (i in 1:16){
    dissic <- r1[[i]]
    so <- r2[[i]]
    tos <- r4[[i]]
    talk <- r3[[i]]
    r <- brick(so,tos,phos,sili,talk,dissic)
    names(r)=c("salt","sst","phos","sili","talk","inorgC")
    sysdf <- as.data.frame(r,xy=TRUE)
    sysdf$pressure <- 0
    sysdf$Pressure <- 0
    sysdf$SST <- sysdf$sst
    cosysdf <- sysdf[,c(1,2,3,4,9,5,6,11,10,7,8)]
    write.csv(cosysdf,csvfnames[i],na="")
  }
}




##empirical

setwd("E:/oa/co2sys_input/empirical")
alk <- raster("Alkalinity.nc")
talk<-flip(alk,direction='x')
talk<-t(talk)
talk<-flip(talk,direction='x')
#talk<-t(talk)
tc <- raster("TCO2.nc")
tc<-flip(tc,direction='x')
tc<-t(tc)
tc<-flip(tc,direction='x')
#temp <- brick("D:/SST/OISST/oisstmasked.nc")
#temp1<- dropLayer(temp,c(1:276))
#temp<-dropLayer(temp1,c(13:174))
#SSTemp<-calc(temp,mean)
#SST <- rotate(SSTemp)
tempWOA <- raster("E:/oa/co2sys_input/empirical/WOA09/temperature_annual_1deg.nc",varname="t_an")
tempW <- rotate(tempWOA)
sal <- raster("E:/oa/co2sys_input/empirical/WOA09/salinity_annual_1deg.nc",varname="s_an")
salt <- rotate(sal)
phos <-raster("E:/oa/co2sys_input/empirical/WOA09/phosphate_annual_1deg.nc",varname="p_an")
phosp <- rotate(phos)
sili <-raster("E:/oa/co2sys_input/empirical/WOA09/silicate_annual_1deg.nc",varname="i_an")
silic <- rotate(sili)

r <- brick(salt,tempW,phosp,silic,talk,tc)
names(r)=c("salt","sst","phos","sili","talk","inorgC")
sysdf <- as.data.frame(r,xy=TRUE)
sysdf$pressure <- 0
sysdf$Pressure <- 0
sysdf$SST <- sysdf$sst
cosysdf <- sysdf[,c(1,2,3,4,9,5,6,11,10,7,8)]
write.csv(cosysdf,"empirical_2005_co2sysin.csv",na="")



#try empirical with GLODAP 2016 talk and tco2
setwd("E:/oa/co2sys_input/empirical")
alk <- raster("D:/omega_arag/co2sys_input/co2sys vars/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TAlk.nc")
alk1 <- crop(alk,extent(180,380,-90,90))
alk2 <- crop(alk,extent(20,180,-90,90))
extent(alk1) <- c(-180,20,-90,90)
talk <- merge(alk1,alk2)
tc <- raster("D:/omega_arag/co2sys_input/co2sys vars/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc")
tc1 <- crop(tc,extent(180,380,-90,90))
tc2 <- crop(tc,extent(20,180,-90,90))
extent(tc1) <- c(-180,20,-90,90)
tco <- merge(tc1,tc2)
tempWOA <- raster("E:/oa/co2sys_input/empirical/WOA09/temperature_annual_1deg.nc",varname="t_an")
tempW <- rotate(tempWOA)
sal <- raster("E:/oa/co2sys_input/empirical/WOA09/salinity_annual_1deg.nc",varname="s_an")
salt <- rotate(sal)
phos <-raster("E:/oa/co2sys_input/empirical/WOA09/phosphate_annual_1deg.nc",varname="p_an")
phosp <- rotate(phos)
sili <-raster("E:/oa/co2sys_input/empirical/WOA09/silicate_annual_1deg.nc",varname="i_an")
silic <- rotate(sili)

r <- brick(salt,tempW,phosp,silic,talk,tco)
names(r)=c("salt","sst","phos","sili","talk","inorgC")
sysdf <- as.data.frame(r,xy=TRUE)
sysdf$pressure <- 0
sysdf$Pressure <- 0
sysdf$SST <- sysdf$sst
cosysdf <- sysdf[,c(1,2,3,4,9,5,6,11,10,7,8)]
write.csv(cosysdf,"empirical_2005new_co2sysin.csv",na="")



################
################
################
#process through CO2Sys_v2.1.xls
################
################
################

#turn xy csv into rasters
#combine rasters into brick
#make model mean

###
#rcp85
parent.folder <- "E:/oa/co2sys_input/rcp85"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
brickname <- paste((substr(sub.folders[1:length(sub.folders)],1,nchar(sub.folders[1:length(sub.folders)]))),"_OAbrick.tif",sep="")
for(i in 1:length(sub.folders)) {
  setwd(sub.folders[i])
  f<-list.files(sub.folders[i],pattern="co2sysout.csv$",full.names=T)
  rastername <- paste((substr(f[1:length(f)],1,nchar(f[1:length(f)])-13)),"OA.tif",sep="")
  for (j in 1:length(f)){
    d <- read.csv(f[j],colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
    r <- rasterFromXYZ(d,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
    writeRaster(r, filename = rastername[j], format="GTiff", overwrite=TRUE)
  }
  l<-list.files(sub.folders[i],pattern="OA.tif$",full.names=T)
  for(k in 1:length(l)) {
    a <- paste("r", k, sep = "")
    r <- raster(l[k])
    assign(a,r)
  }
  br <- brick(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20)
  writeRaster(br, filename = brickname[i], format="GTiff", overwrite=TRUE)
}

m <- list.files(parent.folder, pattern="OAbrick.tif",full.names=T)
for(i in 1:length(m)){
  a <- paste("m", i, sep = "")
  r <- brick(m[i])
  assign(a,r)
}
#mean&median calculation
modelm <- overlay(m1,m2,m3,m4,m5,m6,m7,fun=function(x){ mean(x,na.rm=T)})
writeRaster(modelm, filename = "E:/oa/co2sys_input/rcp85/OArcp85_modelmean.tif", format="GTiff", overwrite=TRUE)
modelmed <- overlay(m1,m2,m3,m4,m5,m6,m7,fun=function(x){ median(x,na.rm=T)})
writeRaster(modelmed, filename = "E:/oa/co2sys_input/rcp85/OArcp85_modelmedian.tif", format="GTiff", overwrite=TRUE)

modelmed <- brick("E:/oa/co2sys_input/rcp85/OArcp85_modelmedian.tif")
for(i in 1:nlayers(modelmed)){
  modelmedf <- focal(modelmed[[i]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                     filename=paste0("E:/oa/co2sys_input/rcp85/OArcp85_modelmedianf",2000+i*5,".tif"))
}
#make all focal median files a single brick
medffiles <- list.files(path="E:/oa/co2sys_input/rcp85/",pattern='modelmedianf',full.names=T)
medfr <- lapply(medffiles, raster::raster)
medfb <- brick(medfr)
writeRaster(medfb,filename="E:/oa/co2sys_input/rcp85/OArcp85_modelmedianf.tif",overwrite=T)



###
#rcp45

parent.folder <- "E:/oa/co2sys_input/rcp45"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
brickname <- paste((substr(sub.folders[1:length(sub.folders)],1,nchar(sub.folders[1:length(sub.folders)]))),"_OAbrick.tif",sep="")
for(i in 1:length(sub.folders)) {
  setwd(sub.folders[i])
  f<-list.files(sub.folders[i],pattern="co2sysout.csv$",full.names=T)
  rastername <- paste((substr(f[1:length(f)],1,nchar(f[1:length(f)])-13)),"OA.tif",sep="")
  for (j in 1:length(f)){
    d <- read.csv(f[j],colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
    r <- rasterFromXYZ(d,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
    writeRaster(r, filename = rastername[j], format="GTiff", overwrite=TRUE)
  }
  l<-list.files(sub.folders[i],pattern="OA.tif$",full.names=T)
  for(k in 1:length(l)) {
    a <- paste("r", k, sep = "")
    r <- raster(l[k])
    assign(a,r)
  }
  br <- brick(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20)
  writeRaster(br, filename = brickname[i], format="GTiff", overwrite=TRUE)
}

m <- list.files(parent.folder, pattern="OAbrick.tif",full.names=T)
for(i in 1:length(m)){
  a <- paste("m", i, sep = "")
  r <- brick(m[i])
  assign(a,r)
}
#mean&median calculation
modelm <- overlay(m1,m2,m3,m4,m5,m6,m7,fun=function(x){ mean(x,na.rm=T)})
writeRaster(modelm, filename = "E:/oa/co2sys_input/rcp45/OArcp45_modelmean.tif", format="GTiff", overwrite=TRUE)
modelmed <- overlay(m1,m2,m3,m4,m5,m6,m7,fun=function(x){ median(x,na.rm=T)})
writeRaster(modelmed, filename = "E:/oa/co2sys_input/rcp45/OArcp45_modelmedian.tif", format="GTiff", overwrite=TRUE)

modelmed <- brick("E:/oa/co2sys_input/rcp45/OArcp45_modelmedian.tif")
for(i in 1:nlayers(modelmed)){
  modelmedf <- focal(modelmed[[i]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                     filename=paste0("E:/oa/co2sys_input/rcp45/OArcp45_modelmedianf",2000+i*5,".tif"))
}
#make all focal median files a single brick
medffiles <- list.files(path="E:/oa/co2sys_input/rcp45/",pattern='modelmedianf',full.names=T)
medfr <- lapply(medffiles, raster::raster)
medfb <- brick(medfr)
writeRaster(medfb,filename="E:/oa/co2sys_input/rcp45/OArcp45_modelmedianf.tif",overwrite=T)





###
#rcp26
parent.folder <- "E:/oa/co2sys_input/rcp26"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
brickname <- paste((substr(sub.folders[1:length(sub.folders)],1,nchar(sub.folders[1:length(sub.folders)]))),"_OAbrick.tif",sep="")
for(i in 1:length(sub.folders)) {
  setwd(sub.folders[i])
  f<-list.files(sub.folders[i],pattern="co2sysout.csv$",full.names=T)
  rastername <- paste((substr(f[1:length(f)],1,nchar(f[1:length(f)])-13)),"OA.tif",sep="")
  for (j in 1:length(f)){
    d <- read.csv(f[j],colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
    r <- rasterFromXYZ(d,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
    writeRaster(r, filename = rastername[j], format="GTiff", overwrite=TRUE)
  }
  l<-list.files(sub.folders[i],pattern="OA.tif$",full.names=T)
  for(k in 1:length(l)) {
    a <- paste("r", k, sep = "")
    r <- raster(l[k])
    assign(a,r)
  }
  br <- brick(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20)
  writeRaster(br, filename = brickname[i], format="GTiff", overwrite=TRUE)
}

m <- list.files(parent.folder, pattern="OAbrick.tif",full.names=T)
for(i in 1:length(m)){
  a <- paste("m", i, sep = "")
  r <- brick(m[i])
  assign(a,r)
}
#mean&median calculation
modelm <- overlay(m1,m2,m3,m4,m5,m6,fun=function(x){ mean(x,na.rm=T)}, filename = "E:/oa/co2sys_input/rcp26/OArcp26_modelmean.tif")
modelmed <- overlay(m1,m2,m3,m4,m5,m6,fun=function(x){ median(x,na.rm=T)}, filename = "E:/oa/co2sys_input/rcp26/OArcp26_modelmedian.tif")
modelmed <- brick("E:/oa/co2sys_input/rcp26/OArcp26_modelmedian.tif")
for(i in 1:nlayers(modelmed)){
  modelmedf <- focal(modelmed[[i]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                     filename=paste0("E:/oa/co2sys_input/rcp26/OArcp26_modelmedianf",2000+i*5,".tif"))
}
#make all focal median files a single brick
medffiles <- list.files(path="E:/oa/co2sys_input/rcp26/",pattern='modelmedianf',full.names=T)
medfr <- lapply(medffiles, raster::raster)
medfb <- brick(medfr)
writeRaster(medfb,filename="E:/oa/co2sys_input/rcp26/OArcp26_modelmedianf.tif",overwrite=T)




###
#historical

#calculate mean and median of hist models
parent.folder <- "E:/oa/co2sys_input/historical"
sub.folders <- list.dirs(parent.folder,recursive=T)[-1]
brickname <- paste((substr(sub.folders[1:length(sub.folders)],1,nchar(sub.folders[1:length(sub.folders)]))),"_OAbrick.tif",sep="")
for(i in 2:10) {
  setwd(sub.folders[i])
  f<-list.files(sub.folders[i],pattern="co2sysout.csv$",full.names=T)
  rastername <- paste((substr(f[1:length(f)],1,nchar(f[1:length(f)])-13)),"OA.tif",sep="")
  for (j in 1:length(f)){
    d <- read.csv(f[j],colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
    r <- rasterFromXYZ(d,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
    writeRaster(r, filename = rastername[j], format="GTiff", overwrite=TRUE)
  }
  l<-list.files(sub.folders[i],pattern="OA.tif$",full.names=T)
  for(k in 1:length(l)) {
    a <- paste("r", k, sep = "")
    r <- raster(l[k])
    assign(a,r)
  }
  br <- brick(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16)
  writeRaster(br, filename = brickname[i], format="GTiff", overwrite=TRUE)
}

m <- list.files(parent.folder, pattern="OAbrick.tif",full.names=T)
for(i in 1:length(m)){
  a <- paste("m", i, sep = "")
  r <- brick(m[i])
  assign(a,r)
}
#add 10 empty rasters to m1 so it matches the rest with 1855-2005
e <- raster(nrow=180,ncol=360,ext=extent(m1),crs=crs(m1))
e[]<-NA
e10<-stack(e,e,e,e,e,e,e,e,e,e)
m1e<-stack(e10,m1)
#get rid of 0 values along coastline for CNRM model (m3)
m3[m3==0] <- NA
#mean&median calculation
modelm <- overlay(m1e,m2,m3,m4,m5,m6,m7,m8,m9,m10,fun=function(x){ mean(x,na.rm=T)}, filename = "E:/oa/co2sys_input/historical/OAhist_modelmean.tif")
modelmed <-overlay(m1e,m2,m3,m4,m5,m6,m7,m8,m9,m10,fun=function(x){ median(x,na.rm=T)},filename = "E:/oa/co2sys_input/historical/OAhist_modelmedian.tif")
#apply focal calculation to give values to the NA coastal areas
for(i in 1:nlayers(modelmed)){
  modelmedf <- focal(modelmed[[i]],w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
                     filename=paste0("E:/oa/co2sys_input/historical/OAhist_modelmedianf",1845+i*10,".tif"))
}
#make all focal median files a single brick
medffiles <- list.files(path="E:/oa/co2sys_input/historical",pattern='modelmedianf',full.names=T)
medfr <- lapply(medffiles, raster::raster)
medfb <- brick(medfr)
writeRaster(medfb,filename="E:/oa/co2sys_input/historical/OAhist_modelmedianf.tif",overwrite=T)


####
##empirical to raster
###

empdf <- read.csv("E:/oa/co2sys_input/empirical/empirical_2005_co2sysout.csv",colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),
                  header=T,row.names=NULL)
emprr <- rasterFromXYZ(empdf,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
writeRaster(emprr, filename = "E:/oa/empirical/OAempirical.tif", format="GTiff", overwrite=TRUE)
emprr[emprr==0] <- NA
emprf <- focal(emprr,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
               filename="E:/oa/empirical/OAempiricalf.tif")

#test with empirical GLODAP2016
empdf <- read.csv("E:/oa/co2sys_input/empirical/empirical_2005new_co2sysout.csv",colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),
                  header=T,row.names=NULL)
emprr <- rasterFromXYZ(empdf,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 ")
writeRaster(emprr, filename = "E:/oa/empirical/OAempiricalnew.tif", format="GTiff", overwrite=TRUE)
emprr[emprr==0] <- NA
emprf <- focal(emprr,w=matrix(1,nrow=3,ncol=3),fun=mean,na.rm=T,NAonly=T,
               filename="E:/oa/empirical/OAempiricalfnew.tif")


###
#empirical & hist taylor diagrams

#calculate accuracy: Taylor diagrams
#get 2005 values for each model to compare to empirical values that were taken in ~2005
empirical <- read.csv("E:/oa/co2sys_input/empirical/empirical_2005new_co2sysout.csv")
emp2005 <- subset(empirical,select=c("x","y","WAr.out"))
names(emp2005) <- c("x","y","empirical")

md1<-read.csv("E:/oa/co2sys_input/historical/CanESM2/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md1) <-c("x","y","CanESM2")
md2<-read.csv("E:/oa/co2sys_input/historical/CMCC-CESM/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md2) <-c("x","y","CMCC-CESM")
md3<-read.csv("E:/oa/co2sys_input/historical/MIROC-ESM/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md3) <-c("x","y","GFDL-ESM2M")
md4<-read.csv("E:/oa/co2sys_input/historical/MIROC-ESM-CHEM/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md4) <-c("x","y","MIROC-ESM-CHEM")
md5<-read.csv("E:/oa/co2sys_input/historical/CNRM-CM5/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md5) <-c("x","y","CNRM-CM5")
md6<-read.csv("E:/oa/co2sys_input/historical/HadGEM2-CC/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md6) <-c("x","y","HadGEM2-CC")
md7<-read.csv("E:/oa/co2sys_input/historical/MPI-ESM-LR/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md7) <-c("x","y","MPI-ESM-LR")
md8<-read.csv("E:/oa/co2sys_input/historical/MPI-ESM-MR/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md8) <-c("x","y","MPI-ESM-MR")
md9<-read.csv("E:/oa/co2sys_input/historical/HadGEM2-ES/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md9) <-c("x","y","HadGEM2-ES")
md10<-read.csv("E:/oa/co2sys_input/historical/NorESM1-ME/hist_2005_co2sysout.csv",
              colClasses=c("NULL",rep("numeric",2),rep("NULL",10),"numeric"),header=T,row.names=NULL)
names(md10) <-c("x","y","NorESM1-ME")

modelm <- brick("E:/oa/co2sys_input/historical/OAhist_modelmean.tif")
modelmed <-brick("E:/oa/co2sys_input/historical/OAhist_modelmedian.tif")

mm <- raster(modelm,layer=16)
mmp <-rasterToPoints(mm)
mm05 <- as.data.frame(mmp)
names(mm05) <-c("x","y","modelmean")

mmd <- raster(modelmed,layer=16)
mmdp <-rasterToPoints(mmd)
mmd05 <- as.data.frame(mmdp)
names(mmd05) <-c("x","y","modelmedian")

alldata <- join_all(list(emp2005,md1,md2,md3,md4,md5,md6,md7,md8,md9,md10,mm05,mmd05),by=c("x","y"))

oldpar<-taylor.diagram(alldata$empirical,alldata$empirical,ref.sd=T,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="blue",gamma.col="blue",
                       main="Taylor Diagram - OA historical")
oldpar<-taylor.diagram(alldata$empirical,alldata$CanESM2,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`CMCC-CESM`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`CNRM-CM5`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`HadGEM2-CC`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`HadGEM2-ES`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`MIROC-ESM`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`MIROC-ESM-CHEM`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`MPI-ESM-LR`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`MPI-ESM-MR`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$`NorESM1-ME`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldata$empirical,alldata$modelmean,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,col="red",pch=19)
oldpar<-taylor.diagram(alldata$empirical,alldata$modelmedian,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,col="green",pch=19)



#taylor diagram of just coral sites
  emp2005 <- read.csv(file="E:/oa/empirical/emp_2005new_suit.csv",row.names=NULL,header=T,colClasses=c(rep("NULL",2),"numeric","numeric","numeric","NULL"))
    names(emp2005) <- c("iscoral","ID","empirical")

  
  ##get mean values for each model, model mean, model median for same test period - 2005
  ##extract only values at coral points
  
  parent.folder <- "E:/oa/co2sys_input/historical"
  
  
  m <- list.files(parent.folder, pattern="_OAbrick.tif",full.names=T)
  m1 <- raster(m[1],band=6) #separate because this one only has 6 layers
  for(i in 2:length(m)){
    a <- paste("m", i, sep = "")
    r <- raster(m[i],band=16)
    assign(a,r)
  }
  
  rasterList <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
  for(i in 1:length(rasterList)){
    r<-rasterList[[i]]
    reefext <- rgis::fast_extract(sf=reefs,ras=r,funct=function(x){max(x,na.rm=T)},parallel=T)
    dfext <- as.data.frame(reefext)
    rdf <- dfext[,c(4,5)]
    names(rdf) <- c("ID",paste("model",i))
    a <- paste("mod", i, sep = "")
    assign(a,rdf)
  }
  
  modelm <- raster("E:/oa/co2sys_input/historical/OAhist_modelmean.tif",band=16)
  reefext <- rgis::fast_extract(sf=reefs,ras=modelm,funct=function(x){max(x,na.rm=T)},parallel=T)
  dfext <- as.data.frame(reefext)
  mndf <- dfext[,c(4,5)]
  names(mndf) <- c("ID","meanmod")
  
  modelmed <- raster("E:/oa/co2sys_input/historical/OAhist_modelmedian.tif",band=16)
    reefext <- rgis::fast_extract(sf=reefs,ras=modelmed,funct=function(x){max(x,na.rm=T)},parallel=T)
  dfext <- as.data.frame(reefext)
  medf <- dfext[,c(4,5)]
  names(medf) <- c("ID","medianmod")
  


alldata <- join_all(list(emp2005,mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mndf,medf),by=c("ID"))
alldf<-alldata[!(alldata$iscoral==0),]
alldf<- alldf[!(alldf$empirical==0),]


oldpar<-taylor.diagram(alldf$empirical,alldf$empirical,ref.sd=T,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,cex.axis=T,col="blue",gamma.col="blue",
                       main="Taylor Diagram - OA historical, coral sites")
oldpar<-taylor.diagram(alldf$empirical,alldf$'model 1',add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,cex.axis=T,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$'model 2',add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,cex.axis=T,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 3`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,cex.axis=T,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 4`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 5`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 6`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 7`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 8`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 9`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$`model 10`,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,pch=19,col="black")
oldpar<-taylor.diagram(alldf$empirical,alldf$meanmod,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,col="red",pch=19)
oldpar<-taylor.diagram(alldf$empirical,alldf$medianmod,add=TRUE,normalize=TRUE,sd.arcs=TRUE,pcex=1.5,col="green",pch=19)


