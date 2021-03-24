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


###RCP26

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp26")


reefs <- st_read("E:/currentcoral/reefs.shp")
m <- list.files("F:/dhw/rcp26/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
rcp26suitcsvout <- paste("rcp26_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
count<-data.frame()
foreach::foreach(i=1:20, .packages=c("raster")) %dopar% {
#for (i in 1:20){
  x <- paste0("dhw_modelmedian_bc_5y",2000+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  write.csv(dfi,rcp26suitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start
#22 min


#get count of suitability
csvfiles <- list.files("E:/dhw/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 




###RCP45

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp45")


reefs <- st_read("E:/currentcoral/reefs.shp")
m <- list.files("F:/dhw/rcp45/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
rcp45suitcsvout <- paste("rcp45_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:20, .packages=c("raster","plyr","tibble")) %dopar% {
  x <- paste0("dhw_modelmedian_bc_5y",2000+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit05[dfi$dhw <=8] <- 1
  dfi$suit05[dfi$dhw > 8] <- 0 #suitable number of dhw is <8
  write.csv(dfi,rcp45suitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 



###RCP85

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/rcp85")


reefs <- st_read("E:/currentcoral/reefs.shp")
m <- list.files("F:/dhw/rcp85/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
rcp85suitcsvout <- paste("rcp85_",2000+(1:20)*5,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:20, .packages=c("raster","plyr","tibble")) %dopar% {
#for (i in 1:20){
  x <- paste0("dhw_modelmedian_bc_5y",2000+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit05[dfi$dhw <=8] <- 1
  dfi$suit05[dfi$dhw > 8] <- 0 #suitable number of dhw is <8
  write.csv(dfi,rcp85suitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start

#get count of suitability
csvfiles <- list.files("E:/dhw/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
} 



###historic oisst

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/historic/empirical")


reefs <- st_read("E:/currentcoral/reefs.shp")
m <- list.files("F:/dhw/emp/", pattern=".tif",full.names=T)
dhw <- raster()
for(i in 1:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
histsuitcsvout <- paste("hist_",1980+(1:7)*5,"_suit.csv",sep="")



start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

count<-data.frame()

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:7, .packages=c("raster","plyr")) %dopar% {
#for (i in 1:7){
  x <- paste0("dhw_emp_oisst_5y",1980+i*5)
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars=dhw_coral)
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/historic/empirical",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+1980),.before=T)
  count <- rbind(y,count)
} 




#mean DHW under all scenarios&empirical

dhwdf<-data.frame()
csvfiles <- list.files("E:/dhw/historic/empirical",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
for(i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  m <- mean(a$dhw,na.rm=T)
  df <- data.frame(year = c(paste0(i*5+1980)),coralmean = c(m),scenario="empirical")
  dhwdf <- rbind(df,dhwdf)
}

csvfiles <- list.files("E:/dhw/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
for(i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  m <- mean(a$dhw,na.rm=T)
  df <- data.frame(year = c(paste0(i*10+1845)),coralmean = c(m),scenario="historical")
  dhwdf <- rbind(df,dhwdf)
}

csvfiles <- list.files("E:/dhw/rcp26",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
for(i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  m <- mean(a$dhw,na.rm=T)
  df <- data.frame(year = c(paste0(i*5+2000)),coralmean = c(m),scenario="rcp26")
  dhwdf <- rbind(df,dhwdf)
}
write.csv(dhwdf,"DHW.csv",row.names=F)

csvfiles <- list.files("E:/dhw/rcp45",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
for(i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  m <- mean(a$dhw,na.rm=T)
  df <- data.frame(year = c(paste0(i*5+2000)),coralmean = c(m),scenario="rcp45")
  dhwdf <- rbind(df,dhwdf)
}
write.csv(dhwdf,"DHW.csv",row.names=F)

csvfiles <- list.files("E:/dhw/rcp85",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
for(i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  m <- mean(a$dhw,na.rm=T)
  df <- data.frame(year = c(paste0(i*5+2000)),coralmean = c(m),scenario="rcp85")
  dhwdf <- rbind(df,dhwdf)
}
write.csv(dhwdf,"DHW.csv",row.names=F)

dhwdf$year <- as.numeric(dhwdf$year)

ggp <- ggplot(dhwdf,aes(year,coralmean, col=scenario))+
  geom_line()
ggp


dhwsuit <- read.csv(file="dhwsuit.csv",row.names=NULL,header=T)
dhwsuit$year <- as.numeric(dhwsuit$year)
dhwsuit$percent.suitable <- as.numeric(dhwsuit$percent.suitable)

p<- ggplot(dhwsuit, aes(x=year, y=suitable, group=scenario, color=scenario)) + 
  geom_line(size=1.5) #+
#geom_point() 


p+labs(title="Percentage of Suitable Sites under DHW")+
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.title = element_blank()) +
  xlab("year")+ylab("% of reef sites")

#historic model

####
#extract reef suit vals
#find dhw for each point
#classify as suitable or not - save as .csv
###

setwd("E:/dhw/historic")


reefs <- st_read("E:/currentcoral/reefs.shp")
m <- list.files("F:/dhw/hist/", pattern="median",full.names=T)
dhw <- raster()
for(i in 2:length(m)){
  r <- raster(m[i])
  dhw <- stack(dhw,r)
}
histsuitcsvout <- paste("hist_",1845+(1:16)*10,"_suit.csv",sep="")


start <- Sys.time()
cores<- detectCores()-1
cl <- makeCluster(cores, output="") 
registerDoParallel(cl)

reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
foreach::foreach(i=1:16, .packages=c("raster","plyr")) %dopar% {
  #for (i in 1:16){
  x <- paste0("dhw_modelmedian_bc_10y",1845+i*10,sep="")
  dfi <- subset(dfext,select=c("pointid","grid_code", "iscoral","ID",x))
  names(dfi) <- c("pointid","grid_code", "iscoral","ID","dhw")
  dfi$suit <- with(dfi, ifelse(dhw <= 8,1, 0)) #suitable number of dhw is <8
  write.csv(dfi,histsuitcsvout[i], row.names=FALSE)
  
  dfi$dhw_coral <- dfi$suit * dfi$iscoral
  c <- count(dfi, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*5+2000),.before=T)
  count <- rbind(y,count)
}

stopCluster(cl)
finish <- Sys.time()
finish-start


#get count of suitability
csvfiles <- list.files("E:/dhw/historic",pattern="suit.csv$",full.names=T)
ncsvfiles <- length(csvfiles)
count1<-data.frame()
for (i in 1:ncsvfiles){
  a <- read.csv(file=csvfiles[i],row.names=NULL,header=T)
  a$dhw_coral <- a$suit * a$iscoral
  c <- count(a, vars="dhw_coral")
  y <- add_column(c,year=paste0(i*10+1845),.before=T)
  count1 <- rbind(y,count1)
} 








#check just raw dhwhist 2005 suitability
dhwall <- brick("F:/dhwhist/ssthist_DHW_annualn.tif")
dhw <- dhwall[[156]]
reefext <- rgis::fast_extract(sf=reefs,ras=dhw,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
dfext$suit <- with(dfext, ifelse(dhw_ssthist_DHW_annualn.156 <= 8,1, 0))
dfext$dhw_coral <- dfext$suit * dfext$iscoral
count(dfext, vars="dhw_coral")

#check just raw dhwhist 2005 suitability
dhwall85 <- brick("E:/dhw/rcp85/sst85_DHW_annual.tif")
dhw85 <- dhwall85[[1]]
reefext <- rgis::fast_extract(sf=reefs,ras=dhw85,funct=function(x){max(x,na.rm=T)},parallel=T)
dfext <- as.data.frame(reefext)
dfext$suit <- with(dfext, ifelse(dhw85_sst85_DHW_annual.1 <= 8,1, 0))
dfext$dhw_coral <- dfext$suit * dfext$iscoral
count(dfext, vars="dhw_coral")
