library(ggplot2)
library(reshape2)
library(jcolors)
library(ggpubr)
library(RColorBrewer)
library(plyr)
library(caret)
library(pROC)



setwd("E:/figures/")

#### overall suitability line plot
suit <- read.csv(file="overallsuit.csv",header=T,sep=",")


p<- ggplot(suit, aes(x=year, y=suitable, group=RCP, color=RCP)) + 
  geom_line(size=1.5) #+
  #geom_point() 


p+labs(title="Projected Percentage of Suitable Sites")+
  theme_minimal() +
  #scale_color_jcolors(palette = "pal2") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.title = element_blank()) +
  xlab("year")+ylab("% of reef sites")





#### suitability each variable
suitvar <- read.csv(file="suitvars.csv",header=T, sep=",")

lp <- ggplot(suitvar, aes(x=year, y=suitable, group=variable))+
  geom_line(aes(color=variable), size=2)

lp + facet_wrap(. ~ RCP, scales="free_x",nrow=1)+
  theme_light() +
  scale_color_brewer(palette='Set2')+
  scale_y_continuous(labels = scales::percent) +
  xlab("year")+ylab("% of reef sites")

  
###
# historic = read.csv(file="historic.csv",header=T, sep=",")
#   
# hist <- ggline(historic, x='year', y='suitable',group = 'variable')
# hist




#### suitability all in one plot
suit <- read.csv(file="suitability.csv",header=T, sep=",")
suit$variables <- factor(suit$variable, levels = c("Storms","Land Use","SST","Population","OA","DHW","Overall"))


lp <- ggplot(suit, aes(x=year, y=suitable, group=variables))+
  geom_line(aes(color=variables,alpha=variables), size=2)+
  scale_alpha_manual(values=c(0.25,0.25,0.25,0.25,0.25,0.25,1)) 

lp + facet_wrap(. ~ RCP, scales="free_x",nrow=1)+
  theme_light() +
  #scale_color_brewer(palette='Set2')+
  scale_color_manual(values=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#9898ab"))+
  scale_y_continuous(labels = scales::percent) +
    scale_x_continuous()+
  xlab("year")+ylab("% of reef sites")+
  scale_fill_discrete(breaks = rev(levels(suit$variables)))



#### suitability all in one plot
suit <- read.csv(file="suitability1.csv",header=T, sep=",")
suit$variables <- factor(suit$variable, levels = c("Storms","Land Use","SST","Population","OA","DHW","Overall",
                                                   "Storms - historic","Land Use - historic","SST - historic",
                                                   "Population - historic","OA - historic","DHW - historic","Overall - historic"))
clr <- c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3", "#A6D854","#FFD92F", "#000000",
         "#66C2A5","#FC8D62","#8DA0CB","#E78AC3", "#A6D854","#FFD92F", "#000000")


lp <- ggplot(suit, aes(x=year, y=suitable, group=variables))+
  geom_line(aes(color=variables,linetype=histproj),size=1.5,alpha=0.5)
  
  #geom_line(aes(color=variables,alpha=variables,linetype=variables,size=variables),size=2)+
  # scale_alpha_manual(values=c(0.3,0.3,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,0.3,0.3,1))+
  # scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid","solid",
  #                                "longdash","longdash","longdash","longdash","longdash","longdash","longdash"))+
  # scale_size_manual(values=c(2,2,2,2,2,2,2,1.3,1.3,1.3,1.3,1.3,1.3,1.3))



lp + facet_wrap(. ~ RCP, scales="free_x",nrow=1)+
  theme_light() +
  #scale_color_brewer(palette='Set2')+
  scale_color_manual(name="Variable",values=clr,breaks = c("Storms","Land Use","SST","Population","OA","DHW","Overall"))+
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks=seq(1980,2100,20))+
  xlab("Year")+ylab("Percentage of reef sites")+
  labs(linetype = " ")+
  ggtitle("Suitability of reef sites")


###
##confusion matrix - SST
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (1985,1995,2005)


#load historical models ("data")
ssthist85 <- read.csv(file="E:/sst/historic/hist_1985_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(ssthist85) <- c("iscoral","ID","suitmod")
ssthist95 <- read.csv(file="E:/sst/historic/hist_1995_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(ssthist95) <- c("iscoral","ID","suitmod")
ssthist05 <- read.csv(file="E:/sst/historic/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(ssthist05) <- c("iscoral","ID","suitmod")


#load empirical data ("reference")
sstemp85 <- read.csv(file="E:/sst/empirical/hist_1985_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(sstemp85) <- c("iscoral","ID","suitemp")
sstemp95 <- read.csv(file="E:/sst/empirical/hist_1995_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(sstemp95) <- c("iscoral","ID","suitemp")
sstemp05 <- read.csv(file="E:/sst/empirical/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(sstemp05) <- c("iscoral","ID","suitemp")

#merge historical and empirical by year with pointid. remove non-coral rows with filter
sst85 <- join_all(list(ssthist85,sstemp85),by=c("ID"))
sst85 <- sst85[!(sst85$iscoral==0),]
sst95 <- join_all(list(ssthist95,sstemp95),by=c("ID"))
sst95 <- sst95[!(sst95$iscoral==0),]
sst05 <- join_all(list(ssthist05,sstemp05),by=c("ID"))
sst05 <- sst05[!(sst05$iscoral==0),]

#merge all 3 dataframes together 
sst859505 <- do.call("rbind",list(sst85,sst95,sst05))


#confusion matrix
cmatrix_sst <- confusionMatrix(sst859505$suitmod,sst859505$suitemp)






###
##confusion matrix - DHW
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (1985,1995,2005)


#load historical models ("data")
dhwhist85 <- read.csv(file="E:/dhw/historic/hist_1985_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(dhwhist85) <- c("iscoral","ID","suitmod")
dhwhist95 <- read.csv(file="E:/dhw/historic/hist_1995_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(dhwhist95) <- c("iscoral","ID","suitmod")
dhwhist05 <- read.csv(file="E:/dhw/historic/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(dhwhist05) <- c("iscoral","ID","suitmod")


#load empirical data ("reference")
dhwemp85 <- read.csv(file="E:/dhw/historic/empirical/hist_1985_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(dhwemp85) <- c("iscoral","ID","suitemp")
dhwemp95 <- read.csv(file="E:/dhw/historic/empirical/hist_1995_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(dhwemp95) <- c("iscoral","ID","suitemp")
dhwemp05 <- read.csv(file="E:/dhw/historic/empirical/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(dhwemp05) <- c("iscoral","ID","suitemp")

#merge historical and empirical by year with pointid. remove non-coral rows with filter
dhw85 <- join_all(list(dhwhist85,dhwemp85),by=c("ID"))
dhw85 <- dhw85[!(dhw85$iscoral==0),]
dhw95 <- join_all(list(dhwhist95,dhwemp95),by=c("ID"))
dhw95 <- dhw95[!(dhw95$iscoral==0),]
dhw05 <- join_all(list(dhwhist05,dhwemp05),by=c("ID"))
dhw05 <- dhw05[!(dhw05$iscoral==0),]

#merge all 3 dataframes together 
dhw859505 <- do.call("rbind",list(dhw85,dhw95,dhw05))


#confusion matrix
cmatrix_dhw <- confusionMatrix(dhw859505$suitmod,dhw859505$suitemp)







###
##confusion matrix - oa
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (2005)

#load historical models ("data")
oahist05 <- read.csv(file="E:/oa/hist/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(oahist05) <- c("iscoral","ID","suitmod")

#load empirical data ("reference")
oaemp05 <- read.csv(file="E:/oa/empirical/emp_2005new_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(oaemp05) <- c("iscoral","ID","suitemp")

#merge historical and empirical by year with pointid. remove non-coral rows with filter. also remove NA sites (there are a lot of them ...78148 sites... 45% coral sites)
oa05 <- join_all(list(oahist05,oaemp05),by=c("ID"))
oa05 <- oa05[!(oa05$iscoral==0),]
oa05 <- oa05[!(is.na(oa05$suitemp)),]

#confusion matrix
cmatrix_oa <- confusionMatrix(oa05$suitmod,oa05$suitemp)





###
##confusion matrix - storms
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (1985,1995,2005)


#load historical models ("data")
stormhist85 <- read.csv(file="E:/storms/historic/hist_1985_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(stormhist85) <- c("iscoral","ID","suitmod")
stormhist95 <- read.csv(file="E:/storms/historic/hist_1995_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(stormhist95) <- c("iscoral","ID","suitmod")
stormhist05 <- read.csv(file="E:/storms/historic/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(stormhist05) <- c("iscoral","ID","suitmod")


#load empirical data ("reference")
stormemp85 <- read.csv(file="E:/storms/empirical/emp_1985_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(stormemp85) <- c("iscoral","ID","suitemp")
stormemp95 <- read.csv(file="E:/storms/empirical/emp_1995_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(stormemp95) <- c("iscoral","ID","suitemp")
stormemp05 <- read.csv(file="E:/storms/empirical/emp_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c(rep("NULL",2),"numeric","numeric","NULL","factor"))
names(stormemp05) <- c("iscoral","ID","suitemp")

#merge historical and empirical by year with pointid. remove non-coral rows with filter
storm85 <- join_all(list(stormhist85,stormemp85),by=c("ID"))
storm85 <- storm85[!(storm85$iscoral==0),]
storm95 <- join_all(list(stormhist95,stormemp95),by=c("ID"))
storm95 <- storm95[!(storm95$iscoral==0),]
storm05 <- join_all(list(stormhist05,stormemp05),by=c("ID"))
storm05 <- storm05[!(storm05$iscoral==0),]

#merge all 3 dataframes together 
storm859505 <- do.call("rbind",list(storm85,storm95,storm05))


#confusion matrix
cmatrix_storm <- confusionMatrix(storm859505$suitmod,storm859505$suitemp)





###
##confusion matrix - land
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (2005)

#load historical model ("data")
landhist05 <- read.csv(file="E:/land_hurtt/historic/hist_2005_suit.csv",row.names=NULL,header=T,
                       colClasses=c(rep("NULL",2),rep("numeric",2),rep("NULL",3),"factor",rep("NULL",2)))
names(landhist05) <- c("iscoral","ID","suitmod")

#load empirical data ("reference)
landemp05 <- read.csv(file="E:/land_hurtt/empirical/emp_2005_suit.csv",row.names=NULL,header=T,
                      colClasses=c(rep("NULL",2),rep("numeric",2),rep("NULL",2),"factor"))
names(landemp05) <- c("iscoral","ID","suitemp")

#merge historical and empirical by year with pointid.oa05 <- join_all(list(oahist05,oaemp05),by=c("ID"))
land05 <- join_all(list(landhist05,landemp05),by=c("ID"))
land05 <- land05[!(land05$iscoral==0),]

#confusion matrix
cmatrix_land <- confusionMatrix(land05$suitmod,land05$suitemp)







###
##confusion matrix - population
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (2005)

#load historical model ("data")
pophist05 <- read.csv(file="E:/pop_hr/hist_hyde/histhyde_2005_suit.csv",row.names=NULL,header=T,
                       colClasses=c(rep("NULL",2),rep("numeric",2),rep("NULL",3),"factor",rep("NULL",2)))
names(pophist05) <- c("iscoral","ID","suitmod")

#load empirical data ("reference)
popemp05 <- read.csv(file="E:/pop_hr/hist_gpw/histgpw_v3_2005_suit.csv",row.names=NULL,header=T,
                      colClasses=c(rep("NULL",3),"numeric",rep("NULL",3),"factor"))
names(popemp05) <- c("ID","suitemp")

#merge historical and empirical by year with pointid
pop05 <- join_all(list(pophist05,popemp05),by=c("ID"))
pop05 <- pop05[!(pop05$iscoral==0),]

#confusion matrix
cmatrix_pop <- confusionMatrix(pop05$suitmod,pop05$suitemp)







###
##confusion matrix - overall
###

#use suitability csv, merge dataframes by pointid, filter to only coral sites. turn suitability column into factor. make confusion matrix
#apply to both empirical and historical model median for overlapping dates (2005)

#load historical models ("data")
suithist05 <- read.csv(file="F:/suit/hist/hist_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c("NULL",rep("numeric",2),rep("NULL",6),"factor","NULL"))
names(suithist05) <- c("iscoral","ID","suitmod")

#load empirical data ("reference")
suitemp05 <- read.csv(file="F:/suit/emp/emp_2005_suit.csv",row.names=NULL,header=F,skip=1,colClasses=c("NULL",rep("numeric",2),rep("NULL",6),"factor","NULL"))
names(suitemp05) <- c("iscoral","ID","suitemp")

#merge historical and empirical by year with pointid. remove non-coral rows with filter. 
suit05 <- join_all(list(suithist05,suitemp05),by=c("ID"))
suit05 <- suit05[!(suit05$iscoral==0),]

#confusion matrix
cmatrix_suit <- confusionMatrix(suit05$suitmod,suit05$suitemp)

suit05$suitmodn <- as.numeric(suit05$suitmod)
rocoverall <- roc(response=suit05$suitemp,predictor=suit05$suitmodn)
plot(rocoverall)




#plot confusion matrices
#inspo from https://www.reddit.com/r/rstats/comments/c6lvg0/confusion_matrix_caret_plotting_superior_to_base/

percent <- function(x, digits = 2, format = "f", ...) {      # display as percentage
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

draw_confusion_matrix <- function(cmtrx,title) {
  total <- sum(cmtrx$table)
  res <- as.numeric(cmtrx$table)
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title, cex.main=2)
  # create the matrix
  classes = colnames(cmtrx$table)
  #rect(150, 430, 240, 370, col=getColor("green", res[4]))
  rect(150, 430, 240, 370, col="#74C476")
  text(195, 435, classes[1], cex=1.2)
  #rect(250, 430, 340, 370, col=getColor("red", res[4]))
  rect(250, 430, 340, 370, col="#FC9272")
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  #rect(150, 305, 240, 365, col=getColor("red", res[4]))
  rect(150, 305, 240, 365, col="#FC9272")
  #rect(250, 305, 340, 365, col=getColor("green", res[4]))
  rect(250, 305, 340, 365, col="#74C476")
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)

  # add in the cmtrx results
  text(195, 400, percent(res[1]/total), cex=1.6, font=2, col='white')
  text(195, 335, percent(res[2]/total), cex=1.6, font=2, col='white')
  text(295, 400, percent(res[3]/total), cex=1.6, font=2, col='white')
  text(295, 335, percent(res[4]/total), cex=1.6, font=2, col='white')
  text(195, 390, res[1], cex=1, font=1, col='white')
  text(195, 325, res[2], cex=1, font=1, col='white')
  text(295, 390, res[3], cex=1, font=1, col='white')
  text(295, 325, res[4], cex=1, font=1, col='white')
  
  # add in the specifics
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cmtrx$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cmtrx$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cmtrx$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cmtrx$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cmtrx$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cmtrx$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cmtrx$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cmtrx$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cmtrx$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cmtrx$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information
  text(30, 35, names(cmtrx$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cmtrx$overall[1]), 3), cex=1.4)
  text(70, 35, names(cmtrx$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cmtrx$overall[2]), 3), cex=1.4)
}



draw_confusion_matrix(cmatrix_sst,title = "Confusion Matrix - SST")

draw_confusion_matrix(cmatrix_sst,title = "Confusion Matrix - DHW")

draw_confusion_matrix(cmatrix_oa,title = expression(paste("Confusion Matrix -",Omega,"arag")))

draw_confusion_matrix(cmatrix_storm,title = "Confusion Matrix - Storms")

draw_confusion_matrix(cmatrix_land,title = expression(paste("Confusion Matrix - Land Use")))

draw_confusion_matrix(cmatrix_pop,title = expression(paste("Confusion Matrix - Population Density")))

draw_confusion_matrix(cmatrix_suit,title= "Confusion Matrix - Overall Suitability")
