
library(plotrix)
library(maps)
library(lattice)
library(PBSmapping)
library(mapdata)    
library(maptools)   
library(mapproj)



#============================================
# BSFRF DATA
#============================================
#==species 1 = opilio; 2 = bairdi; 3 = hybrid
#==sex 1 = male; 2 = female
dat09<-data.frame(read.csv("data/2009.csv"))
dat10<-read.csv("data/2010.csv")
dat16<-read.csv("data/2016.csv")
dat17<-read.csv("data/2017.csv")
dat18<-read.csv("data/2018.csv")
AllStnLoc<-read.csv("data/AllStnLoc.csv")
AllStation<-read.csv("data/AllStation.csv")

#==compare length compositions and densities from each of the areas from 2009 to the length comps from the surveye
par(mfrow=c(3,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
hist(as.numeric(as.character(dat09[dat09$Species==1 & dat09$Sex==1 & dat09$Area=='SE1',]$CarapW)),las=1,xlim=c(0,150),xaxt='n',main='',ylim=c(0,800))
legend("topleft",bty='n',legend="All males; area SE1")
hist(as.numeric(as.character(dat09[dat09$Species==1 & dat09$Sex==1 & dat09$Area=='CN2',]$CarapW)),las=1,xlim=c(0,150),xaxt='n',main='',ylim=c(0,800))
legend("topleft",bty='n',legend="All males; area CN2")
hist(as.numeric(as.character(dat09[dat09$Species==1 & dat09$Sex==1 & dat09$Area=='NW3',]$CarapW)),las=1,xlim=c(0,150),main='',ylim=c(0,800))
legend("topleft",bty='n',legend="All males; area NW3")

#==compare length compositions and densities from each of the areas from 2009 to the length comps from the surveye
par(mfrow=c(3,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
hist(as.numeric(as.character(dat10[dat10$Species==1 & dat10$Sex==1 ,]$CarapW)),las=1,xlim=c(0,150),xaxt='n',main='',ylim=c(0,2500))
legend("topleft",bty='n',legend="All males")
hist(as.numeric(as.character(dat10[dat10$Species==1 & dat10$Sex==2 ,]$CarapW)),las=1,xlim=c(0,150),xaxt='n',main='',ylim=c(0,2500))
legend("topleft",bty='n',legend="All females")


#==plot map of all the side by side with year as color
#==make figure showing the relative size of area swept
#==make figure showing the length compositions by gear by 'experiment'
#==run regression with proportion caught by size against temperature and bottom type 


png("BSFRF_map.png",height=8,width=8,res=1200,units='in')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85')
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE)

axis(side=1)
axis(side=2,las=1)
mtext(side=1,"Longitude",line=2)
mtext(side=2,"Latitude",line=2)


use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)

st_17<-unique(dat17$Station[dat17$Species=='OPI'])
for(x in 1:length(st_17))
{
  temp<-dat17[dat17$Station==st_17[x],]
  points(x=temp$Longitude[1]+rnorm(1,jit_u,jit_sd),y=temp$Latitude[1]+rnorm(1,jit_u,jit_sd),pch=use_pch,col=4)
}

st_16<-unique(dat16$NMFS_Stn[dat16$Species=='OPI'])
for(x in 1:length(st_16))
{
  temp<-dat16[dat16$NMFS_Stn==st_16[x],]
  points(x=as.numeric(as.character(temp$MidTowLongitude[1]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp$MidTowLatitude[1]))+rnorm(1,jit_u,jit_sd),col=3,pch=use_pch)
}

st_18<-unique(dat18$Station[dat18$Species=='OPI'])
for(x in 1:length(st_18))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_18[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=5,pch=use_pch)
}

st_10<-as.character(unique(dat10$NMFSstn[dat10$Species==1]))
for(x in 1:length(st_10))
{
  if(nchar(st_10[x])==3)
  {
    st_10[x]<-paste(substr(st_10[x],1,1),substr(st_10[x],2,3),sep="-")
  }
}


for(x in 1:length(st_10))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_10[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=2,pch=use_pch)
}

st_09<-unique(dat09$NMFS.station[dat09$Species==1])
for(x in 1:length(st_09))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_09[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=1,pch=use_pch)
}

legend('bottomleft',bty='n',col=c(1,2,3,4,5),legend=c(2009,2010,2016,2017,2018),pch=16)
dev.off()


#====================================================================
#==COMPARE NMFS AND BSFRF NUMBERS AT LENGTH WITHIN THE SAME AREAS==
#==plot 5 x 2 of numbers at length in each area in each year
#==this pulls in NMFS survey data
#==to do yearly comparisons and make maps
#==ask me if you want this file, Jim
use_survey<-0
if(use_survey==1)
{
  survDAT<-read.csv("C:/SnowCrab2019/new data 2019/EBSCrab_Haul.csv",header=T,skip=5)
  
  drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
  SurvYR<-unique(drvYear)
  AllStation<-unique(survDAT$GIS_STATION)
  
  #==plot GIS stations
  AllStnLoc<-matrix(ncol=2,nrow=length(AllStation))
  for(w in 1:length(AllStation))
  {
    temp<-survDAT[survDAT$GIS_STATION==AllStation[w],]
    AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
    AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
  }
  
  #write.csv(AllStnLoc,"AllStnLoc.csv",row.names=FALSE)
  #write.csv(AllStation,"AllStation.csv",row.names=FALSE)
  #========================
  #==find survey densities (total)
  #========================
  nmiSurv<-140350
  DensityM<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
  DensityM78<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
  DensityM101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
  DensityF<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
  DensityFmat<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
  StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
  for(y in 1:length(SurvYR))
  {
    yrDAT<-survDAT[drvYear==SurvYR[y],]
    fileyr<-SurvYR[y]
    stationsUNQ<-(unique(yrDAT$GIS_STATION))
    #==density at station
    for(j in 1:length(stationsUNQ))
    {
      stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
      StationYr[y,j]<-as.character(stationsUNQ[j])
      Hauls<-(unique(stationALL$HAUL))
      Males<-stationALL[stationALL$SEX==1,]
      Males78<-stationALL[stationALL$SEX==1 & stationALL$WIDTH>77,]
      Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
      Females<-stationALL[stationALL$SEX==2,]
      FemalesMat<-stationALL[stationALL$SEX==2 & stationALL$EGG_CONDITION>0,]
      #==densities across hauls in crabs per km^2
      tempDensM<-NULL
      tempDensM78<-NULL
      tempDensM101<-NULL
      tempDensF<-NULL
      for(k in 1:length(Hauls))
      {
        SampFactM<-Males$SAMPLING_FACTOR[which(Males$HAUL==Hauls[k])[1]] 
        AreaSweptM<-Males$AREA_SWEPT_VARIABLE[which(Males$HAUL==Hauls[k])[1]]
        tempDensM<-length(Males$HAUL==Hauls[k])*SampFactM/AreaSweptM
        
        SampFactM<-Males78$SAMPLING_FACTOR[which(Males78$HAUL==Hauls[k])[1]] 
        AreaSweptM<-Males78$AREA_SWEPT_VARIABLE[which(Males78$HAUL==Hauls[k])[1]]
        tempDensM78<-length(Males78$HAUL==Hauls[k])*SampFactM/AreaSweptM
        
        SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
        AreaSweptM<-Males101$AREA_SWEPT_VARIABLE[which(Males101$HAUL==Hauls[k])[1]]
        tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
        
        SampFactF<-Females$SAMPLING_FACTOR[which(Females$HAUL==Hauls[k])[1]] 
        AreaSweptF<-Females$AREA_SWEPT_VARIABLE[which(Females$HAUL==Hauls[k])[1]]
        tempDensF<-length(Females$HAUL==Hauls[k])*SampFactF/AreaSweptF
        
        SampFactFm<-FemalesMat$SAMPLING_FACTOR[which(FemalesMat$HAUL==Hauls[k])[1]] 
        AreaSweptF<-FemalesMat$AREA_SWEPT_VARIABLE[which(FemalesMat$HAUL==Hauls[k])[1]]
        tempDensFm<-length(FemalesMat$HAUL==Hauls[k])*SampFactF/AreaSweptF
      }
      DensityM[y,j]<-mean(tempDensM)
      DensityM78[y,j]<-mean(tempDensM78)
      DensityM101[y,j]<-mean(tempDensM101)
      DensityF[y,j]<-mean(tempDensF)
      DensityFmat[y,j]<-mean(tempDensFm)
    }
  }
}

len_bin<-seq(27.5,132.5,5)

NMFS_plot<-function(in_yr,BSFRF_stat)
{
  yrDAT<-survDAT[survDAT$AKFIN_SURVEY_YEAR==in_yr & !is.na(match(survDAT$GIS_STATION,BSFRF_stat)) & survDAT$SEX==1,]
  n_at_len<-rep(0,length(len_bin))
  stationsUNQ<-(unique(yrDAT$GIS_STATION)) 
  dens<-rep(0,length(stationsUNQ))
  
  #==add up N at length
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j] ,]
    Hauls<-(unique(stationALL$HAUL))
    for(k in 1:nrow(stationALL))
    {
      if(!is.na(stationALL[k,]$SAMPLING_FACTOR))
      {
        n_at_len[which(abs(len_bin-stationALL[k,]$WIDTH)==min(abs(len_bin-stationALL[k,]$WIDTH)))] <-
          n_at_len[which(abs(len_bin-stationALL[k,]$WIDTH)==min(abs(len_bin-stationALL[k,]$WIDTH)))] + stationALL[k,]$SAMPLING_FACTOR
      }
    } 
  #==calculate density by station
      tempDensM<-NULL
      for(k in 1:length(Hauls))
      {
        SampFactM<-yrDAT$SAMPLING_FACTOR[which(yrDAT$HAUL==Hauls[k])[1]] 
        AreaSweptM<-yrDAT$AREA_SWEPT_VARIABLE[which(yrDAT$HAUL==Hauls[k])[1]]
        tempDensM<-append(tempDensM,sum(yrDAT$HAUL==Hauls[k],na.rm=T)*SampFactM/AreaSweptM)
      }   
      #if(!is.na(tempDensM))
       dens[j]<-mean(tempDensM,na.rm=T)
  }
  dumb_hist<-NULL
  for(x in 1:length(n_at_len))
    dumb_hist<-append(dumb_hist,rep(len_bin[x],round(n_at_len[x])))
  list(n_at_len,dumb_hist,dens)
}


NMFS_2009<-NMFS_plot(2009,st_09)
NMFS_2010<-NMFS_plot(2010,st_10)
NMFS_2016<-NMFS_plot(2016,st_16)
NMFS_2017<-NMFS_plot(2017,st_17)
NMFS_2018<-NMFS_plot(in_yr=2018,BSFRF_stat=st_18)

#==plot captured numbers at length
png("BSFRF_NMFS_comp.png",height=4,width=8,res=1200,units='in')
input_y<-4000
input_x<-160
par(mfcol=c(2,5),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
hist(as.numeric(as.character(dat09[dat09$Species==1 & dat09$Sex==1,]$CarapW)),las=1,
     xlim=c(0,input_x),xaxt='n',main='',ylim=c(0,input_y),breaks=seq(0,150,5))
hist(NMFS_2009[[2]],las=1,xlim=c(0,input_x),main='',ylim=c(0,input_y),breaks=seq(0,150,5))

hist(as.numeric(as.character(dat10[dat10$Species==1 & dat10$Sex==1,]$CarapW)),las=1,
     xlim=c(0,input_x),xaxt='n',main='',ylim=c(0,input_y),breaks=seq(0,150,5),yaxt='n')
hist(NMFS_2010[[2]],las=1,xlim=c(0,input_x),main='',ylim=c(0,input_y),breaks=seq(0,150,5),yaxt='n')

hist(as.numeric(as.character(dat16[dat16$Sex==1 & dat16$Species=="OPI",]$Carapace)),las=1,
     xlim=c(0,input_x),xaxt='n',main='',ylim=c(0,input_y),breaks=seq(0,150,5),yaxt='n')
hist(NMFS_2016[[2]],las=1,xlim=c(0,input_x),main='',ylim=c(0,input_y),breaks=seq(0,150,5),yaxt='n')

hist(as.numeric(as.character(dat17[dat17$Sex==1 & dat17$Species=="OPI",]$Carapace)),las=1,
     xlim=c(0,input_x),xaxt='n',main='',ylim=c(0,input_y),breaks=seq(0,160,5),yaxt='n')
hist(NMFS_2017[[2]],las=1,xlim=c(0,input_x),main='',ylim=c(0,input_y),breaks=seq(0,150,5),yaxt='n')

hist(as.numeric(as.character(dat18[dat18$Sex==1 & dat18$Species=="OPI",]$Carapace)),las=1,
     xlim=c(0,input_x),xaxt='n',main='',ylim=c(0,input_y),breaks=seq(0,160,5),yaxt='n')
hist(NMFS_2018[[2]],las=1,xlim=c(0,input_x),main='',ylim=c(0,input_y),breaks=seq(0,150,5),yaxt='n')

dev.off()

#=========================================================
# calculate densities to match the length composition data
# for both NMFS and BSFRF areas surveyed
#=========================================================

#+++++++++++++++++++++++++++++++++++++++++
# ADD CALCULATION OF NUMBERS AT LENGTH
# NEED TO HAVE THE SAMPLING FACTOR INFLUENCING THE NUMBERS
# THAT CREATE THE LENGTH COMPSS
# IT IS WRONG UNTIL THAT IS DONE
# (I don't remember if I fixed this, so check it)
#+=========++++++++++++++++++++++++++++++++++++++++

#==density 2017
st_17<-unique(dat17$Station[dat17$Species=='OPI'])
dens_17<-rep(0,length(st_17))  
n_at_len_17<-rep(0,length(len_bin))

for(x in 1:length(st_17))
{
  temp<-dat17[dat17$Station==st_17[x] & dat17$Species=="OPI" & dat17$Sex==1,]
  Hauls<-(unique(temp$Haul))
  tempDensM<-NULL
  for(y in 1:length(Hauls))
  {
  SampFactM<-temp$SampFactor[which(temp$Haul==Hauls[y])][1] 
  AreaSweptM<-temp$Aswept_nm2[which(temp$Haul==Hauls[y])[1]]
  tempDensM<-append(tempDensM,sum(temp$Haul==Hauls[y],na.rm=T)*SampFactM/AreaSweptM)
  }
  if(!is.na(tempDensM))
   dens_17[x]<-mean(tempDensM,na.rm=T)
 
  #==add up N at length
  if(nrow(temp)>0)
  {
  for(k in 1:nrow(temp))
    {
      if(!is.na(temp[k,]$SampFactor))
      {
        n_at_len_17[which(abs(len_bin-temp[k,]$Carapace)==min(abs(len_bin-temp[k,]$Carapace)))] <-
          n_at_len_17[which(abs(len_bin-temp[k,]$Carapace)==min(abs(len_bin-temp[k,]$Carapace)))] + temp[k,]$SampFactor
      }
    } 
  }
}

#==density 2016
st_16<-unique(dat16$NMFS_Stn[dat16$Species=='OPI'])
dens_16<-rep(0,length(st_16))
n_at_len_16<-rep(0,length(len_bin))

for(x in 1:length(st_16))
{
  temp<-dat16[dat16$NMFS_Stn==st_16[x] & dat16$Species=="OPI" & dat16$Sex==1,]
  Hauls<-(unique(temp$Tow))
  tempDensM<-NULL
  for(y in 1:length(Hauls))
  {
    SampFactM<-temp$SampFactor[which(temp$Tow==Hauls[y])][1] 
    AreaSweptM<-temp$Aswept_nm2[which(temp$Tow==Hauls[y])[1]]
    tempDensM<-append(tempDensM,sum(temp$Tow==Hauls[y],na.rm=T)*SampFactM/AreaSweptM)
  }
  if(!is.na(tempDensM))
    dens_16[x]<-mean(tempDensM,na.rm=T)
  
  #==add up N at length
  if(nrow(temp)>0)
  {
  for(k in 1:nrow(temp))
  {
    if(!is.na(temp[k,]$SampFactor))
    {
      n_at_len_16[which(abs(len_bin-temp[k,]$Carapace)==min(abs(len_bin-temp[k,]$Carapace)))] <-
        n_at_len_16[which(abs(len_bin-temp[k,]$Carapace)==min(abs(len_bin-temp[k,]$Carapace)))] + temp[k,]$SampFactor
    }
  } 
 }
}

#==density 2018
st_18<-unique(dat18$Station[dat18$Species=='OPI'])
dens_18<-rep(0,length(st_18))
n_at_len_18<-rep(0,length(len_bin))

for(x in 1:length(st_18))
{
  temp<-dat18[dat18$Station==st_18[x] & dat18$Species=="OPI" & dat18$Sex==1,]
  Hauls<-(unique(temp$Haul))
  tempDensM<-NULL
  for(y in 1:length(Hauls))
  {
    SampFactM<-temp$SUBSAMPLING.FACTOR[which(temp$Haul==Hauls[y])] [1] 
    AreaSweptM<-temp$Aswept_nm2[which(temp$Haul==Hauls[y])[1]]
    tempDensM<-append(tempDensM,sum(temp$Haul==Hauls[y],na.rm=T)*SampFactM/AreaSweptM)
  }
  if(!is.na(tempDensM))
    dens_18[x]<-mean(tempDensM,na.rm=T)
  
  #==add up N at length
  if(nrow(temp)>0)
  {
  for(k in 1:nrow(temp))
  {
    if(!is.na(temp[k,]$SUBSAMPLING.FACTOR))
    {
      n_at_len_18[which(abs(len_bin-temp[k,]$Carapace)==min(abs(len_bin-temp[k,]$Carapace)))] <-
        n_at_len_18[which(abs(len_bin-temp[k,]$Carapace)==min(abs(len_bin-temp[k,]$Carapace)))] + temp[k,]$SUBSAMPLING.FACTOR
    }
  } 
  }
}


#==density 2010
st_10<-as.character(unique(dat10$NMFSstn[dat10$Species==1]))
dens_10<-rep(0,length(st_10))
n_at_len_10<-rep(0,length(len_bin))

for(x in 1:length(st_10))
{
  temp<-dat10[dat10$NMFSstn==st_10[x] & dat10$Species==1 & dat10$Sex==1,]
  Hauls<-(unique(temp$Haul.))
  tempDensM<-NULL
  for(y in 1:length(Hauls))
  {
    SampFactM<-temp$NumCrab[1] 
    AreaSweptM<-temp$AreaSweptnm2[which(temp$Haul.==Hauls[y])[1]]
    tempDensM<-append(tempDensM,sum(temp$Haul.==Hauls[y],na.rm=T)*SampFactM/AreaSweptM)
  }
  if(!is.na(tempDensM))
    dens_10[x]<-mean(tempDensM,na.rm=T)
  #==add up N at length
  if(nrow(temp)>0)
  {
  for(k in 1:nrow(temp))
  {
    if(!is.na(temp[k,]$NumCrab))
    {
      n_at_len_10[which(abs(len_bin-temp[k,]$CarapW)==min(abs(len_bin-temp[k,]$CarapW)))] <-
        n_at_len_10[which(abs(len_bin-temp[k,]$CarapW)==min(abs(len_bin-temp[k,]$CarapW)))] + temp[k,]$NumCrab
    }
  } 
 }  
}

#==density 2009
st_09<-unique(dat09$NMFS.station[dat09$Species==1])
dens_09<-rep(0,length(st_09))
n_at_len_09<-rep(0,length(len_bin))
dat09$CarapW<-as.numeric(as.character(dat09$CarapW))

for(x in 1:length(st_09))
{
  temp<-dat09[dat09$NMFS.station==st_09[x] & dat09$Species==1 & dat09$Sex==1,]
  Hauls<-(unique(temp$Haul))
  tempDensM<-NULL
  for(y in 1:length(Hauls))
  {
    SampFactM<-temp$Num.Crab[1] 
    AreaSweptM<-temp$Aswept[which(temp$Haul==Hauls[y])[1]]
    tempDensM<-append(tempDensM,sum(temp$Haul==Hauls[y],na.rm=T)*SampFactM/AreaSweptM)
  }
  if(!is.na(tempDensM))
    dens_09[x]<-mean(tempDensM,na.rm=T)
  
  #==add up N at length
  if(nrow(temp)>0)
  {
  for(k in 1:nrow(temp))
  {
    if(!is.na(temp[k,]$Num.Crab))
    {
      n_at_len_09[which(abs(len_bin-temp[k,]$CarapW)==min(abs(len_bin-temp[k,]$CarapW)))] <-
        n_at_len_09[which(abs(len_bin-temp[k,]$CarapW)==min(abs(len_bin-temp[k,]$CarapW)))] + temp[k,]$Num.Crab
    }
  } 
  }  
}


avg_den_09_BSFRF<-mean(dens_09,na.rm=T)
avg_den_09_NMFS<-mean(NMFS_2009[[3]][NMFS_2009[[3]]!="NaN"],na.rm=T)

avg_den_10_BSFRF<-mean(dens_10,na.rm=T)
avg_den_10_NMFS<-mean(NMFS_2010[[3]][NMFS_2010[[3]]!="NaN"],na.rm=T)

avg_den_16_BSFRF<-mean(dens_16,na.rm=T)
avg_den_16_NMFS<-mean(NMFS_2016[[3]][NMFS_2016[[3]]!="NaN"],na.rm=T)

avg_den_17_BSFRF<-mean(dens_17,na.rm=T)
avg_den_17_NMFS<-mean(NMFS_2017[[3]][NMFS_2017[[3]]!="NaN"],na.rm=T)

avg_den_18_BSFRF<-mean(dens_18,na.rm=T)
avg_den_18_NMFS<-mean(NMFS_2018[[3]][NMFS_2018[[3]]!="NaN"],na.rm=T)


#==plot proportion at length
#==another place whewre I should have used a function, but copied codee!!!
inylim<-1

make_dens<-function(input)
{
  min_div<-1
  max_div<-0
  div<-(min_div+max_div)/2
 for(x in 1:25)
 {
  temp<-input*((max_div+min_div)/2)
  if(sum(temp)>1)
    min_div<-(min_div+max_div)/2
  if(sum(temp)<1)
    max_div<-(max_div+min_div)/2
 }
  return(min_div)
}


tot_num_BSFRF<-c(sum(n_at_len_09),
                 sum(n_at_len_10),
                 sum(n_at_len_16),
                 sum(n_at_len_17),
                 sum(n_at_len_18))

tot_num_NMFS<-c(sum(NMFS_2009[[1]]),
                sum(NMFS_2010[[1]]),
                sum(NMFS_2016[[1]]),
                sum(NMFS_2017[[1]]),
                sum(NMFS_2018[[1]]))
input<-matrix(c(tot_num_BSFRF,tot_num_NMFS),nrow=2,byrow=T)

png("Sample_sizes.png",height=4,width=5,res=500,units='in')
par(mar=c(4,5,1,1))
col1<-'light grey'
col2<-'dark grey'
barplot(input,beside=T,names.arg = c(2009,2010,2016,2017,2018),las=1,col=c(col1,col2))
mtext(side=1,"Year",line=2.5)
mtext(side=2,"Total crab captured",line=3.5)
legend("topright",bty='n',col=c(col1,col2),legend=c("NMFS","BSFRF"),pch=15)
dev.off()
#================================
# Plotting inferred selectivities
#================================
png("BSFRF_NMFS_select.png",height=8,width=5,res=500,units='in')
xlim_n<-27.5
xlim_up<-132.5


#==2009
par(mfrow=c(5,2),mar=c(.1,.1,.1,.1),oma=c(4,4,1,4))
ab<-make_dens(n_at_len_09)*n_at_len_09
dc1<-hist(NMFS_2009[[2]],breaks=seq(0,160,5),plot=F)
dc<-make_dens(dc1$counts)*dc1$counts

BS_plot<-ab*avg_den_09_BSFRF
NM_plot<-dc*avg_den_09_NMFS

plot(BS_plot~len_bin,type='l',las=1,xaxt='n',xlim=c(xlim_n,xlim_up))
lines(NM_plot~dc1$mids,lty=2)
legend("topright",legend=c(2009,"BSFRF","NMFS"),lty=c(NA,1,2),col=c(NA,1,1),bty='n')

sel_09<-NM_plot[which(!is.na(match(dc1$mids,len_bin)))]/BS_plot
plot(sel_09~len_bin,type='l',ylim=c(0,inylim),las=1,yaxt='n',xaxt='n')
axis(side=4,las=1)


#==2010
ab<-make_dens(n_at_len_10)*n_at_len_10
dc1<-hist(NMFS_2010[[2]],breaks=seq(0,160,5),plot=F)
dc<-make_dens(dc1$counts)*dc1$counts

BS_plot<-ab*avg_den_10_BSFRF
NM_plot<-dc*avg_den_10_NMFS

plot(BS_plot~len_bin,type='l',las=1,xaxt='n',xlim=c(xlim_n,xlim_up),col=2)
lines(NM_plot~dc1$mids,lty=2,col=2)
legend("topright",legend=2010,bty='n')

sel_10<-NM_plot[which(!is.na(match(dc1$mids,len_bin)))]/BS_plot
plot(sel_10~len_bin,type='l',ylim=c(0,inylim),las=1,yaxt='n',xaxt='n',col=2)
axis(side=4,las=1)

#==2016
ab<-make_dens(n_at_len_16)*n_at_len_16
dc1<-hist(NMFS_2016[[2]],breaks=seq(0,160,5),plot=F)
dc<-make_dens(dc1$counts)*dc1$counts

BS_plot<-ab*avg_den_16_BSFRF
NM_plot<-dc*avg_den_16_NMFS

plot(BS_plot~len_bin,type='l',las=1,xaxt='n',xlim=c(xlim_n,xlim_up),col=3)
lines(NM_plot~dc1$mids,lty=2,col=3)
legend("topright",legend=2016,bty='n')

sel_16<-NM_plot[which(!is.na(match(dc1$mids,len_bin)))]/BS_plot
plot(sel_16~len_bin,type='l',ylim=c(0,inylim),las=1,yaxt='n',xaxt='n',col=3)
axis(side=4,las=1)

#==2017
ab<-make_dens(n_at_len_17)*n_at_len_17
dc1<-hist(NMFS_2017[[2]],breaks=seq(0,160,5),plot=F)
dc<-make_dens(dc1$counts)*dc1$counts

BS_plot<-ab*avg_den_17_BSFRF
NM_plot<-dc*avg_den_17_NMFS

plot(BS_plot~len_bin,type='l',las=1,xaxt='n',xlim=c(xlim_n,xlim_up),col=4)
lines(NM_plot~dc1$mids,lty=2,col=4)
legend("topright",legend=2017,bty='n')

sel_17<-NM_plot[which(!is.na(match(dc1$mids,len_bin)))]/BS_plot
plot(sel_17~len_bin,type='l',ylim=c(0,inylim),las=1,yaxt='n',xaxt='n',col=4)
axis(side=4,las=1)

#==2018
ab<-make_dens(n_at_len_18)*n_at_len_18
dc1<-hist(NMFS_2018[[2]],breaks=seq(0,160,5),plot=F)
dc<-make_dens(dc1$counts)*dc1$counts

BS_plot<-ab*avg_den_18_BSFRF
NM_plot<-dc*avg_den_18_NMFS

plot(BS_plot~len_bin,type='l',las=1,xaxt='n',xlim=c(xlim_n,xlim_up),col=5)
axis(side=1)
lines(NM_plot~dc1$mids,lty=2,col=5)
legend("topright",legend=2018,bty='n')

sel_18<-NM_plot[which(!is.na(match(dc1$mids,len_bin)))]/BS_plot
plot(sel_18~len_bin,type='l',ylim=c(0,inylim),las=1,yaxt='n',xaxt='n',col=5)
axis(side=4,las=1)
axis(side=1)
mtext(side=4,outer=T,"Selectivity",line=2.5)
mtext(side=2,outer=T,"Density at length (crab per nm^2)",line=3)
mtext(side=1,outer=T,"Carapace width (mm)",line=2.5)
dev.off()

png("BSFRF_NMFS_select_allyrs.png",height=4,width=5,res=500,units='in')

plot(sel_09~len_bin,type='l',ylim=c(0,1.50),las=1,ylab="Selectivity",xlab="Carapace width (mm)")
lines(sel_10~len_bin,col=2)
lines(sel_16~len_bin,col=3)
lines(sel_17~len_bin,col=4)
lines(sel_18~len_bin,col=5)
legend("top",col=seq(1,5),bty='n',legend=c(2009,2010,2016,2017,2018),pch=15,cex=.8)
dev.off()

#=========================================================
# calculate the empirical availabilities for each year
#=========================================================


#=========================
# Zeroes experiment
#========================
nobs<-400
obs_dens<-exp(rnorm(nobs,3,1))
hist(obs_dens,breaks=seq(0,max(obs_dens),length.out=100))
area_obs<-500

pop_est<-mean(obs_dens)*area_obs*nobs

mod_obs<-obs_dens[-obs_dens<1]
mod_pop_est<-mean(mod_obs)*length(mod_obs)*area_obs

#====================
#==individual maps
#====================
png("BSFRF_map_09.png",height=8,width=8,res=1200,units='in')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85')
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE)

axis(side=1)
axis(side=2,las=1)
mtext(side=1,"Longitude",line=2)
mtext(side=2,"Latitude",line=2)

use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)

st_09<-unique(dat09$NMFS.station[dat09$Species==1])
for(x in 1:length(st_09))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_09[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=5,pch=use_pch)
}

legend('bottomleft',bty='n',col=c(5,2,1,3,4),legend=c(2009,2010,2016,2017,2018),pch=16)
dev.off()


#==2010

png("BSFRF_map_10.png",height=8,width=8,res=1200,units='in')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85')
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE)

axis(side=1)
axis(side=2,las=1)
mtext(side=1,"Longitude",line=2)
mtext(side=2,"Latitude",line=2)

use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)

st_10<-as.character(unique(dat10$NMFSstn[dat10$Species==1]))
for(x in 1:length(st_10))
{
  if(nchar(st_10[x])==3)
  {
    st_10[x]<-paste(substr(st_10[x],1,1),substr(st_10[x],2,3),sep="-")
  }
}


for(x in 1:length(st_10))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_10[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=2,pch=use_pch)
}

legend('bottomleft',bty='n',col=c(5,2,1,3,4),legend=c(2009,2010,2016,2017,2018),pch=16)
dev.off()


#==2016
png("BSFRF_map_16.png",height=8,width=8,res=1200,units='in')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85')
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE)

axis(side=1)
axis(side=2,las=1)
mtext(side=1,"Longitude",line=2)
mtext(side=2,"Latitude",line=2)

use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)

st_16<-unique(dat16$NMFS_Stn[dat16$Species=='OPI'])
for(x in 1:length(st_16))
{
  temp<-dat16[dat16$NMFS_Stn==st_16[x],]
  points(x=as.numeric(as.character(temp$MidTowLongitude[1]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp$MidTowLatitude[1]))+rnorm(1,jit_u,jit_sd),col=1,pch=use_pch)
}


legend('bottomleft',bty='n',col=c(5,2,1,3,4),legend=c(2009,2010,2016,2017,2018),pch=16)
dev.off()

#==2017
png("BSFRF_map_17.png",height=8,width=8,res=1200,units='in')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85')
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE)

axis(side=1)
axis(side=2,las=1)
mtext(side=1,"Longitude",line=2)
mtext(side=2,"Latitude",line=2)

use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)


st_17<-unique(dat17$Station[dat17$Species=='OPI'])
for(x in 1:length(st_17))
{
  temp<-dat17[dat17$Station==st_17[x],]
  points(x=temp$Longitude[1]+rnorm(1,jit_u,jit_sd),y=temp$Latitude[1]+rnorm(1,jit_u,jit_sd),pch=use_pch,col=3)
}


legend('bottomleft',bty='n',col=c(5,2,1,3,4),legend=c(2009,2010,2016,2017,2018),pch=16)
dev.off()


#==2018
png("BSFRF_map_18.png",height=8,width=8,res=1200,units='in')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85')
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE)

axis(side=1)
axis(side=2,las=1)
mtext(side=1,"Longitude",line=2)
mtext(side=2,"Latitude",line=2)

use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)


st_18<-unique(dat18$Station[dat18$Species=='OPI'])
for(x in 1:length(st_18))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_18[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=4,pch=use_pch)
}


legend('bottomleft',bty='n',col=c(5,2,1,3,4),legend=c(2009,2010,2016,2017,2018),pch=16)
dev.off()


#########################################
# single figure with maps
##########################################
inmar<-c(rep(0,4))
png("BSFRF_map_plots.png",height=8,width=8,res=500,units='in')
par(mfrow=c(2,3),xpd=FALSE,oma=c(3,3,1,1))
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',mar=inmar)
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE,mar=inmar)


axis(side=2,las=1)

use_pch<-16
jit_u<-0
jit_sd<-.01
set.seed(1)

st_09<-unique(dat09$NMFS.station[dat09$Species==1])
for(x in 1:length(st_09))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_09[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=1,pch=use_pch)
}


#==2010

state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',mar=inmar)
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE,mar=inmar)

st_10<-as.character(unique(dat10$NMFSstn[dat10$Species==1]))
for(x in 1:length(st_10))
{
  if(nchar(st_10[x])==3)
  {
    st_10[x]<-paste(substr(st_10[x],1,1),substr(st_10[x],2,3),sep="-")
  }
}


for(x in 1:length(st_10))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_10[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=2,pch=use_pch)
}

#==2016
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',mar=inmar)
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE,mar=inmar)


st_16<-unique(dat16$NMFS_Stn[dat16$Species=='OPI'])
for(x in 1:length(st_16))
{
  temp<-dat16[dat16$NMFS_Stn==st_16[x],]
  points(x=as.numeric(as.character(temp$MidTowLongitude[1]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp$MidTowLatitude[1]))+rnorm(1,jit_u,jit_sd),col=3,pch=use_pch)
}


#==2017
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',mar=inmar)
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE,mar=inmar)

axis(side=1)
axis(side=2,las=1)
mtext(side=2,"Latitude",line=0,outer=T)

st_17<-unique(dat17$Station[dat17$Species=='OPI'])
for(x in 1:length(st_17))
{
  temp<-dat17[dat17$Station==st_17[x],]
  points(x=temp$Longitude[1]+rnorm(1,jit_u,jit_sd),y=temp$Latitude[1]+rnorm(1,jit_u,jit_sd),pch=use_pch,col=4)
}


#==2018
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',mar=inmar)
abline(h=seq(52,65,1/3),v=seq(-177,-155,.6),col='grey')
state.map <- map('worldHires', xlim=c(-175, -160), ylim=c(53, 62),
                 plot = TRUE, fill = TRUE,col='grey85',add=TRUE,mar=inmar)

axis(side=1)
mtext(side=1,"Longitude",line=0,outer=T)


st_18<-unique(dat18$Station[dat18$Species=='OPI'])
for(x in 1:length(st_18))
{
  temp<-AllStnLoc[which(!is.na(match(AllStation,st_18[x]))),]
  points(x=as.numeric(as.character(temp[2]))+rnorm(1,jit_u,jit_sd),y=as.numeric(as.character(temp[1]))+rnorm(1,jit_u,jit_sd),col=5,pch=use_pch)
}

plot.new()
legend('center',bty='n',col=c(1,2,3,4,5),legend=c(2009,2010,2016,2017,2018),pch=16,cex=1.5)


dev.off()



