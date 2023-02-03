#FINDING low rank and sparse features in hourly, third-octave band spectra in SanctSound project

# purpose: separate transients and better estimates of ambient

#-----------------------------------------------------------------------------------------
#NOTES
#tested on SB01 and CI01 to run rrpca (robust PCA)
# low-rank spectra- how does this relate to ships vs wind speed
# sparse spectra-   how does this relate to biological events
# Then run SVD to find representative specta- how does this relate to percentile spectra

#UPDATES
# 2021.06.24-- updated to use with SB02 to quantify major shift in conditions, 
# can we determine is shift was in transient or chronic?
#NOTE: hourly SPL data are combined with weather, and only daily metrics of AIS ships and biological detentions

#-----------------------------------------------------------------------------------------
rm(list=ls())

library(data.table)
library(rpca)
library(rsvd)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(dplyr)
library(plyr)

#-----------------------------------------------------------------------------------------
#INPUT FILES
#-----------------------------------------------------------------------------------------
tDir  =  "E:\\RESEARCH\\SanctSound\\"
inDir = paste0(tDir,"data2\\combineFiles3_Detections")
cDir  = "E:\\CODE\\GitHub\\MM_SanctSound_SourceSeparation\\"

#eventually, create a loop through the sites...
sanct = "SB"   #"SB"
site  = "SB02" #"SB02"

#1/3 octave band data
#-----------------------------------------------------------------------------------------
o3b = fread(input=paste0(cDir, "O3OB_SanctSound.csv"))

#TOL files by site-deployment
#-----------------------------------------------------------------------------------------
sdir     = paste0(tDir, "data\\",sanct,"\\", site,sep="\\")
TOLFiles = list.files(sdir,pattern = "TOL_1h.csv",recursive = T)
PSDFiles = list.files(sdir,pattern = "PSD_1h.csv",recursive = T)
uSites   = unique( sapply( strsplit(basename(TOLFiles),"_") , `[`, 2) )


#summary of detection data to compare to results (output of 3_combineFiles_detections)
#-----------------------------------------------------------------------------------------
if (site == "SB02") {
  inFiles  = list.files(inDir,sanct,full.names = T)
  load(inFiles[5])#HOUR data
  load(inFiles[1])#DAY data
  
  # select time periods of interest
  timesInterst = read.csv( paste0(tDir,"data\\AnalysisPeriods.csv") )
  tmpAP   = timesInterst[timesInterst$Site == site,]
  dataAhrt  = dataAhr[ (as.Date(dataAhr$Day)  >= as.Date(tmpAP$DateStart,format = "%m/%d/%Y") 
                        & as.Date(dataAhr$Day) <= as.Date(tmpAP$DateEnd,format = "%m/%d/%Y") ),]
  dataAhrt  = dataAhrt[!is.na(dataAhrt$Deployment),] #had to remove NAs, effort control, but does not work in the analysis
  
  dataAdyt  = iData[ (as.Date(iData$Day)  >= as.Date(tmpAP$DateStart,format = "%m/%d/%Y") 
                      & as.Date(iData$Day) <= as.Date(tmpAP$DateEnd,format = "%m/%d/%Y") ),]
  dataAdyt  = dataAdyt[!is.na(dataAdyt$Deployment),] #had to remove NAs, effort control, but does not work in the analysis
  rm(tmpAP)
  
}
#narw
length( dataAdyt [ dataAdyt$northatlanticrightwhaleP == 1 & dataAdyt$Yr == 2020 & dataAdyt$Month == 4, 2] )
length( dataAdyt [ dataAdyt$northatlanticrightwhaleP == 1 & dataAdyt$Yr == 2019 & dataAdyt$Month == 4, 2] )
#fin whales
length( dataAdyt [ dataAdyt$finwhaleP == 1 & dataAdyt$Yr == 2020 & dataAdyt$Month == 4, 2] )
length( dataAdyt [ dataAdyt$finwhaleP == 1 & dataAdyt$Yr == 2019 & dataAdyt$Month == 4, 2] )

#SET UP PARAMETERS-- plotting
#-----------------------------------------------------------------------------------------
clrs= colorRampPalette(rev(c("orange","purple")))(100)

# currently no data for December 2020 in the comparison
#-----------------------------------------------------------------------------------------
#ANALYSIS: run as separate sites
#-----------------------------------------------------------------------------------------
# Robust principal components analysis separates a matrix into a low-rank plus sparse component
#a method for the robust separation of a a rectangular (m, n) matrix A into a low-rank component L and a sparse component S:

setwd(sdir)

#FORMAT DATA: TOL variation-- hourly one third octave bands, 
#-----------------------------------------------
df = NULL
for (ii in 1:length(TOLFiles)){
  dftmp = as.data.frame ( fread(input=TOLFiles[ii]))
  deployment = sapply( strsplit(basename(TOLFiles[ii]),"_") , `[`, 3)
  cat(names(dftmp),deployment,"\n")
  df = rbind(df,cbind(dftmp,site,deployment))
}
head(df)

#format the date column 
df$`yyyy-mm-ddTHH:MM:SSZ`
df$DateF     = as.POSIXct( gsub(".000Z", "", gsub("T", " ", df$`yyyy-mm-ddTHH:MM:SSZ`)), tz = "GMT" )
df$DateFday  = as.Date(df$DateF)
df$Hr  = hour(df$DateF)
df$Mth = month(df$DateFday)

#truncate to time period of interest
if (site == "SB02") {
  tmpAP   = timesInterst[timesInterst$Site == site,]
  nv   = df[ (df$DateFday >= as.Date(tmpAP$DateStart,format = "%m/%d/%Y") &
                df$DateFday <= as.Date(tmpAP$DateEnd,format = "%m/%d/%Y") ),]
} else {
  nv   = df
}

#format the data  
as.data.frame(colnames(nv))
idx = grep("^TOL", colnames(nv))
hix  = as.numeric( gsub("TOL_","",names(nv)[idx]) )
Nv   = as.matrix( ( nv[,idx]) ) #dB
hist(as.numeric(Nv), main = "check TOL values")
NvP  = as.matrix(10^(Nv/20))     #pressure units
nvDate = nv$DateF

#RRPCA method for source separation
#-----------------------------------------------------------------------------------------
input = NvP # or Nv (dB units)
lamd = max(input)^-0.5 #default
nvpcaTOL = rrpca(input)
#NOTE: pressure units avoid a negative change in ambient when transients removed

#OUTPUT formatting
#-----------------------------------------------------------------------------------------
#input data 
Am = as.data.frame(input)   

#low rank
Lr = as.data.frame(nvpcaTOL$L) 
colnames(Lr) = hix
LrDB = 10*log10( Lr^2 )  #CHECK: min(LrDB$`63`), no negative values, just values without transients
colnames(LrDB) = hix
# LrDB$DateTimeF = nv$DateF
LrDB =as.matrix(LrDB)

#sparse matrix
Sp = as.data.frame(nvpcaTOL$S) 
colnames(Sp) = hix
# negative and zero values-- does not make sense to convert back to dB
SpDB = 10*log10( (Sp)^2 ) 
colnames(SpDB) = hix
#Sp$DateTimeF = nv$DateF

# visualize results-- LOW RANK (works for all sites)
#-----------------------------------------------------------------------------------------
#1) spectra for each hour
NvMlr   = reshape :: melt(t(LrDB)  )

#only plot specific days to simplify spectra plots
uhrs = unique( NvMlr$X2 )
dateInfo = as.data.frame( cbind((uhrs), (as.character(nvDate) )) )
dateInfo$DateF     = as.POSIXct( gsub(".000Z", "", gsub("T", " ", dateInfo$V2)), tz = "GMT" )
dateInfo$yr  = year(dateInfo$DateF )
dateInfo$mth = month(dateInfo$DateF )
dateInfo$hr  = hour(dateInfo$DateF )
dateInfo$jday  = yday(dateInfo$DateF )
dateInfo$wday  = weekdays(dateInfo$DateF )
dateInfo$day  = as.Date(dateInfo$DateF )
idx    = uhrs[which(dateInfo$mth  == 4)] #only select April days
smth = "April"


#1) spectra-- truncated to specific days- 
tLabel  = paste0( site, "-", as.character(min(nv$DateF)), "-", as.character( max(nv$DateF) ), " (hours = ", length(nv$DateF), ")")

# original matrix
NvMO   = reshape :: melt(t(Nv)  )
NvMOt  = NvMlr[NvMO$X2 %in% idx,]
pO = ggplot(NvMOt, aes(X1, value, group = as.factor(X2)))+ 
  geom_line(alpha=.05)+ scale_x_continuous(trans = "log10")+ 
  labs(title = paste(site, " Original (", smth, ")"))+
  xlab("") + ylab("")+
  theme_minimal()

# low rank matrix-- color by year?


NvMlr   = reshape :: melt(t(LrDB)  )
#add year column in NvMlr based on V2 value-- this takes way too long!!
NvMlr$yr = 0
idx2     = uhrs[which(dateInfo$yr  == 2019)]
for (ii in 1:length(idx2)){
  NvMlr$yr[ NvMlr$X2  == idx2[ii] ] = 2019
}

idx2     = uhrs[which(dateInfo$yr  == 2020)]
for (ii in 1:length(idx2)){
  NvMlr$yr[ NvMlr$X2  == idx2[ii] ] = 2020
}

#only select April dates
NvMlrt  = NvMlr[NvMlr$X2 %in% idx,] 
pL = ggplot(NvMlrt, aes(X1, value, group = as.factor(X2), color = as.factor(yr) ) ) + 
  geom_line(alpha=.05)+
  scale_x_continuous(trans = "log10")+ 
  labs(title = "Low Rank")+ 
  scale_color_manual(values=c('#E69F00','#56B4E9')) + #c("gray",'#E69F00','#56B4E9')
  xlab("") + ylab("")+
  ylim(c(60,114))+
  theme_minimal()+
  theme(legend.position="none") 
# 
#E69F00 = orange (2019)
#56B4E9 = blue   (2020)

# sparse matrix
NvMsp <- reshape :: melt(t(Sp)  )
NvMspt  = NvMsp[NvMsp$X2 %in% idx,]
pS = ggplot(NvMspt, aes(X1, (value), group = as.factor(X2)))+ 
  geom_line(alpha=.07)+ scale_x_continuous(trans = "log10")+ 
  labs(title = "Sparse", caption = tLabel )+
  xlab("Frequency [Hz]") + ylab("")+
   theme_minimal()

grid.arrange(pL, pS)


#SB02-- figure out if off at same windspeed-- yes it is an offset in dB
# as.data.frame(colnames(dataAdyt))
# dataAdyt$BioP2 = rowSums( dataAdyt[, c(50, 53:57)] ) 
# 
# dataAdytt = dataAdyt[ dataAdyt$BioP2 == 0, ]
# #dataAdytt = dataAdyt[ dataAdyt$avgWSPDL > 20, ]
# 
# ggplot(dataAdytt,aes( avgWSPDL, OL_125, color=as.factor(Yr) ) ) +
# geom_point(size = dataAdytt$PropVessel_daily/10) 
# 
# dataAdyttt = dataAdytt[, c(60, 41, 6:15)]
# dataAdytttm <- reshape :: melt(( dataAdyttt)  )

#2) Spectrograms/images (not used)
# oldpar <- par(no.readonly=T)
# par(mfrow=c(2,1),oma=c(4,2,1,1),mar=c(1,4.5,1,1)+.1,mgp=c(3.5,0.75,0),las=1)
# zl=c(min(min(LrDB),min(LrDB)), max(Nv))
# clrs=colorRampPalette(rev(c("orange","purple")))(100)
# 
# image(Nv[dim(Nv)[1]:1,],col=clrs,axes=F,zlim=zl,main = tLabelo)
# axis(side=2,at=1:11/11,labels=o3b$CentHz[seq(from=3,to=33,by=3)])
# 
# image(LrDB[dim(LrDB)[1]:1,],col=clrs,axes=F,zlim=zl,main = tLabel)
# axis(side=2,at=1:11/11,labels=o3b$CentHz[seq(from=3,to=33,by=3)])
# 
# mtext(text="time in hours",side=1,outer=T,las=0,line=2.5)
# mtext(text="one-third octave band (Hz)",side=2,outer=T,las=0)
# par(oldpar)

#DIFFERENCE in LR from "ambient" in each frequency band and each hour (difference matrix)
# how much did ambient change when transients removed?
#-----------------------------------------------------------------------------------------
# Difference in LR to ambient, if positive, transients removed, so ambient (LR) is actually XX lower 
diffAmb = ( Nv - LrDB ) # difference in dB
diffAmb = ( Nv - LrDB ) # energy from the transients... difference in dB

colnames(diffAmb) = hix
diffAmbt = reshape :: melt(t(diffAmb)  )
diffAmb = as.data.frame(( diffAmb))

#what was the average change for each year- in specific frequency 
diffAmb_date    =  cbind(diffAmb, nv$DateF) 
diffAmb_date$Yr = year(diffAmb_date$`nv$DateF`)
#aggregate(diffAmb$`20`,  by=list(diffAmb_date$Yr), mean)
#aggregate(diffAmb$`125`, by=list(diffAmb_date$Yr), min)

#PLOT: what does the difference look like, summarized by frequency, each point is an hourly measure
# ggplot(diffAmbt, aes(value, as.factor(X1), group = X1) )+
#   geom_boxplot()+
#   xlab("Difference in ambient")+
#   ylab("Frequency (1/3 octave bands)")+
#   ggtitle("How much transient sounds add to ambient levels")

#just plot at bars for mean different with error bars
diffAmbtFQ = aggregate(diffAmbt, by=list(diffAmbt$X1), mean)
tmp = aggregate(diffAmbt, by=list(diffAmbt$X1), sd)
diffAmbtFQ$sd = tmp$value

ggplot(diffAmbtFQ, aes(as.factor(X1), value, group = X1) )+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin=value, ymax=value+sd), position=position_dodge(.5), width=.2) +
  ylab("Average Difference Origional Sound Levels and Low Rank")+
  xlab("Frequency (1/3 octave bands)")+
  ggtitle(paste0(site, ": How much transient sounds add to ambient levels"))+
  ylim(c(-.2,6))+
  theme_minimal()
#save out to illustrator and add biological sound bands to the graphic, included in ESOMM

#now want to zoom into April 2020 vs 2019... evidence for more "effect" of transient signals?
#-----------------------------------------------------------------------------------------
diffAmb_date$Mth = month(diffAmb_date$`nv$DateF`)
diffAmb_dateFeb  = diffAmb_date[diffAmb_date$Mth == 4,]

diffAmbtFQ2 = aggregate(diffAmb_dateFeb,  by=list(diffAmb_dateFeb$Yr), mean)
diffAmbtFQ2m = reshape :: melt(diffAmbtFQ2, id.vars = "Yr", 
                               measure.vars = names(diffAmbtFQ2)[2:34])

tmp = aggregate(diffAmb_dateFeb, by=list(diffAmb_dateFeb$Yr), sd)
tmp = reshape :: melt(tmp, id.vars = "Group.1", 
                      measure.vars = names(diffAmbtFQ2)[2:34])
diffAmbtFQ2m$sd = tmp$value

ggplot(diffAmbtFQ2m, aes(as.factor(variable), value, fill = as.factor(Yr)) )+
  geom_bar(stat = "identity",color="black", position=position_dodge(), width=0.65, size=0.3) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), position=position_dodge(.5), width=.2) +
  ylab("Average Difference Origional Sound Levels and Low Rank")+
  xlab("")+
  scale_fill_manual(values=c( "#E69F00","#56B4E9") )+
  ggtitle("")+
  theme_minimal()

#difference in low rank
LrDBt = reshape :: melt(t(LrDB)  )
#add date information
LrDB = as.data.frame(( LrDB))
LrDB_date     =  cbind(LrDB, nv$DateF) 
LrDB_date$Mth = month(LrDB_date$`nv$DateF`)
LrDB_date$Yr  = year(LrDB_date$`nv$DateF`)
LrDB_date_Apr  = LrDB_date[LrDB_date$Mth == 4,]


LrDB_date_Apr2 = aggregate(LrDB_date_Apr,  by=list(LrDB_date_Apr$Yr), mean)
LrDB_date_Apr2m = reshape :: melt(LrDB_date_Apr2, id.vars = "Yr", 
                               measure.vars = names(LrDB_date_Apr2)[2:34])
tmp = aggregate(LrDB_date_Apr, by=list(LrDB_date_Apr$Yr), sd)
tmp = reshape :: melt(tmp, id.vars = "Group.1", 
                      measure.vars = names(diffAmbtFQ2)[2:34])
LrDB_date_Apr2m$sd = tmp$value

ggplot(LrDB_date_Apr2m, aes( ((variable)), value, fill = as.factor(Yr)) )+
  geom_bar(stat = "identity",color="black", position=position_dodge(), width=0.65, size=0.3) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), position=position_dodge(.5), width=.2) +
  ylab("Average Difference Origional Sound Levels and Low Rank")+
  xlab("Low Rank Levels ")+
  scale_fill_manual(values=c( "#E69F00","#56B4E9") )+
  ggtitle("")+
  theme_minimal()
#use this one!
ggplot(LrDB_date_Apr2m, aes(((variable)), value, color = as.factor(Yr)) )+
  geom_point( position = position_dodge(.5) ) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), position=position_dodge(.5), width=.2) +
  ylab("Low Rank")+
  xlab("Frequency 1/3 octave band")+
  scale_color_manual(values=c( "#E69F00","#56B4E9") )+
  ggtitle("")+
  theme_minimal()

x1 = LrDB_date_Apr2m[(LrDB_date_Apr2m$Yr== 2019),]
x2 = LrDB_date_Apr2m[(LrDB_date_Apr2m$Yr== 2020),3]

xd = LrDB_date_Apr2m[(LrDB_date_Apr2m$Yr== 2020),3] -  LrDB_date_Apr2m[(LrDB_date_Apr2m$Yr== 2019),3]
cbind(x1, x2,xd)
mean(xd) ## 2 dB offset??

ggplot(LrDB_date_Apr2m, aes(as.numeric(as.character((variable)) ), value, color = as.factor(Yr)) )+
  geom_point()+
  geom_line()
  geom_ribbon(aes(ymin=value-sd, ymax=value+sd), alpha = .1) +
  ylab("Low Rank")+
  xlab("Frequency 1/3 octave band ")+
  scale_fill_manual(values=c( "#E69F00","#56B4E9") )+
  ggtitle("")+
  theme_minimal()

# just look at sparse matrix
Sp = as.data.frame(( Sp))
Sp_date     =  cbind(Sp, nv$DateF) 
Sp_date$Mth = month(Sp_date$`nv$DateF`)
Sp_date$Yr  = year(Sp_date$`nv$DateF`)
Sp_date_Apr = diffAmb_date[Sp_date$Mth == 4,]

Sp_date_Apr2 = aggregate(Sp_date_Apr,  by=list(LrDB_date_Apr$Yr), mean)
Sp_date_Apr2m = reshape :: melt(Sp_date_Apr2, id.vars = "Yr", 
                                  measure.vars = names(Sp_date_Apr2)[2:34])
tmp = aggregate(Sp_date_Apr, by=list(Sp_date_Apr$Yr), sd)
tmp = reshape :: melt(tmp, id.vars = "Group.1", 
                      measure.vars = names(diffAmbtFQ2)[2:34])
Sp_date_Apr2m$sd = tmp$value


ggplot(Sp_date_Apr2m, aes(((variable)), value, color = as.factor(Yr)) )+
  geom_point( position = position_dodge(.5) ) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), position=position_dodge(.5), width=.2) +
  ylab("Low Rank")+
  xlab("Frequency 1/3 octave band")+
  scale_fill_manual(values=c( "#E69F00","#56B4E9") )+
  ggtitle("")+
  theme_minimal()
#use this one!
ggplot(Sp_date_Apr2m, aes( ((variable)), value, fill = as.factor(Yr)) )+
  geom_bar(stat = "identity",color="black", position=position_dodge(), width=0.65, size=0.3) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), position=position_dodge(.5), width=.2) +
  ylab("Monthly Average Sparse Matrix Values")+
  xlab("Frequency 1/3 octave band")+
  scale_fill_manual(values=c( "#E69F00","#56B4E9") )+
  ggtitle("")+
  theme_minimal()

# what are the features of the low rank matrix? SVD of the LR spectra
#SB02: can we separate wind vs ship dominated spectra?
#-----------------------------------------------------------------------------------------
# SPL in row = hour, columns = frequency
indf = nvpcaTOL$L 

# standardize the matrix, used math of correspondence: 
#https://www.displayr.com/math-correspondence-analysis/
B = as.matrix( scale(indf) )
B[1,1]

# find a new set of variables/features that are uncorrelated and explain as much variance in the data as possible
Lsvd = svd(B) 
length(Lsvd$d)  #singular values, 33 values... for each frequency 1 x 33
dim(Lsvd$u)     #left singular vectors, rows of inputs data (hour)  721 x 30
dim(Lsvd$v)     #right singular vectors, matrix columns (frequency) 33  X 33

#PLOT: summary of the SVD results
par(mfrow = c(1, 3))
image(t(B)[, nrow(B):1], main = "Original Data")
plot(Lsvd$u[, 1], 17521:1, ylab = "Hour index",      xlab = "First left singular vector", pch = 19)
plot(Lsvd$v[, 1],          xlab = "Frequency index", ylab = "First right singular vector", pch = 19)

#variance explained by the columns- frequency in low-rank matrix
prop.table = round_any( prop.table(svd(B)$d^2),.001, floor)
plot(Lsvd$d)
sum(prop.table[1:3]) #variance explained by first 3 components
sum(prop.table[1:2]) #variance explained by first 2 components

#PLOT: How do I know what these components represent, what frequency bands??
dim( data.frame(Lsvd$u[ ,  1:2]) )
par(mfrow = c(1, 3))
plot(Lsvd$u[ ,  1]) #looks like windspeed plot, except at the end with fish
plot(Lsvd$u[ ,  2])
plot(Lsvd$u[ ,  3])

#frequency bands-- 
par(mfrow = c(1, 3))
plot(hix, Lsvd$v[ ,  1]) 
plot(hix, Lsvd$v[ ,  2])
plot(hix, Lsvd$v[ ,  3])
FQ_LRLF = hix[11:18]
FQ_LRHF = hix[24:30]

#PLOT: approximate the original data with outer product of the 1st and 2nd singular vectors
# can see what frequency bands are contributing to the SVD
approx  = with(Lsvd, outer(u[, 1], v[, 1]) )
approx2 = with(Lsvd, outer(u[, 2], v[, 2]) )
approx3 = with(Lsvd, outer(u[, 3], v[, 3]) )
par(mfrow = c(1, 4))
image(t(B)[, nrow(B):1], main = "Original Data")
image(t(approx)[, nrow(approx):1],   main = "Approximated Matrix-1st ")

#ugh did not separate out the fish and wind!!!
image(t(approx2)[, nrow(approx2):1], main = "Approximated Matrix-2nd ")
image(t(approx3)[, nrow(approx2):1], main = "Approximated Matrix-2nd ")
#?? so, we can separate out the low frequency (wind) and high frequency (shrimp),
# and recreate the conditions... and see the magnitude change over time, accurate estimate of ambient in 500






# visualize results-- Sparse
#-----------------------------------------------------------------------------------------
NvMsp <- reshape :: melt(t(SpDB)  )
hist(NvMsp$value)
ggplot(NvMsp, aes(X1, (value), group = as.factor(X2)))+ 
  geom_line(alpha=.05)+ scale_x_continuous(trans = "log10")+ 
  ggtitle("Sparse")+ xlab("")+ theme_minimal()

#?? still do not understand why levels are negative-- 
# how could ambient levels be lower when transients removed?
# copied to PPT "Source Separation"

### Separate by year/month to see differences by frequency
# Is the observed change in ambient from transients?
nvRc = cbind(as.character(nv$DateF), diffAmb)
colnames(nvRc)[1] = "DateF"
nvRc$DateF = as.POSIXct(nvRc$DateF , tz = "GMT") #head(nvRc) #check
nvRc$Yr   = year(nvRc$DateF)
nvRc$Mth  = month(nvRc$DateF)
nvRcMar   = nvRc[nvRc$Mth == 3,]

nvRcMar19 = nvRcMar[nvRcMar$Yr == 2019,]
nvRcMar20 = nvRcMar[nvRcMar$Yr == 2020,]
diffAmb19 <- reshape :: melt(t(nvRcMar19[,2:34]) )
diffAmb20 <- reshape :: melt(t(nvRcMar20[,2:34]) )

ggplot(diffAmb19, aes(value, as.factor(X1), group = X1) )+
  geom_boxplot()+
  xlab("Difference in ambient")+
  ylab("Frequency (1/3 octave bands)")+
  ggtitle("March 2019")
ggplot(diffAmb20, aes(value, as.factor(X1), group = X1) )+
  geom_boxplot()+
  xlab("Difference in ambient (pressure units)")+
  ylab("Frequency (1/3 octave bands)")+
  ggtitle("March 2020")

# Combine vessel detections with difference...



#PLOTS
#-----------------------------------------------------------------------------------------
#low-rank spectra-- spectra for each hour-- too much data.. all data
NvM <- reshape :: melt(t(nvpcaTOL$L)  )
p2 = ggplot(NvM, aes(X1, value, group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("Low rank")+ xlab("")+ theme_minimal()
#sparse spectra-- spectra for each hour
NvMs <- reshape :: melt(t(nvpcaTOL$S)  )
p3=ggplot(NvMs, aes(X1, value,group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("Sparse")+ xlab("")+ theme_minimal()
#low rank + sparse gets back to origional data
NvMs <- reshape :: melt(t(nvpcaTOL$S) +t(nvpcaTOL$L)  )
p4=ggplot(NvMs, aes(X1, value,group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("L+S")+ xlab("")+ theme_minimal()
#spectra from original data
NvM1 <- reshape :: melt(t(Nv)  )
NvM1$X1 = as.numeric( gsub("TOL_","",NvM1$X1) )
p1=ggplot(NvM1, aes(log10(X1), value, group = as.factor(X2)) ) + geom_line(alpha=.1) + xlab("")+ ggtitle("Original Data")+ theme_minimal()
#grid.arrange(p1,p2,p3,p4)
#grid.arrange(p2,p3,nrow=1)
grid.arrange(p1,p2,p3,ncol=1)

#-----------------------------------------------------------------------------------------
#compare 2019 vs 2020 spectra-- 
#-----------------------------------------------------------------------------------------
nv$Yr = year(nv$DateF)
idx19 = which(nv$Yr =="2019")
idx20 = which(nv$Yr =="2020")

idxMarch = which( nv$Mth == 4 )
nvMarch = nv[idxMarch,]
idxM19 = which(nvMarch$Yr =="2019")
idxM20 = which(nvMarch$Yr =="2020")

#LOW RANK SPECTRA-- separate by year
#-----------------------------------------------------------------------------------------
#select only specific rows--- year and months
Lrpc19 = as.data.frame( nvpcaTOL$L[idx19,] )#2019 data
Lrpc20 = as.data.frame( nvpcaTOL$L[idx20,] )#2020 data
dim(Lrpc19)+ dim(Lrpc20)
NvM19 <- reshape :: melt(t(Lrpc19)  )
Freq= rep(hix, length(NvM19$X1)/length(hix))
NvM19 = cbind(NvM19,(Freq))
p2019 = ggplot(NvM19, aes(Freq, value, group = as.factor(X2))) +
  geom_line(alpha=.1) + scale_x_continuous(trans = "log10") +
  ylim(c(60,110))+
  ggtitle("2019- Low rank")+ xlab("")+ theme_minimal()
#2020
NvM20 <- reshape :: melt(t(Lrpc20)  )
Freq  = rep(hix, length(NvM20$X1)/length(hix))
NvM20 = cbind(NvM20,(Freq))
p2020 = ggplot(NvM20, aes(Freq, value, group = as.factor(X2))) +
  geom_line(alpha=.1) + scale_x_continuous(trans = "log10") +
  ylim(c(60,110))+
  ggtitle("2020- Low rank")+ xlab("")+ theme_minimal()

#quantiles for the low rank spectra
quants=c(.05,.25,.5,.75,.95)
dim(Lrpc20)
q2020 = apply(Lrpc20,2,quantile,probs=quants)
colnames(q2020) = hix
q2019 = as.data.frame( apply(Lrpc19,2,quantile,probs=quants))
colnames(q2019) = hix
#how do I plot these....
tst = reshape :: melt(t(q2019)  )
tst$YR = 2019
tst2 = reshape :: melt(t(q2020)  )
tst2$YR = 2020
ggplot(tst,aes(X1,value,color=X2))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans = "log10") +
    ggtitle("2020- Low rank")+ xlab("")+ theme_minimal()
#sparse spectra-- spectra for each hour
tst3 = rbind(tst,tst2)
ggplot(tst3,aes(X1,value,color=X2))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans = "log10") +
  ggtitle("Low rank compare")+ xlab("")+ theme_minimal()+
  facet_wrap(~YR)
# 3 dB higher in 2020.... in low rank


#SPARse RANK SPECTRA-- separate by year
#-----------------------------------------------------------------------------------------
#select only specific rows--- year and months
Srpc19 = as.data.frame( nvpcaTOL$S[idx19,] )#2019 data
Srpc20 = as.data.frame( nvpcaTOL$S[idx20,] )#2020 data

NvM19 = reshape :: melt(t(Srpc19)  )
Freq= rep(hix, length(NvM19$X1)/length(hix))
NvM19 = cbind(NvM19,(Freq))
ggplot(NvM19, aes(Freq, value, group = as.factor(X2))) +
  geom_line(alpha=.1) + scale_x_continuous(trans = "log10") +
  ylim(c(-10,40))+
  ggtitle("2019- Sparse")+ xlab("")+ theme_minimal()

#2020
NvM20 <- reshape :: melt(t(Srpc20)  )
Freq  = rep(hix, length(NvM20$X1)/length(hix))
NvM20 = cbind(NvM20,(Freq))
ggplot(NvM20, aes(Freq, value, group = as.factor(X2))) +
  geom_line(alpha=.1) + scale_x_continuous(trans = "log10") +
  ylim(c(-10,40))+
  ggtitle("2020- Sparse")+ xlab("")+ theme_minimal()

#quantiles for the sparse
q2020 = apply(Srpc20,2,quantile,probs=quants)
colnames(q2020) = hix
q2019 = as.data.frame( apply(Srpc19,2,quantile,probs=quants))
colnames(q2019) = hix
tst = reshape :: melt(t(q2019)  )
tst$YR = 2019
tst2 = reshape :: melt(t(q2020)  )
tst2$YR = 2020
tst3 = rbind(tst,tst2)
ggplot(tst3,aes(X1,value,color=X2))+
  geom_point()+
  geom_line()+
  scale_x_continuous(trans = "log10") +
  ggtitle("Sparse compare")+ xlab("")+ theme_minimal()+
  facet_wrap(~YR)
# infrequeny 95th percentile events in low frequency 



# NOT WORKING YET......
#what do these low-rank spectra represent- looks like different wind speeds
#LOW RANK is really same spectra with sparse variation removed.... how can I collapse these? 
c1  = as.data.frame ( cbind( nvpcaTOL$L, (cdf$avgWSPD) ))
colnames(c1) = c(hix,"wspd")
wspdR = do.call(rbind, replicate(ncol(c1)-1, c1$wspd, simplify=FALSE)) #replicate windspeed for each frequency band- because same!
wspdRm <- reshape :: melt((wspdR)  )
colnames(wspdRm)=c("FQ","TIMESTEP","WSPD")
NvMw = cbind(NvM, wspdRm[,3])
colnames(NvMw)=c("FQ","TIMESTEP","SPL", "WSPD")


#PLOT: Low-rank spectra colored by windspeed
ggplot(NvMw, aes(FQ, SPL, group = as.factor(TIMESTEP), colour=as.factor( round_any(WSPD, 5, f = ceiling) ) ) ) + 
  geom_line(alpha = .1) + 
  ggtitle("Low rank") + 
  scale_color_discrete(name = "Windspeed m/s") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30]))
#PLOT: the first three Low Rank to visualize different spectra-- not re-ordered in any way
NvMw3 = NvMw[NvMw$TIMESTEP <=5,]
#ggplot(NvMw3, aes(FQ, SPL, group = as.factor(TIMESTEP)) ) + 
  #geom_line(alpha = 1) 

#NEXT: create average "low rank" spectra for each windspeed and plot,
summary(c1$wspd )
c1$wspdR = round_any(c1$wspd, 5, f = floor) #add average windspeed
n = ncol(nvpcaTOL$L) # number of frequency bands
st = nrow(nvpcaTOL$L) # number of hours, timesteps
quants = c(.25,.5,.75)
quants = c(.5)
w5  = c1[ c1$wspdR <5,] #pick windspeed
w5p = apply( w5[1:30] , 2 , quantile , probs = quants , na.rm = TRUE )
w5p = reshape::melt(w5p)
n5 = nrow(w5)
w5p$FQ <- row.names(w5p)

w10  = c1[ c1$wspdR >=5 & c1$wspdR <10,]
w10p = apply( w10[1:30] , 2 , quantile , probs = quants , na.rm = TRUE )
w10p = reshape::melt(w10p)
n10 = nrow(w10)
w10p$FQ <- row.names(w10p)

w15  = c1[ c1$wspdR >=10 & c1$wspdR <15,]
w15p = apply( w15[1:30] , 2 , quantile , probs = quants , na.rm = TRUE )
w15p = reshape::melt(w15p)
n15 = nrow(w15)
w15p$FQ <- row.names(w15p)

w20  = c1[ c1$wspdR >=15 & c1$wspdR <21,]
w20p = apply( w20[1:30] , 2 , quantile , probs = quants , na.rm = TRUE )
w20p = reshape::melt(w20p)
n20 = nrow(w20)
w20p$FQ <- row.names(w20p)

ggplot()+
  geom_line(data = w5p,  aes(log10(as.numeric(as.character((FQ)))),value) ,color = "#CD5C5C" ) +    #red
  geom_line(data = w10p, aes(log10(as.numeric(as.character((FQ)))),value) ,color = "#9ACD32"  ) +   #green-light
  geom_line(data = w15p, aes(log10(as.numeric(as.character((FQ)))),value) ,color = "#2E8B57"  ) +   #green-forest
  geom_line(data = w20p, aes(log10(as.numeric(as.character((FQ)))),value) ,color = "dodgerblue"  )+ #blue
  xlab("Frequency")+
  ylab("Low-rank SPL")+
  theme_minimal()
#copy this to PPT and label

#sparse what do the sparse represent? 
#---------------------------------------------------
#can I match with presence of biological sounds.. hard to do all at once
#plot for each source, but does not deal with possible overlap in sources present

#create a few summary rows for detection variables....
unique( cdetft$midshipmanDet)
unique( cdetft$bocaccioDet) #number of calls present
unique( cdetft$FishCallsDet) #same as bocaccioDet
unique( cdetft$bluewhaleDet)
unique( cdetft$explosionsDet)
unique( cdetft$VesselDet_count)

cdetft$Fish   = rowSums(cbind( cdetft$midshipmanDet,cdetft$bocaccioDet ) )
cdetft$FishP  = as.numeric( lapply(cdetft$Fish, function(x) ifelse(x > 1, 1, x)) )
cdetft$VesP   = as.numeric( lapply(cdetft$VesselDet_count, function(x) ifelse(x > 1, 1, x)) )
cdetft$ALL    = rowSums(cbind( cdetft$FishP, cdetft$bluewhaleDet, cdetft$explosionsDet,cdetft$VesselDet_count) )
#only MS
cdetft$ALL1    = rowSums(cbind( cdetft$bocaccioDet, cdetft$bluewhaleDet, cdetft$explosionsDet, cdetft$VesselDet_count) )
cdetft$onlyMS = 0
idx =  which (cdetft$ALL1 == 0 & cdetft$midshipmanDet >=1 )
cdetft$onlyMS[idx] = 1
#only VES
cdetft$ALL2    = rowSums(cbind( cdetft$FishP, cdetft$bluewhaleDet,cdetft$explosionsDet,cdetft$midshipmanDet) )
cdetft$onlyVES = 0
idx =  which (cdetft$ALL2 == 0 & cdetft$VesP >=1 )
cdetft$onlyVES[idx] = 1

#create appended data frame
c1S = as.data.frame ( cbind( nvpcaTOL$S, (cdetft) ))
colnames(c1S)[1:30] = hix

#PLOT: sparse matrix with MIDSHIPMAN presence
detR = do.call(rbind, replicate(ncol(nvpcaTOL$S), c1S$onlyMS, simplify=FALSE)) #only samples with MS
samps = nrow(c1S[c1S$onlyMS>0,])
#detR = do.call(rbind, replicate(ncol(nvpcaTOL$S), c1S$midshipmanDet, simplify=FALSE))
detRm <- reshape :: melt((detR)  )
colnames(detRm)=c("FQ","TIMESTEP","DET")
NvM <- reshape :: melt(t(nvpcaTOL$S)  )
NvLdet = cbind(NvM, detRm[,3])
colnames(NvLdet)=c("FQ","TIMESTEP","SPL", "DET")
#unique ( as.numeric( as.character(NvLdet$DET) ) )
NvLdett = NvLdet[as.numeric( as.character(NvLdet$DET)) > 0,]
pMS = ggplot(NvLdett, aes(FQ, SPL, group = as.factor(TIMESTEP), colour=as.factor(DET ) ) ) + 
  geom_line(alpha = .5) + 
  geom_point(size=.5) +
  ggtitle(paste0("Sparse- midshipman only (N=",samps,")")) + 
  scale_color_discrete(name = "Detections") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30])) +
  theme(legend.position = "none")
#PLOT: sparse matrix with BLUE present
samps = nrow(c1S[c1S$bluewhaleDet>0,])
detR = do.call(rbind, replicate(ncol(nvpcaTOL$S), c1S$bluewhaleDet, simplify=FALSE))
detRm <- reshape :: melt((detR)  )
colnames(detRm)=c("FQ","TIMESTEP","DET")
NvM <- reshape :: melt(t(nvpcaTOL$S)  )
NvLdet = cbind(NvM, detRm[,3])
colnames(NvLdet)=c("FQ","TIMESTEP","SPL", "DET")
NvLdett = NvLdet[as.numeric( as.character(NvLdet$DET)) >= 1,]
pBW = ggplot(NvLdett, aes(FQ, SPL, group = as.factor(TIMESTEP), colour=as.factor(DET ) ) ) + 
  geom_line(alpha = .5) + 
  geom_point(size=.5) +
  ggtitle(paste0("Sparse- blue whale (N=",samps,")")) + 
  scale_color_discrete(name = "Detections") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30]))+
  theme(legend.position = "none")

#PLOT: sparse matrix with VESSEL present, only
samps = nrow(c1S[c1S$VesP>0,])
detR = do.call(rbind, replicate(ncol(nvpcaTOL$S), c1S$VesP, simplify=FALSE))
detRm <- reshape :: melt((detR)  )
colnames(detRm)=c("FQ","TIMESTEP","DET")
NvM <- reshape :: melt(t(nvpcaTOL$S)  )
NvLdet = cbind(NvM, detRm[,3])
colnames(NvLdet)=c("FQ","TIMESTEP","SPL", "DET")
NvLdett = NvLdet[as.numeric( as.character(NvLdet$DET)) >= 1,]
pVE = ggplot(NvLdett, aes(FQ, SPL, group = as.factor(TIMESTEP), colour=as.factor(DET ) ) ) + 
  geom_line(alpha = .5) + 
  geom_point(size=.5) +
  ggtitle(paste0("Sparse- vessel (N=",samps,")")) + 
  scale_color_discrete(name = "Detections") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30]))+
  theme(legend.position = "none")

#PLOT: sparse matric with VESSEL present
detR = do.call(rbind, replicate(ncol(nvpcaTOL$S), c1S$explosionsDet, simplify=FALSE))
detRm <- reshape :: melt((detR)  )
colnames(detRm)=c("FQ","TIMESTEP","DET")
NvM <- reshape :: melt(t(nvpcaTOL$S)  )
NvLdet = cbind(NvM, detRm[,3])
colnames(NvLdet)=c("FQ","TIMESTEP","SPL", "DET")
NvLdett = NvLdet[as.numeric( as.character(NvLdet$DET)) >= 1,]
pEX = ggplot(NvLdett, aes(FQ, SPL, group = as.factor(TIMESTEP), colour=as.factor(DET ) ) ) + 
  geom_line(alpha = .5) + 
  geom_point(size=.5) +
  ggtitle("Sparse- explosions present") + 
  scale_color_discrete(name = "Detections") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30]))+
  theme(legend.position = "none")


grid.arrange(pMS,pBW,pVE,pEX, nrow = 2, ncol = 2)
#copy this to PPT and label


#----------------------------------------------------------------------------------------------
#what does this all mean??!!: 
# now my spectra are the low-rank background-- can I classify these?
# can I use SVD to determine spectra that captures most of the variation?
#LOW RANK == AMBIENT interpretations
#----------------------------------------------------------------------------------------------
#SVD tutorials
# https://www.displayr.com/singular-value-decomposition-in-r/
#----------------------------------------------------------------------------------------------
# input matrix-- want to find spectra that represents data
# SPL in row = hour, columns = frequency
indf = nvpcaTOL$L 
dim(indf) 

#standardize matrix, used math of correspondance: #https://www.displayr.com/math-correspondence-analysis/
# and Dan's code from NOAA data (S1_mergeRED_v7.m)
X = as.matrix(indf)
M = dim(X)[1]
N = dim(X)[2]
uX = colMeans(X)
sX = apply(X, 2, sd)
MusX = 1; NusX = length(sX); # stupid dim(), for repmat
B = (X - matrix(t(matrix(uX,MusX,NusX*1)),MusX*M,NusX*1,byrow=T)) / matrix(t(matrix(sX,MusX,NusX*1)),MusX*M,NusX*1,byrow=T) 
dim(B) 
#? can the scale funtion just do this? YES!!!
Bs = as.matrix( scale(indf) )
Bs[1,1]
B[1,1]

# https://bookdown.org/rdpeng/exdata/dimension-reduction.html
image(1:30, 1:721, t(indf)[, nrow(indf):1]) #low rank matrix
image(1:30, 1:721, t(B)[, nrow(B):1])       #scaled low rank matrix

# hierarchical clustering algorithm, 
# NOTE: does not account for just mean shifted data, and multiple patterns layered on top of each
heatmap(B)
heatmap(indf) # clusters of data in frequency (wind, low frequency crap, snapping shrimp, and clusters of data in time)
hix #1/3 otave bands for reference to graphic

hh = dist(indf) %>% hclust
dataMatrixOrdered <- indf[hh$order, ]
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(rowMeans(dataMatrixOrdered), 721:1, xlab = "Row Mean", ylab = "Row", pch = 19)
plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column Mean", pch = 19)
hix[11:18]

# multivariate observations (hours) with m features (FQ)
# find a new set of variables/features that are uncorrelated and explain as much variance in the data as possible
Lsvd = svd(B) 
length(Lsvd$d) #singular values, 30 values... for each frequency 1 x 30
dim(Lsvd$u) #left dingular vectors, rows of inputs data (hour)  721 x 30
dim(Lsvd$v) #righg dingular vectors, matrix columns (frequency) 30  X 30

#PLOT: summary of the SVD
par(mfrow = c(1, 3))
image(t(B)[, nrow(B):1], main = "Original Data")
plot(Lsvd$u[, 1], 721:1, ylab = "Row", xlab = "First left singular vector", pch = 19)
plot(Lsvd$v[, 1],        xlab = "Column", ylab = "First right singular vector", pch = 19)
#variance explained by the columns- frequency in low-rank matrix
prop.table = round_any( prop.table(svd(B)$d^2),.001, floor)
plot(Lsvd$d)
sum(prop.table[1:3]) #variance explained by first 3 components
sum(prop.table[1:2]) #variance explained by first 3 components

#PLOT: How do I know what these components represent, what frequency bands??
dim( data.frame(Lsvd$u[ ,  1:2]) )
par(mfrow = c(1, 3))
plot(Lsvd$u[ ,  1]) #looks like windspeed plot, except at the end with fish
plot(Lsvd$u[ ,  2])
plot(Lsvd$u[ ,  3])
#frequency bands
par(mfrow = c(1, 3))
plot(Lsvd$v[ ,  1]) 
plot(Lsvd$v[ ,  2])
plot(Lsvd$v[ ,  3])
FQ_LRLF = hix[11:18]
FQ_LRHF = hix[24:30]

#PLOT: approximate the original data with outer product of the 1st and 2nd singular vectors
# can see what frequency bands are contributing to the SVD
approx  = with(Lsvd, outer(u[, 1], v[, 1]) )
approx2 = with(Lsvd, outer(u[, 2], v[, 2]) )
approx3 = with(Lsvd, outer(u[, 3], v[, 3]) )
par(mfrow = c(1, 4))
image(t(B)[, nrow(B):1], main = "Original Data")
image(t(approx)[, nrow(approx):1],   main = "Approximated Matrix-1st (wind,fish)")
#ugh did not separate out the fish and wind!!!
image(t(approx2)[, nrow(approx2):1], main = "Approximated Matrix-2nd (shrimp)")
image(t(approx3)[, nrow(approx2):1], main = "Approximated Matrix-2nd (low frequency crap)")
#?? so, we can separate out the low frequency (wind) and high frequency (shrimp),
# and recreate the conditions... and see the magnitude change over time, accurate estimate of ambient in 500

#INTERPRETATIONS PLOTS...
#PLOT: SPL and LR-SVD1-- and add running average
plot(Lsvd$u[ ,  1]) 
plot(Nv[ ,  1]) 
plot(indf[ ,  1]) 
length(Lsvd$u[ ,  1])
length(cdetft$DateTime)
#add timeStamp, windspeed, SPL in same band, to LR-SVDLF
FQ_LRLF = hix[11:18]
SPL_LRLF = 10*log10( rowMeans(10^(indf[,11:18]/10)) )
SPL_LR = 10*log10( rowMeans(10^(Nv[,11:18]/10)) )
# CHECK: SPL_LRLF[1]  indf[1,11:18]

dfLL = cbind(cdetft$DateTime, cdetft$WSPDms_mean,  SPL_LR, SPL_LRLF,  as.data.frame( Lsvd$u[ ,  1]) )
colnames(dfLL) = c("DateTime","WSPD", "SPL", "LRSPL", "LRSVD1")
dfLL$Diff = dfLL$LRSPL -dfLL$SPL 
dfLL$DateTimeF = as.POSIXct(dfLL$DateTime, tz="GMT")
dfLL$pos = dfLL$Diff < 0

par(mfrow = c(3, 1))
plot(dfLL$DateTimeF, dfLL$SPL)
plot(dfLL$DateTimeF, dfLL$LRSPL)
plot(dfLL$DateTimeF, dfLL$LRSVD1)
plot(dfLL$DateTimeF,dfLL$Diff)

pSVD = ggplot(dfLL,aes(DateTimeF,LRSVD1 ) )+
  #geom_point(alpha = .2)+
  geom_line()
pSPL = ggplot(dfLL,aes(DateTimeF,SPL))+
  #geom_point(alpha = .2)+
  geom_line()
pLRSPL = ggplot(dfLL,aes(DateTimeF,LRSPL))+
  #geom_point(alpha = .2)+
  geom_line()
pWSPD = ggplot(dfLL,aes(DateTimeF,WSPD))+
  #geom_point(alpha = .2)+
  geom_line()
pDiff= ggplot(dfLL,aes(DateTimeF,Diff))+
  #geom_point(alpha = .2)+
  geom_line() +
  xlab("")+
  
  theme_minimal()
grid.arrange(pSPL,pLRSPL,pDiff,nrow = 3,ncol=1)

pComp = ggplot(dfLL,aes(x = DateTimeF ) )+
  geom_line(aes(y =SPL)   ,color="black") + 
  geom_line(aes(y =LRSPL) ,color="red"  ) +
  xlab("")+
  ylab ("SPL [250-1250 Hz]") +
  theme_minimal()
grid.arrange(pComp,pDiff,nrow = 2,ncol=1)

#plot positive negative values diferent colr
pDiff2= ggplot()+
  geom_bar(data=dfLL, aes(x=DateTimeF,y=Diff,fill=pos),size=1.55,stat="identity",position='identity') +
  xlab("")+
  theme_minimal()+
  ylim(c(-2,2)) +
  theme(legend.position = "none")
  
grid.arrange(pComp,pDiff2,nrow = 2,ncol=1)


par(mar = c(5, 5, 3, 5))
plot(dfLL$DateTimeF, dfLL$SPL, type ="l", ylab = "SPL [250-1250 Hz]",xlab = "",
     col = "black")
par(new = TRUE)
plot(dfLL$DateTimeF, dfLL$LRSVD1, type = "l", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",col = "red")
axis(side = 4)
mtext("LR-SVDlow", side = 4, line = 3)
legend("topleft", c("SPL", "LR-SVDlow"),
       col = c("black", "red"), lty = c(1, 1))

par(mar = c(5, 5, 3, 5))
plot(dfLL$DateTimeF, dfLL$LRSVD1, type ="l", ylab = "LR-SVDlow",xlab = "",
     col = "red")
par(new = TRUE)
plot(dfLL$DateTimeF, dfLL$WSPD, type = "l", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",col = "blue")
axis(side = 4)
mtext("WSPD", side = 4, line = 3)
legend("topleft", c("LR-SVDlow", "WSPD"),
       col = c("red", "blue"), lty = c(1, 1))

par(mar = c(5, 5, 3, 5))
plot(dfLL$DateTimeF, dfLL$SPL, type ="l", ylab = "SPL [250-1250 Hz]",xlab = "",
     col = "blue")
par(new = TRUE)
plot(dfLL$DateTimeF, dfLL$WSPD, type = "l", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",col = "red", lty = 2)
axis(side = 4)
mtext("WSPD", side = 4, line = 3)
legend("topleft", c("SPL", "WSPD"),
       col = c("blue", "red"), lty = c(1, 2))

#find breakpoints or transitions in conditions in the month
# https://lindeloev.github.io/mcp/articles/packages.html
library(bcp)
fit_bcp = bcp(dfLL$LRSVD1, d = 1000)
plot(fit_bcp)
length( which(fit_bcp$posterior.prob >.8) )
idxbp = which(fit_bcp$posterior.prob >.8)

par(mar = c(5, 5, 3, 5))
plot(dfLL$DateTimeF, dfLL$LRSVD1, type ="l", ylab = "LR-SVDlow",xlab = "",
     col = "red")
abline(v = dfLL$DateTimeF[idxbp],col="black",lty = 2)

#----------------------------------------------------------------------------
# sparse MATRIX: hierarchical clustering algorithm, sparse matrix 
#----------------------------------------------------------------------------
# input matrix-- want to find spectra that represents data
# SPL in row = hour, columns = frequency
indf = nvpcaTOL$S 
dim(indf) 

#standardize matrix, used math of correspondance: #https://www.displayr.com/math-correspondence-analysis/
# and Dan's code from NOAA data (S1_mergeRED_v7.m)
X = as.matrix(indf)
M = dim(X)[1]
N = dim(X)[2]
uX = colMeans(X)
sX = apply(X, 2, sd)
MusX = 1; NusX = length(sX); # stupid dim(), for repmat
B = (X - matrix(t(matrix(uX,MusX,NusX*1)),MusX*M,NusX*1,byrow=T)) / matrix(t(matrix(sX,MusX,NusX*1)),MusX*M,NusX*1,byrow=T) 
dim(B) 
#? can the scale funtion just do this? YES!!!
Bs = as.matrix( scale(indf) )
Bs[1,1]
B[1,1]
# https://bookdown.org/rdpeng/exdata/dimension-reduction.html
image(1:30, 1:721, t(indf)[, nrow(indf):1]) #lsparse matrix
image(1:30, 1:721, t(B)[, nrow(B):1])       #scaled sparse matrix

#CLUSTER: hierarchical clustering algorithm, 
# NOTE: does not account for just mean shifted data, and multiple patterns layered on top of each
heatmap(B)
heatmap(indf) # cluters of data in frequency (wind, low frequency crap, snapping shrimp, and clusters of data in time)
hh = dist(indf) %>% hclust
dataMatrixOrdered <- indf[hh$order, ]
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(rowMeans(dataMatrixOrdered), 721:1, , xlab = "Row Mean", ylab = "Row", pch = 19)
plot(colMeans(dataMatrixOrdered), xlab = "Column", ylab = "Column Mean", pch = 19)

#SVD:  multivariate observations (hours) with m features (FQ)
# find a new set of variables/features that are uncorrelated and explain as much variance in the data as possible
Ssvd = svd(B) 

#PLOT: summary of the SVD
par(mfrow = c(1, 3))
image(t(B)[, nrow(B):1], main = "Original Data")
plot(Ssvd$u[, 1], 721:1, , ylab = "Row", xlab = "First left singular vector", pch = 19)
plot(Ssvd$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)
#variance explained by the columns- frequency in low-rank matrix
prop.table = round_any( prop.table(svd(B)$d^2),.001, floor)
plot(Ssvd$d)
sum(prop.table[1:10]) #variance explained by first 10 components
sum(prop.table[1:3]) #variance explained by first 3 components

#PLOT: How do I know what these components represent, what frequency bands??
dim( data.frame(Ssvd$u[ ,  1:2]) )
par(mfrow = c(1, 3))
plot(Ssvd$u[ ,  1]) #looks like windspeed plot, except at the end with fish
plot(Ssvd$u[ ,  2])
plot(Ssvd$u[ ,  3])
#frequency bands
par(mfrow = c(1, 3))
plot(Ssvd$v[ ,  1]) 
plot(Ssvd$v[ ,  2])
plot(Ssvd$v[ ,  3])
approx  = with(Ssvd, outer(u[, 1], v[, 1]) )
approx2 = with(Ssvd, outer(u[, 2], v[, 2]) )
approx3 = with(Ssvd, outer(u[, 3], v[, 3]) )
par(mfrow = c(1, 4))
image(t(B)[, nrow(B):1], main = "Original Data")
image(t(approx)[, nrow(approx):1],   main = "Approximated Matrix-1st")
#ugh did not separate out the fish and wind!!!
image(t(approx2)[, nrow(approx2):1], main = "Approximated Matrix-2nd (fish??)")
image(t(approx3)[, nrow(approx2):1], main = "Approximated Matrix-2nd (??)")
#?? so, we can separate out the low frequency (wind) and high frequency (shrimp),
# and recreate the conditions... and see the magnitude change over time, accurate estimate of ambient in 500


#----------------------------------------------------------------------------
# NOW WHATT!!! Thinking.... I am not sure how this is useful, do you compare the input spectra to the re-created SVD options
# and look at differences to understand what is explaining the variatiton. 
# in other words, for each time step (30 days??) calculate the difference from these "typical conditions"
#----------------------------------------------------------------------------
# only in the frequency bands the SVD represents?
svd1.compare = abs (B[,6:17] - approx[,6:17])
sum( svd1.compare )
image(t(svd1.compare)[, nrow(svd1.compare):1], main = "")
# or for each 720 hours, compare SVD to multiple known SVDs? 
# NO!!! because this would mean the entire period would need to look the same

# compare each time-step spectra (hour) to v? but then magnitude is unknown? 
# correlation, not difference https://arxiv.org/ftp/arxiv/papers/1211/1211.7102.pdf
# Step through data to correlate with known SVDs to see if changes occur. hourly steps??
#thinking about CT scans to detect changes in the slices.... maybe daily SVDs would be better for this
# or do you just keep re-calculating and seeing how he SVDs change??

#re-created input matrix with multiplication
sum(Lsvd$d * Lsvd$u[1,]* Lsvd$v[1,])
B[1,1]
Lsvd$u %*% diag(Lsvd$d) %*% t(Lsvd$v)


#----------------------------------------------------------------------------
#EXTRA: Math of correspondance: #https://www.displayr.com/math-correspondence-analysis
#----------------------------------------------------------------------------
n = sum(nvpcaTOL$L)
P = nvpcaTOL$L/n #propotion of total
column.masses = colSums(P)
row.masses    = rowSums(P)
E = row.masses %o% column.masses
R = P - E # residuals
I = R / E # indexed residuals, further the value from table, the large the observed proportion to expected

#sparse
Ssvd = svd(nvpcaTOL$S) 
plot(Ssvd$d)
round_any( prop.table(svd(nvpcaTOL$S)$d^2),.001, floor)






#NEXT: PSD variation-- hourly Hz bins
#-----------------------------------------------
df   = as.data.frame ( fread(input=PSDFiles[2]))
stmp = sapply( strsplit(basename(PSDFiles[2]),"_") , `[`, 2)
#format the data column: 
df$yyyy_mm_ddTHH_MM_SSZ[1:10]
df$DateF     = as.POSIXct( gsub(".000Z", "", gsub("T", " ", df$yyyy_mm_ddTHH_MM_SSZ)), tz = "GMT" )
df$DateFday  = as.Date(df$DateF)
#truncate to time period of interest
ttmp = timesInterst[timesInterst$Site == stmp,]
nv   = df[ (df$DateFday >= as.Date(ttmp$DateStart,format = "%m/%d/%Y") & df$DateFday  <= as.Date(ttmp$DateEnd,format = "%m/%d/%Y") ),]
hix  = as.numeric( gsub("TOL_","",names(nv)[2:31]) )
Nv   = as.matrix( nv[,2:31] )





#KURT CODE:
#-----------------------------------------------
#visualize-- input
rpcaplot <- function(Nv,L,S,kolors){ # packaged code to avoid repetition below, produces a 3 panel plot
  oldpar <- par(no.readonly=T)
  
  par(mfrow=c(3,1),oma=c(4,2,0,0),mar=c(0,4.5,0,0)+.1,mgp=c(3.5,0.75,0),las=1)
  zl <- c(min(min(L),min(S)),max(Nv))
  image(Nv[dim(Nv)[1]:1,],col=kolors,axes=F,zlim=zl)
  axis(side=2,at=1:10/11,labels=o3b$CentHz[seq(from=3,to=30,by=3)])
  
  image(L[dim(L)[1]:1,],col=kolors,axes=F,zlim=zl)
  axis(side=2,at=1:10/11,labels=o3b$CentHz[seq(from=3,to=30,by=3)])
  
  image(S[dim(S)[1]:1,],col=kolors,axes=F,zlim=zl)
  axis(side=2,at=1:10/11,labels=o3b$CentHz[seq(from=3,to=30,by=3)])
  axis(side=1,at=seq(from=300,to=3300,by=300)/3600,outer=T,labels=seq(from=5,to=55,by=5))
  
  mtext(text="time in minutes",side=1,outer=T,las=0,line=2.5)
  mtext(text="one-third octave band (Hz)",side=2,outer=T,las=0)
  
  par(oldpar)
}
rpcaplot(Nv = as.matrix(nv[,hix,with=F]), L=nvpca$L, S=nvpca$S, clrs)
# RPCA() using decibel data -- does not separate as well
system.time(nvpca <- rpca(as.matrix(nv[,hix,with=F]),max.iter=20000))
rpcaplot(Nv=as.matrix(nv[,hix,with=F]),L=nvpca$L,S=nvpca$S,clrs)
# RPCA() using energy data
nvp2 <- rpca(as.matrix(10^(nv[,hix,with=F]/10)),max.iter=20000)

#dealing with negative S values
minS <- 1e-2 # hard coded after reviewing the histogram of positive values for this file
nix <- nvp2$S < minS
nvp2$L[nix] <- nvp2$L[nix]-(minS-nvp2$S[nix])
nvp2$S[nix] <- minS

#dealing with negative L values
minL <- 1 # hard coded after reviewing the histogram of positive values for this file
nix <- nvp2$L < minL
nvp2$S[nix] <- nvp2$S[nix]-(minL-nvp2$L[nix])
nvp2$L[nix] <- minL
rpcaplot(Nv=as.matrix(nv[,hix,with=F]),L=10*log10(nvp2$L),S=10*log10(nvp2$S),clrs)

if(0){
  hist(log10(nvp2$L[0<nvp2$L]))
  hist(log10(nvp2$S[0<nvp2$S]))
  plot(-nvp2$S[0>nvp2$S],nvp2$L[0>nvp2$S],type="p",log="xy",cex=0.1)
  abline(0,1)
  plot(-nvp2$L[0>nvp2$L],nvp2$S[0>nvp2$L],type="p",log="xy")
  abline(0,1)
}

