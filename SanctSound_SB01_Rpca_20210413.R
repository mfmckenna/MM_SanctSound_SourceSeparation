#finding low rank and sparse features in hourly, third octave band spectra in SanctSound project

#use SB01 and CI01 to run rpca (robust PCA):
# low-rank spectra- how does this relate to ships vs windspeed
# sparse spectra-  how does this relate to biological events
# Then run SVD to find representative specta- how does this relate to percentile spectra

rm(list=ls())

library(data.table)
library(rpca)
library(rsvd)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(dplyr)
library(plyr)

#-----------------------------------------------
#INPUT FILES
#-----------------------------------------------

#SET UP PARAMETERS-- plotting
#-----------------------------------------------
o3b =fread(input="E:\\RESEARCH\\SanctSound\\code\\Kurt_NVSPL\\O3OB_SanctSound.csv")
clrs=colorRampPalette(rev(c("orange","purple")))(100)

#TOL files by site-deployment
#-----------------------------------------------
sdir     = "E:\\RESEARCH\\SanctSound\\data"
TOLFiles = list.files(sdir,pattern = "TOL_1h.csv",recursive = T)
PSDFiles = list.files(sdir,pattern = "PSD_1h.csv",recursive = T)
uSites   = sapply( strsplit(basename(TOLFiles),"_") , `[`, 2)

#Time periods of interest for Frontiers paper
#-----------------------------------------------
timesInterst = read.csv( "E:\\RESEARCH\\SanctSound\\data\\AnalysisPeriods.csv")

#Combined files with OB levels, windspeed, AIS
#-----------------------------------------------
dirSPL  = paste0("E:\\RESEARCH\\SanctSound\\data2\\combineFiles")
nFiles  = length( list.files(path=dirSPL, pattern = "FrontiersDataHR_2", full.names=TRUE, recursive = FALSE))
cFiles  =         list.files(path=dirSPL, pattern  = "FrontiersDataHR_2", full.names=TRUE, recursive = FALSE)
#need to automate which file to choose based on usite!!
cdf = read.csv ( cFiles[10] )
head(cdf)

#-----------------------------------------------
#Combined files with detections, tide, windspeed (no AIS)
#-----------------------------------------------
dirSPL    = paste0("E:\\RESEARCH\\SanctSound\\data2\\combineFiles_bySite")
nFiles    = length( list.files(path=dirSPL, pattern = "CombinedDataDets_Hourly", full.names=TRUE, recursive = FALSE))
detFiles  = list.files(path=dirSPL,pattern  = "CombinedDataDets_Hourly", full.names=TRUE, recursive = FALSE)
#need to automate which file to choose based on usite!!
# ???? DO NOT HAVE THIS FILE ????
#-----------------------------------------------
#ANALYSIS: run as separate sites
#-----------------------------------------------
# Robust principal components analysis separates a matrix into a low-rank plus sparse component
#a method for the robust seperation of a a rectangular (m, n) matrix A into a low-rank component L and a sparse comonent S:

setwd(sdir)
#TOL variation-- hourly one third octave bands
#-----------------------------------------------
dfCI01=as.data.frame ( fread(input=TOLFiles[1])) 
(names(dfCI01))
df   = as.data.frame ( fread(input=TOLFiles[47]))
stmp = sapply( strsplit(basename(TOLFiles[47]),"_") , `[`, 2)
#format the date column 
df$`yyyy-mm-ddTHH:MM:SSZ`
df$DateF     = as.POSIXct( gsub(".000Z", "", gsub("T", " ", df$`yyyy-mm-ddTHH:MM:SSZ`)), tz = "GMT" )
df$DateFday  = as.Date(df$DateF)
#truncate to time period of interest
ttmp = timesInterst[timesInterst$Site == stmp,]
nv   = df[ (df$DateFday >= as.Date(ttmp$DateStart,format = "%m/%d/%Y") & df$DateFday  <= as.Date(ttmp$DateEnd,format = "%m/%d/%Y") ),]
#format the data  
hix  = as.numeric( gsub("TOL_","",names(nv)[5:34]) )
Nv   = as.matrix( nv[,2:31] ) #dB
NvP  = as.matrix(10^(Nv/10)) #pressure units

#rpca
nvpcaTOL = rrpca(Nv,trace=T)
zl =  c(min(min(nvpcaTOL$L),min(nvpcaTOL$S)),max(Nv))

#are these all the low-rank spectra?
NvM <- reshape :: melt(t(nvpcaTOL$L)  )
p2=ggplot(NvM, aes(X1, value, group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("Low rank")
#are these all the sparse spectra
NvMs <- reshape :: melt(t(nvpcaTOL$S)  )
p3=ggplot(NvMs, aes(X1, value,group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("sparse")
#low rank + sparse gets back to origional data
NvMs <- reshape :: melt(t(nvpcaTOL$S) +t(nvpcaTOL$L)  )
p4=ggplot(NvMs, aes(X1, value,group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("L+S")
#spectra from origional data
NvM1 <- reshape :: melt(t(Nv)  )
NvM1$X1 = as.numeric( gsub("TOL_","",NvM1$X1) )
p1=ggplot(NvM1, aes(log10(X1), value, group = as.factor(X2)) ) + geom_line(alpha=.1) + ggtitle("Original Data")
grid.arrange(p1,p2,p3,p4)
grid.arrange(p2,p3,nrow=1)

#add timeStamp, windspeed, SPL in same band, to LR-SVDLF
FQ_LRLF = hix[11:18] # frequency range of interest- 125- 630
indf = nvpcaTOL$L   #low rank data

SPLLRLF = 10*log10( rowMeans(10^(indf[,11:18]/10)) ) #low rank data
SPLLR   = 10*log10( rowMeans(10^(Nv[,11:18]/10)) )   #original data

dfLL = as.data.frame( cbind(as.character(df$DateF), SPLLR, SPLLRLF) )
colnames(dfLL) = c("DateTime", "SPL", "LRSPL")
dfLL$SPL = as.numeric(as.character(dfLL$SPL ))
dfLL$LRSPL = as.numeric(as.character(dfLL$LRSPL ))
dfLL$Diff = dfLL$LRSPL -dfLL$SPL 
dfLL$DateTimeF = as.POSIXct(dfLL$DateTime, tz="GMT")
dfLL$pos = dfLL$Diff < 0

par(mfrow = c(3, 1))
plot(dfLL$DateTimeF, dfLL$SPL)
plot(dfLL$DateTimeF, dfLL$LRSPL)
plot(dfLL$DateTimeF, dfLL$Diff)

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

