#finding low rank and sparce features in hourly, third octave band spectra in SanctSound project

#use SB01 and CI01 to run rpca (robust PCA):
# low-rank spectra- how does this relate to ships vs windspeed
# sparce spectra-  how does this relate to biological events
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
cdf = read.csv ( cFiles[1] )
head(cdf)

#Combined files with detections, tide, windspeed (no AIS)
#-----------------------------------------------
dirSPL    = paste0("E:\\RESEARCH\\SanctSound\\data2\\combineFiles_bySite")
nFiles    = length( list.files(path=dirSPL, pattern = "CombinedDataDets_Hourly", full.names=TRUE, recursive = FALSE))
detFiles  = list.files(path=dirSPL,pattern  = "CombinedDataDets_Hourly", full.names=TRUE, recursive = FALSE)
#need to automate which file to choose based on usite!!
cdetf = read.csv (detFiles[1] )
head(cdetf)
cdetf = cdetf[cdetf$Site == uSites[1],]
#Truncate to Frontiers period of interest
tmpAP   = timesInterst[timesInterst$Site == uSites[1],]
cdetft  = cdetf[ (as.Date(cdetf$Date) >= as.Date(tmpAP$DateStart,format = "%m/%d/%Y") & as.Date(cdetf$Date)  <= as.Date(tmpAP$DateEnd,format = "%m/%d/%Y") ),]
cdetft  = cdetft[!is.na(cdetft$Deployment),] #had to remove NAs, effort control

rm(tmpAP,cdetf)
#-----------------------------------------------
#ANALYSIS: run as separate sites
#-----------------------------------------------
# Robust principal components analysis separates a matrix into a low-rank plus sparse component
#a method for the robust seperation of a a rectangular (m, n) matrix A into a low-rank component L and a sparse comonent S:

setwd(sdir)
#TOL variation-- hourly one third octave bands
#-----------------------------------------------
df   = as.data.frame ( fread(input=TOLFiles[2]))
stmp = sapply( strsplit(basename(TOLFiles[2]),"_") , `[`, 2)
#format the date column 
df$yyyy_mm_ddTHH_MM_SSZ
df$DateF     = as.POSIXct( gsub(".000Z", "", gsub("T", " ", df$yyyy_mm_ddTHH_MM_SSZ)), tz = "GMT" )
df$DateFday  = as.Date(df$DateF)
#truncate to time period of interest
ttmp = timesInterst[timesInterst$Site == stmp,]
nv   = df[ (df$DateFday >= as.Date(ttmp$DateStart,format = "%m/%d/%Y") & df$DateFday  <= as.Date(ttmp$DateEnd,format = "%m/%d/%Y") ),]
#format the data  
hix  = as.numeric( gsub("TOL_","",names(nv)[2:31]) )
Nv   = as.matrix( nv[,2:31] ) #dB
NvP  = as.matrix(10^(Nv/10)) #pressure units

#rpca
nvpcaTOL = rrpca(Nv,trace=T)
zl =  c(min(min(nvpcaTOL$L),min(nvpcaTOL$S)),max(Nv))
par(mfrow=c(3,1),oma=c(4,2,0,0),mar=c(0,4.5,0,0)+.1,mgp=c(3.5,0.75,0),las=1)
image(1:dim(Nv)[1],log(hix),Nv)
#image(1:dim(NvP)[1],log(hix),NvP)
#image(1:dim(nvpcaTOL$L)[1],log(hix),as.matrix(nvpcaTOL$L))  #low rank
#image(1:dim(nvpcaTOL$S)[1],log(hix),as.matrix(nvpcaTOL$S)) #sparce matrix


#are these all the low-rank spectra?
NvM <- reshape :: melt(t(nvpcaTOL$L)  )
p2=ggplot(NvM, aes(X1, value, group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("Low rank")
#are these all the sparce spectra
NvMs <- reshape :: melt(t(nvpcaTOL$S)  )
p3=ggplot(NvMs, aes(X1, value,group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("Sparce")
#low rank + sparce gets back to origional data
NvMs <- reshape :: melt(t(nvpcaTOL$S) +t(nvpcaTOL$L)  )
p4=ggplot(NvMs, aes(X1, value,group = as.factor(X2))) + geom_line(alpha=.1) + ggtitle("L+S")
#spectra from origional data
NvM1 <- reshape :: melt(t(Nv)  )
NvM1$X1 = as.numeric( gsub("TOL_","",NvM1$X1) )
p1=ggplot(NvM1, aes(log10(X1), value, group = as.factor(X2)) ) + geom_line(alpha=.1) + ggtitle("Original Data")
grid.arrange(p1,p2,p3,p4)
grid.arrange(p2,p3,nrow=1)

#what do these low-rank spectra represent- looks like different wind speeds
#LOW RANK is really same spectra with sparce variation removed.... how can I collapse these? 
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

#SPARCE what do the sparce represent? 
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

#PLOT: Sparce matrix with MIDSHIPMAN presence
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
  ggtitle(paste0("Sparce rank- midshipman only (N=",samps,")")) + 
  scale_color_discrete(name = "Detections") +
  #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30])) +
  theme(legend.position = "none")
#PLOT: Sparce matric with BLUE present
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
  ggtitle(paste0("Sparce rank- blue whale (N=",samps,")")) + 
  scale_color_discrete(name = "Detections") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30]))+
  theme(legend.position = "none")

#PLOT: Sparce matric with VESSEL present, only
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
  ggtitle(paste0("Sparce rank- vessel (N=",samps,")")) + 
  scale_color_discrete(name = "Detections") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(name = "Frequency 1/3 octave bands", breaks = c(1,10,20,30), labels = c(hix[1], hix[10], hix[20],hix[30]))+
  theme(legend.position = "none")

#PLOT: Sparce matric with VESSEL present
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
  ggtitle("Sparce rank- explosions present") + 
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
heatmap(indf) # cluters of data in frequency (wind, low frequency crap, snapping shrimp, and clusters of data in time)
hix #1/3 otave bands for reference to graphic

hh = dist(indf) %>% hclust
dataMatrixOrdered <- indf[hh$order, ]
par(mfrow = c(1, 3))
image(t(dataMatrixOrdered)[, nrow(dataMatrixOrdered):1])
plot(rowMeans(dataMatrixOrdered), 721:1, , xlab = "Row Mean", ylab = "Row", pch = 19)
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
plot(Lsvd$u[, 1], 721:1, , ylab = "Row", xlab = "First left singular vector", pch = 19)
plot(Lsvd$v[, 1], xlab = "Column", ylab = "First right singular vector", pch = 19)
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
# SPARCE MATRIX: hierarchical clustering algorithm, sparce matrix 
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
image(1:30, 1:721, t(indf)[, nrow(indf):1]) #lsparce matrix
image(1:30, 1:721, t(B)[, nrow(B):1])       #scaled sparce matrix

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

#sparce
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

