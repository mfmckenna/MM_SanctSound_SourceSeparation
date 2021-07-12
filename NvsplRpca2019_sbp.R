rm(list=ls())
library(data.table)
library(rpca)

o3b <- fread(input="I:/My Drive/share_files/Chronic-Transient_Discrimination/code/NVSPL_Kurt/O3OB.csv")
sdir <- "./" # using Tools..Global Options..Default Working Directory to be in the right place
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

nv <- fread(input="I:/My Drive/share_files/Chronic-Transient_Discrimination/code/NVSPL_Kurt/NVSPL_OLYM001_2010_06_11_07.csv")
hix <- grep("^H[^u]",names(nv))
clrs <- colorRampPalette(rev(c("orange","purple")))(100)
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
