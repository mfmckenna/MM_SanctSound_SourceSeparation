#FINDING low rank and sparse features in hourly, third-octave band spectra in SanctSound project

# purpose: separate transients and better estimates of ambient
rm(list=ls())

#-----------------------------------------------------------------------------------------
#NOTES ####

#-----------------------------------------------------------------------------------------
#SET UP ####
#-----------------------------------------------------------------------------------------
library(data.table)
library(rpca)
library(rsvd)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(dplyr)
library(plyr)
library(svglite)

# SET directory set up (unique to MFM computer) 
#top directory with all SanctSound data products
tDir  =  "E:\\RESEARCH\\SanctSound\\data"
#directory with combined files for detentions
inDir ="E:\\RESEARCH\\SanctSound\\data2\\combineFiles3_Detections\\"
#directory with code
cDir  = "E:\\CODE\\GitHub\\MM_SanctSound_SourceSeparation\\" 

# GET sanctuary directories loop through sanctuaries
dataDir = list.dirs(path = tDir, full.names = F, recursive = F)
dirIdx  = which(nchar(dataDir)==2 )
sanct   = dataDir[dirIdx]
dataDir = list.dirs(path = tDir, full.names = T, recursive = F)
sDirs   = dataDir[dirIdx]

# LOAD input files and settings
o3b = fread(input=paste0(cDir, "O3OB_SanctSound.csv")) #1/3 octave band data details for plotting
clrs= colorRampPalette(rev(c("orange","purple")))(100) #SET UP PARAMETERS-- plotting
smth = "April" #month to plot

#-----------------------------------------------------------------------------------------
# Process data by sanctuary sites ####
#-----------------------------------------------------------------------------------------
for (ss in 1:length(sDirs)){ # START: loop through sanctuaries ss = 8
 
  indir = sDirs[ss]
  sdir  = list.dirs(indir, full.names = F, recursive = F)
  dirIdx = which(nchar(sdir) == 4 )
  sites = sdir[dirIdx]
  dataDir = list.dirs(path = indir, full.names = T, recursive = F)
  sitesDirs  = dataDir[dirIdx]
  
  for (jj in 1:length(sites)) {#START: loop through sites jj = 2
    
    sdir = sitesDirs[jj]
    
    #TOL FILES: by site-deployment
    #-----------------------------------------------------------------------------------------
    TOLFiles = list.files(sdir,pattern = "TOL_1h.csv",recursive = T)
    PSDFiles = list.files(sdir,pattern = "PSD_1h.csv",recursive = T)
    site   = unique( sapply( strsplit(basename(TOLFiles),"_") , `[`, 2) )
    cat("Processing sanctuary ", ss, " of ", length(sDirs), " (site: ", site, ") \n"  )
    setwd(sdir)
    
    
    #FORMAT DATA: TOL variation-- hourly one third octave bands, 
    #-----------------------------------------------------------------------------------------
    df = NULL
    #different headers for SB sites
    for (ii in 1:length(TOLFiles)) {
      
      dftmp = as.data.frame ( fread(input=TOLFiles[ii]))
      deployment = sapply( strsplit(basename(TOLFiles[ii]),"_") , `[`, 3)
     
      
      if (sanct[ss] == "SB") {
        names(dftmp)[1] = "yyyy-mm-ddTHH:MM:SSZ"
        names(dftmp)[6] = "TOL_31.5"  #sometimes 31.5 is different from 31_5--
        #truncate to only 
        dftmp = dftmp[,c(1, 5:34)]}
      
      if (sanct[ss] != "SB") {
        names(dftmp)[1] = "yyyy-mm-ddTHH:MM:SSZ"
        names(dftmp)[3] = "TOL_31.5" }
      
      df = rbind(df, cbind(dftmp,site,deployment))   
      
    }
    #format the date column 
    df$`yyyy-mm-ddTHH:MM:SSZ`
    df$DateF     = as.POSIXct( gsub(".000Z", "", gsub("T", " ", df$`yyyy-mm-ddTHH:MM:SSZ`)), tz = "GMT" )
    df$DateFday  = as.Date(df$DateF)
    df$Hr  = hour(df$DateF)
    df$Mth = month(df$DateFday)
    
    #remove rows with NA--- only came up with PM08
    #df[df$deployment == 13,]
    idNA = which(is.na(df))
    df =  df[complete.cases(df[ , 5:6]),]
    
    #truncate to time period of interest
    nv   = df[ (df$DateFday   >= as.Date("01/01/2019",format = "%m/%d/%Y") &
                  df$DateFday <= as.Date("12/31/2020",format = "%m/%d/%Y") ),]
    
    #prep data the data  as.data.frame(colnames(nv))
    idx  = grep("^TOL", colnames(nv))
    hix  = as.numeric( gsub("TOL_","",names(nv)[idx]) )
    Nv   = as.matrix( ( nv[,idx]) )  #dB-- why is this convert to character???
    NvP  = as.matrix(10^(Nv/20))     #pressure units
    nvDate = nv$DateF
    
   
    #ANALYSIS: run as separate sites
    #-----------------------------------------------------------------------------------------
    # Robust principal components analysis separates a matrix into a low-rank plus sparse component
    #a method for the robust separation of a a rectangular (m, n) matrix A into a low-rank component L and a sparse component S:
    input = ( NvP )# or Nv (dB units)
    lamd = max(input)^-0.5 #default 
    nvpcaTOL = rrpca(input)
    sampleHours = nrow(input)
    
    #SUMMARIZE OUTPUT:
    #-----------------------------------------------------------------------------------------
    #input data 
    Am = as.data.frame(input)   
    #low rank
    Lr = as.data.frame(nvpcaTOL$L) 
    colnames(Lr) = hix
    LrDB = 10*log10( Lr^2 )  #CHECK: min(LrDB$`63`), no negative values, just values without transients
    colnames(LrDB) = hix
    LrDB =as.matrix(LrDB)
    #sparse matrix
    Sp = as.data.frame(nvpcaTOL$S) 
    colnames(Sp) = hix
    SpDB = 10*log10( (Sp)^2 ) # negative and zero values-- does not make sense to convert back to dB
    colnames(SpDB) = hix
    
    # VISULAZE DATA: LOW RANK 
    #-----------------------------------------------------------------------------------------
    #spectra for each hour-- lots of data, so truncate to month of interest
    NvMlr   = reshape :: melt(t(LrDB)  )
    tLabel  = paste0( site, "-", as.character(min(nv$DateF)), "-", as.character( max(nv$DateF) ), " (hours = ", length(nv$DateF), ")")
    
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
    
    #original matrix
    NvMO   = reshape :: melt(t(Nv)  )
    NvMOt  = NvMlr[NvMO$X2 %in% idx,]
    pO = ggplot(NvMOt, aes(X1, value, group = as.factor(X2)))+ 
      geom_line(alpha=.05)+ scale_x_continuous(trans = "log10")+ 
      labs(title = paste0(site, ": Original (plot for ", smth, ") \n Total Hours : ",sampleHours ))+
      xlab("Frequency (1/3 octave bands)") + ylab("Hourly median SPL")+
      ylim(c(60,114))+
      theme_minimal()
    
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
      labs(title = "Low Rank") + 
      scale_color_manual(values=c('#E69F00','#56B4E9'), name = "") + #E69F00 = orange (2019), #56B4E9 = blue   (2020)
      xlab("Frequency (1/3 octave bands)") + ylab("Hourly median SPL")+
      ylim(c(60,114))+
      annotate("text", x = 50, y = 65, label = ("2019"), color = '#E69F00' )+
      annotate("text", x = 50, y = 63, label = ("2020"), color = '#56B4E9' )+
      theme_minimal()+
      theme(legend.position="none") 
    
    # sparse matrix
    NvMsp <- reshape :: melt(t(Sp)  )
    NvMspt  = NvMsp[NvMsp$X2 %in% idx,]
    pS = ggplot(NvMspt, aes(X1, (value), group = as.factor(X2)))+ 
      geom_line(alpha=.07)+ scale_x_continuous(trans = "log10")+ 
      labs(title = "Sparse", caption = tLabel )+
      xlab("Frequency (1/3 octave bands)") + ylab("Pressure units")+
      #ylim(c(-2,10)) +
      theme_minimal()
    
    grid.arrange(pO, pL, pS)
    g = arrangeGrob(pO, pL, pS)
    
    #save out image
    filename =  paste0(tDir, "\\", site, "_RRPCAplots_", sampleHours, "hrs_", smth,".png" )
    ggsave(filename, g, width = 11.5, height = 12, units = "in",dpi = 300)
    
    filename =  paste0(tDir, "\\", site, "_RRPCAplots_", sampleHours, "hrs_", smth,".svg" )
    ggsave(filename, g, width = 11.5, height = 12, units = "in",dpi = 300 )
    
    
    #can I do some sort of tile graph showing difference from ambient by day
    #original matrix (Day/hour, Frequency, dB value)
    Nv2 = as.data.frame(Nv)
    Nv2$DateF = nvDate
    names(Nv2) = c(hix,"DateF")
    Nv2m = reshape :: melt(Nv2, id.vars = "DateF",  measure.vars =  names(Nv2)[1:30])
    Nv2m$day = as.Date(Nv2$DateF)
    
    
    LrDB2 = as.data.frame(LrDB)
    LrDB2$DateF = nvDate
    LrDB2m = reshape :: melt(LrDB2, id.vars = "DateF",  measure.vars =  names(LrDB2)[1:30])
    Nv2m$Diff =Nv2m$value - LrDB2m$value
    Nv2mDay = aggregate( Nv2m$Diff, by=list(Nv2m$day,Nv2m$variable), mean, na.rm=T) 
    names(Nv2mDay) = c("Day","Frequency","Difference")
    Nv2mDay$Difference2 = Nv2mDay$Difference
    Nv2mDay$Difference2[Nv2mDay$Difference2 < 0] = 0 
    
    ggplot(Nv2mDay, aes(Day, Frequency, fill= (Difference2))) +
      geom_tile() +
      scale_fill_gradient(low="gray", high="blue") +
      labs( fill = "Transient Energy \n daily average") +
      geom_vline(xintercept=c( as.Date("2019-04-01"),as.Date("2020-04-01")) ,
                 linetype=4, colour="black") +
      annotate("text", x = as.Date("2019-04-11"), y = 29 , label = "April 2019") +
      annotate("text", x = as.Date("2020-04-11"), y = 29 , label = "April 2020") +
      xlab("") +
      ylab("") +
      theme_minimal()+
      theme(plot.margin=unit(c(0,5.75,0,2.25), "cm"),
            axis.title = element_text(family="Times",face="bold", size=10),
            axis.text  = element_text(family="Times",face="bold", size=10),
            legend.text = element_text(face="bold", size=10),
            legend.title = element_text(face="bold", size=12),
            legend.box.background = element_rect(color="white", size=1,fill="white"),
            legend.spacing.y = unit(.1, "cm"),legend.spacing.x = unit(1, "cm"),
            plot.title = element_text(face="bold", size=24, hjust = 0.5),
            plot.subtitle = element_text(face="bold", size=20, hjust = 0.5) )
    
    
    #DETECTIONS: of detection data to compare to results (output of 3_combineFiles_detections)
    #-----------------------------------------------------------------------------------------
    #not currently working...
    if (site == "SB02") {
      inFiles  = list.files(inDir, site, full.names = T)
      load(inFiles[5])#HOUR data
      load(inFiles[1])#DAY data
      
      # select time periods of interest
      dataAhrt  = dataAhr[ as.Date(dataAhr$Day)  >=  as.Date("01/01/2019",format = "%m/%d/%Y" )
                            & as.Date(dataAhr$Day) <= as.Date("12/31/2020",format = "%m/%d/%Y" ),]
      dataAhrt  = dataAhrt[!is.na(dataAhrt$Deployment),] #had to remove NAs, effort control, but does not work in the analysis
      
      dataAdyt  = iData[ (as.Date(iData$Day)  >= as.Date("01/01/2019",format = "%m/%d/%Y" ) 
                          & as.Date(iData$Day) <= as.Date("12/31/2020",format = "%m/%d/%Y" ) ),]
      dataAdyt  = dataAdyt[!is.na(dataAdyt$Deployment),] #had to remove NAs, effort control, but does not work in the analysis
      rm(tmpAP)
      
      #narw
      length( dataAdyt [ dataAdyt$northatlanticrightwhaleP == 1 & dataAdyt$Yr == 2020 & dataAdyt$Month == 4, 2] )
      length( dataAdyt [ dataAdyt$northatlanticrightwhaleP == 1 & dataAdyt$Yr == 2019 & dataAdyt$Month == 4, 2] )
      #fin whales
      length( dataAdyt [ dataAdyt$finwhaleP == 1 & dataAdyt$Yr == 2020 & dataAdyt$Month == 4, 2] )
      length( dataAdyt [ dataAdyt$finwhaleP == 1 & dataAdyt$Yr == 2019 & dataAdyt$Month == 4, 2] )
      
      #add detections to see if more "transient signal when calls present
      Nv2mDay$northatlanticrightwhaleP = 0
      idxDay = dataAdyt$Day[ dataAdyt$northatlanticrightwhaleP > 0 ]
      idxDay = idxDay[!is.na(idxDay)]
      Nv2mDay$northatlanticrightwhaleP[which(Nv2mDay$Day %in% idxDay)] = 1
      Nv2mDays = Nv2mDay[which(Nv2mDay$Frequency %in% c("50","63","80","100","125","160","200","250") ),] 
      ggplot(Nv2mDays, aes(x=Frequency, y=Difference, fill= as.factor(northatlanticrightwhaleP) )) + 
        geom_boxplot()
      
      #add prob of bio present with 
      Nv2mDay$Bio = 0
      udays = unique( dataAdyt$Day )
      for (bb in 1:length(udays)) {
        tmp = dataAdyt$PercentVessel_dailyP[dataAdyt$Day == udays[bb]]
        Nv2mDay$Bio[which(Nv2mDay$Day == udays[bb])] = tmp
      }
      
      Nv2mDays = Nv2mDay[which(Nv2mDay$Frequency %in% c("50","63","80","100","125","160","200","250") ),] 
      ggplot(Nv2mDays, aes(x=Bio, y=Difference,color = as.factor( northatlanticrightwhaleP)  )) + 
        geom_point() +
        geom_smooth(method = "glm")
      
     
      }
    
    
  }
}





LrDB$DateF = nvDate

NvMlr

NvMOt$Diff = NvMOt$value - NvMlrt$value
NvMO   = reshape :: melt(t(Nv)  )


ggplot(NvMOt, aes(as.factor(X2), as.factor(X1), fill= (Diff))) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  labs(title = paste0(site,": summary of sounds"),  fill = "Detections") +
  xlab("") +
  ylab("")
