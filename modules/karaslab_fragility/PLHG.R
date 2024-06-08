# updated 022924
# Anne-Cecile Lesage
# UTMB

# % Implementation of PLHG as described in: 
# % Ictal high frequency oscillations distinguish two types of seizure territories in humans
# % Shennan A. Weiss, Garrett P. Banks, Guy M. McKhann, Jr, Robert R. Goodman, Ronald G. Emerson, Andrew J. Trevelyan, and Catherine A. Schevon 
# % Brain. 2013 Dec; 136(12): 3796?3808.
# 
# %%%%%% Inputs: 
#   
# % ts: time series data, with format: samples x channels. 
# 
# % tsBaseline: a baseline clip retrieved from the 30 seconds preceding
# % seizure onset. Again, format is samples x channels. 
# 
# %%%%%%% Outputs: 
#   
# % plhgMaster: PLHG values in format time x channels. Time values correspond
# % to the values provided in timeVals. 
# 
# % timeVals: a vector of length size(PLHG,1), corresponding to the time of
# % each row of plhgMaster. 
# 
# % sigLeads: 
# 

# library to read Matlab data formats into R
library(R.matlab)
library(pracma)

library(writexl)
library(reticulate)
library(ggplot2)
library(gsignal)

# importing python libraries
#use_python("/Library/Frameworks/Python.framework/Versions/3.11/bin/python3") to be edited for specific computer python 3 path
nump <- import("numpy")


## Read patient Rave preprocessed data
patient='pt01' # patient name
ir=1 # ictal run index
ravepreproc='/Users/aclesage/rave_data/data_dir/retrospective/'
raveraw='/Users/aclesage/rave_data/raw_dir/'

if(patient=='pt01'){
  
  goodelec=c(1:4,7:36,42:43,46:69,72:95)
  pathvolt=paste(ravepreproc,'pt01/rave/preprocess/voltage/',sep="")
  carfile=paste(ravepreproc,'pt01/rave/data/reference/ref_1-4,7-36,42-43,46-69,72-95.h5',sep="")
  channelfile=paste(raveraw,'pt01/run1/sub-pt01_ses-presurgery_task-ictal_acq-ecog_run-01_channels.xls',sep="")
  fs=1000
  tonsetmultirun=c(75.95,93,108.81,127.72)
  soz=c(33,34,62:69)
  resect=c(33,34,62:69)
  displayelec<- c(33,34,62:69, 88:91)
  #dispbrainelec<-c(33:36,48:69,84:91)
  scaling=1000000 
  
}

#insoz=goodelec%in%soz

library(readxl)
channelDefs <- read_xls(channelfile)
channelDefs <- channelDefs[goodelec,]

# set directory
setwd(pathvolt)
# show 
getwd()




tonset=tonsetmultirun[ir]
signalfile='electrode_1.h5'
groupfile=paste('/notch/run',as.character(ir),sep="")

ie=goodelec[1]
signalfile=paste('electrode_',as.character(ie),'.h5',sep="")
h5ls(signalfile)
en <- h5read(signalfile, groupfile)

h5ls(carfile)
groupcar=paste('/voltage/run',as.character(ir),sep="")
car <- h5read(carfile, groupcar)


# Number of time points
nt=length(en) 
## Number of electrodes
nel=length(goodelec) 

# Complete time serie of good electrode
# each column is an electrode
timeSeries=matrix(0,nt,nel)


for(ii in 1:nel){
  ie=goodelec[ii]
  signalfile=paste('electrode_',as.character(ie),'.h5',sep="")
  en <- h5read(signalfile, groupfile)
  timeSeries[1:nt,ii]=en[1:nt]-car[1:nt]
  #print(signalfile)
}


###########################
#Epoch signal -30:30s around seizure
##############################
# Seizure onset time point
itonset=tonset*fs
epochdur=20
epochm=5
epochp=25
# epoch window starting time point
its=itonset-epochm*fs
# epoch window ending time point
ite=itonset+epochp*fs

# Number of epoched signal time points
ntep=ite-its+1 
# starting time in seconds
ts=(its)/fs
tsepoch=matrix(0,ntep,nel)
tsepochv=matrix(0,ntep,nel)

scalingdata=max(timeSeries)

for(ii in 1:nel){
  tsepoch[1:ntep,1:nel]=timeSeries[its:ite,1:nel]
  #tsepochv[1:ntep,1:nel]=timeSeries[its:ite,1:nel]/scalingdata
}

scalingdata=max(tsepoch)/5

for(ii in 1:nel){
  #tsepoch[1:ntep,1:nel]=timeSeries[its:ite,1:nel]
  tsepochv[1:ntep,1:nel]=timeSeries[its:ite,1:nel]/scalingdata
}

tsEpoch <- data.frame(tsepoch)
colnames(tsEpoch) <- channelDefs$name

tsEpochv <- data.frame(tsepochv)
colnames(tsEpochv) <- channelDefs$name

#####################
# Signal visualization
#####################
## select channels and standardize signal

gaps <- 2

displayelec<-soz


indisplay=goodelec%in%displayelec
displayChannels=which(indisplay==TRUE)
displayNames <- channelDefs$name[displayChannels]
displayNames

indisplays=goodelec%in%soz
displayChannelsoz=which(indisplays==TRUE)
displayNamesoz <- channelDefs$name[displayChannelsoz]
displayNamesoz

indisplayr=goodelec%in%resect
displayChannelresect=which(indisplayr==TRUE)
displayNameresect <- channelDefs$name[displayChannelresect]
displayNameresect

plotData <- tsEpochv[,displayChannels]
for(i in seq_along(plotData)){
  plotData[, i] <- (plotData[, i]- mean(plotData[, i]))+
    (ncol(plotData)-i)*gaps
}

plot(plotData[, 1],type="l" ,cex=0.1,
     ylim = range(plotData), yaxt = "n")
for(i in 2:ncol(plotData)){
  lines(plotData[, i])
}
axis(2, at = rev(seq_along(displayChannels) - 1)*gaps,
     labels = displayNames,las=1)

# ending time in seconds
te=(ite)/fs

sizeWindow <- 3000 
sizeSkip <- 333


#Define tsBaseline
ntb=30*fs
tsBaseline=matrix(0,ntb,nel)

its=5*fs
ite=its+ntb-1

for(ii in 1:nel){
  tsBaseline[1:ntb,1:nel]=timeSeries[its:ite,1:nel]
}


# number of windows for PLHG analysis  
nw<-floor(ntep-sizeWindow)/sizeSkip

stepsBuffer=zeros(ntb,nw)

for(ii in 1:nw){
  stepsBuffer[1:ntb,ii]<-(ii-1)*sizeSkip+1:sizeWindow
}
for(ii in 1:sizeWindow){
  if(stepsBuffer[ii,nw]>ntep) stepsBuffer[ii,nw]<-0
}

#check<-stepsBuffer[,nw]
timeVals<-stepsBuffer[1,]

PLVMaster<-matrix(0,nw,nel)

nyquist <- fs/2
transitionWidth <- 0.1 
dt<-1/fs
##########################################
# Filter signal low frequency

filterwindow<-c(4,30)
fvec   <- vector(mode="numeric", length=6)
fvec[1] <- 0.0
fvec[2] <- (1 - transitionWidth) * filterwindow[1]/nyquist
fvec[3] <- filterwindow[1]/nyquist
fvec[4] <- filterwindow[2]/nyquist
fvec[5] <- (1 + transitionWidth) * filterwindow[2]/nyquist
fvec[6] <- 1.0
idealresponse<-c(0, 0, 1, 1, 0, 0)

sprintf(" Filter Data in low Frequency Band 4-30 Hz")
start <- Sys.time()
print(start)
# build firls filter
fir_4_30<-firls(499,fvec,idealresponse)
filter_4_30 <- filtfilt(fir_4_30,tsepoch)
hilbert_4_30<-hilbert(filter_4_30)
phi_4_30<-nump$angle(hilbert_4_30)
end <- Sys.time()
print(end - start)

##########################################
# Filter signal high gamma

filterwindow<-c(80,150)
fvec   <- vector(mode="numeric", length=6)
fvec[1] <- 0.0
fvec[2] <- (1 - transitionWidth) * filterwindow[1]/nyquist
fvec[3] <- filterwindow[1]/nyquist
fvec[4] <- filterwindow[2]/nyquist
fvec[5] <- (1 + transitionWidth) * filterwindow[2]/nyquist
fvec[6] <- 1.0
idealresponse<-c(0, 0, 1, 1, 0, 0)

sprintf(" Filter Data in Frequency Band 80-150 Hz")
start <- Sys.time()
print(start)

# build firls filter
fir_80_150<-firls(499,fvec,idealresponse)
filter_80_150 <- filtfilt(fir_80_150,tsepoch)
hilbert_80_150<-hilbert(filter_80_150)

a_80_150<-abs(hilbert_80_150)
hilbert_a_80_150<-hilbert(a_80_150)
phi_a_80_150<-nump$angle(hilbert_a_80_150)

end <- Sys.time()
print(end - start)


##########################################
# Filter pre seizure baseline

sprintf(" Filter pre seizure baseline in Frequency Band 80-150 Hz")
start <- Sys.time()
print(start)
filter_80_150_baseline <- filtfilt(fir_80_150,tsBaseline)
hilbert_80_150_baseline<-hilbert(filter_80_150_baseline)
a_80_150_baseline<-abs(hilbert_80_150_baseline)
end <- Sys.time()
print(end - start)

a_80_150_baseline<-colMeans(a_80_150_baseline)


plhgMaster<-matrix(0,nw,nel)


for(jj in 1:nw){
  #jj<-1
  #print(jj)
  currentTime<-stepsBuffer[,jj]
  phi_4_30_jj<-phi_4_30[currentTime,]
  phi_a_80_150_jj<-phi_a_80_150[currentTime,]
  
  PLV<-abs(colMeans(exp(1i*(phi_4_30_jj-phi_a_80_150_jj))))
  
  a_80_150_jj<-a_80_150[currentTime,]
  a_80_150_norm_jj<-colMeans(a_80_150_jj)/a_80_150_baseline
  PLHG=a_80_150_norm_jj*PLV
  
  PLVMaster[jj,]<-PLV
  plhgMaster[jj,]<-PLHG
}
# 
# mplhg<-colMeans(plhgMaster)
# sdplhg<-apply(plhgMaster,2,sd)
# 
# # compute Zscore
# for(jj in 1:nel){
#   mc<-mplhg[jj]
#   sdc<-sdplhg[jj]
#   plhgMaster[,jj]=(plhgMaster[,jj]-mc)/sdc
# }
# 
# maxVal<-apply(plhgMaster,2,max)
# mu<-mean(maxVal)
# sigma<-sd(maxVal)
# coeffVar<-sigma/mu

mu=mean(plhgMaster)
sigma<-sd(plhgMaster)
sigThres<-mu+sigma*2.5
sigLeads<-maxVal>sigThres

sigTime<-matrix(NaN,1,nel)

votethres   <- vector(mode="numeric", length=nel)

jj<-3
for(jj in 1:nel){
  if(sigLeads[jj]==TRUE){
    votethres[jj]=1
  }
  currentInd<-which(plhgMaster[,jj]>=sigThres)
  if(length(currentInd)>0){
    sigTime[1,jj]=timeVals[currentInd[1]]/1000
  }
}

ResPLHG<-data.frame(channelDefs$name,votethres,sigTime[1,])   

write_xlsx(ResPLHG,"/Users/aclesage/Documents/grant_writing/R21March2024/FigPreliminaryWork/ResPLHG.xlsx")