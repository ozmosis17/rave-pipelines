library(glmnet)

source("./modules/karaslab_fragility/jiefei-fragility.R")

library(rhdf5)

## Determine the signal strength (larger scaling -> weaker signal)
signalScaling <- 1000000

library(readxl)


## Read signal data, each column will be one electrode
library(R.matlab)

## Read patient Rave preprocessed data
patient='pt01' # patient name

if(patient=='pt01'){

  goodelec=c(1:4,7:24,26:36,42:43,46:54,56:70,72:95)
  runs=3
  pathvolt='/Volumes/OFZ1_T7/karaslab/rave_data/data_dir/FragilityEEGDataset/PT01/rave/preprocess/voltage/'
  carfile='/Volumes/OFZ1_T7/karaslab/rave_data/data_dir/FragilityEEGDataset/pt01/rave/data/reference/ref_1-4,7-24,26-36,42-43,46-54,56-70,72-95.h5'
  channelfile='/Volumes/OFZ1_T7/karaslab/rave_data/raw_dir/PT01/run01/sub-pt01_ses-presurgery_task-ictal_acq-ecog_run-01_channels.xls'
  fs=1000
  tonsetmultirun=c(75.95,93,108.81,127.72)
  soz=c(33,34,62:69)
  displayelec<- c(33,34,62:69, 88:91)


}else if(patient=='jh103'){

  runs=3
  pathvolt='/Volumes/OFZ1_T7/karaslab/rave_data/data_dir/retrospective/jh103/rave/preprocess/voltage/'
  carfile='/Volumes/OFZ1_T7/karaslab/rave_data/data_dir/retrospective/jh103/rave/data/reference/ref_1-4,7-12,15-23,25-33,47-63,65-66,69-71,73-110.h5'
  channelfile='/Volumes/OFZ1_T7/karaslab/Documents/RAVEProjects/OpenNeuroFragility/jh103channels.xls'
  fs=1000
  goodelec=c(1:4,7:12,15:23,25:33,47:63,65,66,69:71,73:110)
  soz<-c(86,94,15:21,24,25:33)
  tonsetmultirun=c(59.89,59.89)
  displayelec<-c(1,2,15:21,25:33,75:78,86,94)

}

## Read signal data from RAVE preprocessed raw signal
# Epoch signal
# each column will be one electrode


insoz=goodelec%in%soz

x <- read.delim('/Volumes/OFZ1_T7/karaslab/rave_data/raw_dir/PT01/run01/sub-pt01_ses-presurgery_task-ictal_acq-ecog_run-01_channels.tsv')
x <- x[goodelec,]

library(readxl)
channelDefs <- read_xls(channelfile)
channelDefs <- channelDefs[goodelec,]

# set directory
setwd(pathvolt)
# show
getwd()



ir=1 # ictal run index
tonset=tonsetmultirun[1] # choose trial's epoch start time
signalfile='electrode_1.h5'
groupfile='/notch/presurgery_ictal_ecog_01'

ie=goodelec[1]
signalfile=paste('electrode_',as.character(ie),'.h5',sep="")
h5ls(signalfile)
en <- h5read(signalfile, groupfile)

h5ls(carfile)
groupcar='/voltage/presurgery_ictal_ecog_01'
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

# outputfile=paste('TimeSeriesRun01',patient,'.mat',sep='')
# writeMat(outputfile,timeSeries=timeSeries)

###########################
#Epoch signal -10:10s around seizure
##############################
# Seizure onset time point
itonset=tonset*fs
# epoch window starting time point
its=itonset-10*fs
# epoch window ending time point
ite=itonset+10*fs

# Number of epoched signal time points
ntep=ite-its+1
# starting time in seconds
ts=(its)/fs
# ending time in seconds
te=(ite)/fs


tsepoch=matrix(0,ntep,nel)

for(ii in 1:nel){
  tsepoch[1:ntep,1:nel]=timeSeries[its:ite,1:nel]
}

# tsepoch is a mxn matrix where m = # timepoints in epoch, n = # electrodes

tsEpoch <- data.frame(tsepoch/signalScaling)
colnames(tsEpoch) <- channelDefs$name

#####################
# Signal visualization
#####################
## select channels and standardize signal

gaps <- 2


indisplay=goodelec%in%displayelec
displayChannels=which(indisplay==TRUE)
displayNames <- channelDefs$name[displayChannels]
displayNames


plotData <- tsEpoch[,displayChannels]
for(i in seq_along(plotData)){
  #print(i)
  # plotData[, i] <- (plotData[, i]- mean(plotData[, i]))*4 +
  #   (ncol(plotData)-i)*gaps
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


#####################
# Compute fragility at a time window
#####################
## The number of samples in each window
ntw <- 250
## sliding time increment
ntsl <- 125
## The index of the window we are going to calculate fragility
iw <- 1
## Sample indices for the selected window
si <- seq_len(ntw) + (iw-1)*ntsl


## measurements at time point t
xt <- tsEpoch[si,]
## measurements at time point t plus 1
xtp1 <- tsEpoch[si + 1,]


## Coefficient matrix A (adjacency matrix)
## each column is coefficients from a linear regression
## formula: xtp1 = xt*A
A <- ridge(xt,xtp1, intercept = F, lambda = 0.001)
R2 <- ridgeR2(xt,xtp1,A)

## Check heatmap
## set diagonal to 0 to see if there is any outlier
## (diagonal is always outlier)
heatA <- A
diag(heatA) <- 0
heatmap(heatA, Colv = NA, Rowv = NA, scale="none")

## Check Goodness of Fit
hist(R2)

## checking adjacency matrix stability
AEigen <- eigen(A)
e <- Mod(AEigen$values)
me <- max(e)
if(me <1){
  message("Stable solution")
}else{
  message("Solution is not stable")
}


## fragility for each electrode
fragrow <- fragilityRow(A)

fragrow



#####################
# Compute Fragility matrix
#####################
## parallel computing
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores())

## The number of samples in each window
ntw <- 250
## sliding time increment
ntsl <- 125
## total number of possible windows
nWins <- floor(nrow(tsEpoch)/ntsl-1)
## iw: The index of the window we are going to calculate fragility
res <- foreach(iw = seq_len(nWins),
               .combine = cbind,
               .packages = "glmnet")%dopar%{
                 ## Sample indices for the selected window
                 si <- seq_len(ntw) + (iw-1)*ntsl
                 ## measurements at time point t
                 xt <- tsEpoch[si,]
                 ## measurements at time point t plus 1
                 xtp1 <- tsEpoch[si + 1,]

                 ## Coefficient matrix A (adjacency matrix)
                 ## each column is coefficients from a linear regression
                 ## formula: xtp1 = xt*A + E
                 A <- ridge(xt,xtp1, intercept = F)

                 ## fragility vector
                 fragilityRow(A)
               }
resnaked=res
rownames(res) <- channelDefs$name
colnames(res) <- as.character(round((seq_len(nWins)-1)*ntsl/1000-10))

dfres<-data.frame(res)
library("writexl")
write_xlsx(dfres,"pt01run01lambda1e-4.xlsx")

## use ranking to increase contrast
res2 <- matrix(rank(res), nrow(res), ncol(res))
attributes(res2) <- attributes(res)

# visualize entire fragility map (all electrodes)
library(RColorBrewer)
coul <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(256)
heatmap(res2, Colv = NA, Rowv = NA,
        scale="none", col = coul, xlab = "time")

# visualize on display electrodes
displayChannelsrev=rev(displayChannels)
resdisplay=res[displayChannelsrev,]
res2d <- matrix(rank(resdisplay), nrow(resdisplay), ncol(resdisplay))
attributes(res2d) <- attributes(resdisplay)
# visualize entire fragility map (all electrodes)
library(RColorBrewer)
coul <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(256)
heatmap(res2d, Colv = NA, Rowv = NA,
        scale="none", col = coul, xlab = "time")

# # visualize electrodes organized as in soz and out soz
# sozelec=which(insoz==TRUE)
# sozelecc=which(insoz==FALSE)
# ressoz=resnaked[sozelec,]
# ressozc=resnaked[sozelecc,]
# namesoz=channelDefs$name[sozelec]
# namesozc=channelDefs$name[sozelecc]
