# updated 043024
# Anne-Cecile Lesage
# UTMB
# Look at fragility maps


library(readxl)
library(stringr)
library(pracma)
library(png)
library(ggplot2)
library(viridis)


thresholdmean_fragility <- function(fragrank, stimes, t_step, threshold_start, threshold_end, threshold = 0.5) {

  n_windows <- dim(fragrank)[2]

  # convert from input t_start and t_end to timewindow indices
  tw_start <- floor(which.min(abs(threshold_start-stimes)/t_step))
  tw_end <- floor(which.min(abs(threshold_end-stimes)/t_step))
  if (tw_end > n_windows) { tw_end <- n_windows }

  # subset fragility matrix to specified timewindows
  mat <- fragrank[,tw_start:tw_end]

  avg_f <- rowMeans(mat)
  elecfm <- which(avg_f > threshold)

  return(elecfm)
}

epoch_time_window = c(-10,30)
t_window = 250
t_step = 125
sz_onset = 0

resthresholdtot<-data.frame()


pts <- dipsaus::parse_svec("1")

path='/Volumes/bigbrain/Fragility2024/FragilityResultsACLNoNorm/'
patient_xls <- readxl::read_xlsx("/Volumes/bigbrain/Fragility2024/FragilityEEGDataset_pipeline_update_042624.xlsx")
#pipeline_xls <- readxl::read_xlsx("/Volumes/OFZ1_T7/karaslab/rave_data/bids_dir/FragilityEEGDataset/FragilityEEGDataset_pipeline.xlsx")
patient_xls$subject[pts]


threshold_start = -10
threshold_end = 20
threshold=0.6

ploton=TRUE

irvec=c()

for(i in pts){

#i=1

  subject_code <- patient_xls$subject[i]
  project <- patient_xls$project[i]
  electrodes <- dipsaus::parse_svec(patient_xls$good_electrodes[i])
  soz<-dipsaus::parse_svec(patient_xls$soz[i])
  fs<-patient_xls$sample_rate[i]
  if(class(fs)=='character'){
    fs=as.numeric(fs)
  }
  ictal_runs <- dipsaus::parse_svec(patient_xls$ictal_runs[i])
  n_tps=(epoch_time_window[2]-epoch_time_window[1])*fs
  channelfile=paste(path,subject_code,'/',subject_code,'_ses-presurgery_task-ictal_acq-ecog_channels.xls',sep="")
  channelDefs<-read_xls(channelfile)
  channelDefs<-channelDefs[electrodes,]

  n_elec=length(electrodes)
  insoz=electrodes%in%soz
  displayChannelsoz=which(insoz==TRUE)
  displayNamesoz <- channelDefs$name[insoz]
  elecsoz=displayChannelsoz
  elecsozc=which(insoz==FALSE)
  elecsozsozc=c(elecsoz,elecsozc)
  elecnum <- channelDefs$name[elecsozsozc]
  colorelec<-elecnum
  nsoz=length(elecsoz)
  colorelec[1:n_elec]="black"
  colorelec[1:nsoz]="blue"

  # Number of steps
  n_steps <- floor((n_tps - t_window) / t_step) + 1
  n_elec<-length(electrodes)
  stimes=(seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]
  irun=1

  quantilesave=c()
  for(j in ictal_runs){

    file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_fragilitynorank.csv',sep="")
    fragilitymap=read.csv(file)

    tag=paste(subject_code,'_seizure',as.character(j))
    print(tag)

    ncolfrag=ncol(fragilitymap)-1

    elecfragfile<-fragilitymap[,1]
    fragcol=seq_len(ncolfrag)+1
    # if(fragilitymap>0.1){
    #
    # }

    fragilitymap<-fragilitymap[,fragcol]

    fragilitymap<-fragilitymap[elecsozsozc,]
    stimes=stimes[1:ncolfrag]
    frag_data <- expand.grid(Time = stimes, Electrode = elecnum)
    frag_data$Value <- c(t(fragilitymap))
    titlepng=paste(subject_code,'Seizure',as.character(j),sep=" ")

    if(ploton==TRUE){
      #plot organized fragility map with soz electrodes on bottom
      ggplot(frag_data, aes(x = Time, y = Electrode, fill = Value)) +
        geom_tile() +
        ggtitle(titlepng)+
        scale_fill_viridis(option = "turbo") +
        labs(x = "Time (s)", y = "Electrode") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 5,colour=colorelec),     # Adjust depending on electrodes
        )

      file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_fragilitylnonorm.png',sep="")
      ggsave(file)
    }


    for(ii in 1:n_elec){
      for(jj in 1:n_steps){
        if(fragilitymap[ii,jj]>0.01) fragilitymap[ii,jj]=0.01
      }
    }


    maxfrag=max(fragilitymap)
    fragilitymap=(maxfrag-fragilitymap)/maxfrag

    stimes=stimes[1:ncolfrag]
    frag_data <- expand.grid(Time = stimes, Electrode = elecnum)
    frag_data$Value <- c(t(fragilitymap))
    titlepng=paste(subject_code,'Seizure',as.character(j),sep=" ")

    if(ploton==TRUE){
      #plot organized fragility map with soz electrodes on bottom
      ggplot(frag_data, aes(x = Time, y = Electrode, fill = Value)) +
        geom_tile() +
        ggtitle(titlepng)+
        scale_fill_viridis(option = "turbo") +
        labs(x = "Time (s)", y = "Electrode") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 5,colour=colorelec),     # Adjust depending on electrodes
        )

      file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_fragilitylow.png',sep="")
      ggsave(file)
    }

#
#     dum=rank(fragilitymap)
#     fragrank=matrix(dum, nrow(fragilitymap), ncol(fragilitymap))
#
#     maxfrank=max(fragrank)
#     fragrank=(maxfrank-fragrank)/maxfrank
#     frag_data$Value <- c(t(fragrank))
#
#     if(ploton==TRUE){
#       #plot organized fragility map with soz electrodes on bottom
#       ggplot(frag_data, aes(x = Time, y = Electrode, fill = Value)) +
#         geom_tile() +
#         ggtitle(titlepng)+
#         scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +
#         #scale_fill_viridis(option = "turbo") +
#         labs(x = "Time (s)", y = "Electrode") +  #
#         theme_minimal() +
#         theme(
#           axis.text.y = element_text(size = 5,colour=colorelec),     # Adjust depending on electrodes
#         )
#
#       file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_fragilityrank.png',sep="")
#       ggsave(file)
#     }
  }



}


