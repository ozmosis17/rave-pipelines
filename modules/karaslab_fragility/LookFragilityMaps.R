# updated 051524
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

epoch_time_window = c(-10,20)
t_window = 250
t_step = 125
sz_onset = 0

resthresholdtot<-data.frame()


pts <- dipsaus::parse_svec("2")

path='/Volumes/bigbrain/Fragility2024/FragilityResultsACLNoRank/'
patient_xls <- readxl::read_xlsx("/Volumes/bigbrain/Fragility2024/FragilityEEGDataset_pipeline_update_042624.xlsx")
#pipeline_xls <- readxl::read_xlsx("/Volumes/OFZ1_T7/karaslab/rave_data/bids_dir/FragilityEEGDataset/FragilityEEGDataset_pipeline.xlsx")
patient_xls$subject[pts]


threshold_start = -10
threshold_end = 20
threshold=0.6

ploton=TRUE

irvec=c()

for(i in pts){

#i=2

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
  displayNamesoz <- channelDefs$name[soz]
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
        labs(x = "Time (s)", y = "Electrode") +
        scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5)+  #
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 5,colour=colorelec),     # Adjust depending on electrodes
        )

      file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_fragilitynorank.png',sep="")
      ggsave(file)
    }

    dum=rank(fragilitymap)
    fragrank=matrix(dum, nrow(fragilitymap), ncol(fragilitymap))

    maxfrank=max(fragrank)
    fragrank=fragrank/maxfrank
    frag_data$Value <- c(t(fragrank))

    if(ploton==TRUE){
      #plot organized fragility map with soz electrodes on bottom
      ggplot(frag_data, aes(x = Time, y = Electrode, fill = Value)) +
        geom_tile() +
        ggtitle(titlepng)+
        labs(x = "Time (s)", y = "Electrode") +
        scale_fill_viridis(option = "turbo") +
        #scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +  #
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 5,colour=colorelec),     # Adjust depending on electrodes
        )

      file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_fragilityrank.png',sep="")
      ggsave(file)
    }

    elecfm<-thresholdmean_fragility (fragrank, stimes, t_step, threshold_start, threshold_end, threshold)
    esoz<-c(1:length(elecsoz))
    startsozc=length(elecsoz)+1
    esozc<-c(startsozc:n_elec)
    tpv<-intersect(esoz,elecfm)
    tp=length(tpv)
    fn<-length(esoz)-tp
    fpv<-intersect(esozc,elecfm)
    fp=length(fpv)
    tn<-length(esozc)-fp
    ppv=tp/(tp+fp)
    npv=tn/(tn+fn)
    sensitivity=tp/(tp+fn)
    specificity=tn/(tn+fp)

    p20s<-which(stimes<=20)
    stimes20<-stimes[p20s]
    botsoz=1:nsoz
    heatmapsoz=fragrank[botsoz,p20s]
    heatmapsozc=fragrank[-botsoz,p20s]

    f90soz=quantile(heatmapsoz, probs=c(0.9))
    f90sozc=quantile(heatmapsozc,probs=c(0.9))

    interpretabilityratiosoz=f90soz/f90sozc

    irvec=rbind(irvec,interpretabilityratiosoz)

    resthreshold<-data.frame(subject_code,j,interpretabilityratiosoz,ppv,npv,specificity,sensitivity)
    resthresholdtot<-rbind(resthresholdtot,resthreshold)

    thresholdfrag=0.92


    testfrag<-vector(mode="numeric", length=n_elec)
    lengthfrag<-vector(mode="numeric", length=n_elec)
    maxfrag<-vector(mode="numeric", length=n_elec)
    startfrag<-vector(mode="numeric", length=n_elec)

    startfrag[1:n_elec]=1000000

    ltfrag=2

    for(ie in 1:n_elec){
    #  ie=1

      sige=fragrank[ie,]


      seq<-which(sige>thresholdfrag)

      resultfrag<- rle(diff(seq))
      fraglengthp<-max(resultfrag$length)*0.125

      if(fraglengthp>ltfrag){
        testfrag[ie]=1
        lengthfrag[ie]=fraglengthp
        maxfrag[ie]=max(sige)
        #startfrag[ie]=seq[1]*0.125-10
      }

    }


    dffrag=data.frame(elecnum,testfrag,lengthfrag,startfrag)

    # tpv<-intersect(esoz,elecfm)
    # tp=length(tpv)
    # fn<-length(esoz)-tp
    # fpv<-intersect(esozc,elecfm)
    # fp=length(fpv)
    # tn<-length(esozc)-fp
    # ppv=tp/(tp+fp)
    # npv=tn/(tn+fn)
    # sensitivity=tp/(tp+fn)
    # specificity=tn/(tn+fp)

    # quantilematrixsozsozc=matrix(0,20,length(stimes20))
    # cmeansoz=c(1:length(stimes20))*0
    # cmeansozc=c(1:length(stimes20))*0
    # csdsoz=c(1:length(stimes20))*0
    # csdsozc=c(1:length(stimes20))*0
    #
    #
    #
    # for(i in 1:length(stimes20)){
    #
    #   colsoz=heatmapsoz[,i]
    #   colsozc=heatmapsozc[,i]
    #
    #   meansoz=mean(colsoz)
    #   sdsoz=sd(colsoz)
    #   meansozc=mean(colsozc)
    #   sdsozc=sd(colsozc)
    #
    #   cmeansoz[i]=meansoz
    #   cmeansozc[i]=meansozc
    #   csdsoz[i]=sdsoz
    #   csdsozc[i]=sdsozc
    #
    #   f10colsoz<-quantile(colsoz,probs=c(0.1))
    #   f20colsoz<-quantile(colsoz,probs=c(0.2))
    #   f30colsoz<-quantile(colsoz,probs=c(0.3))
    #   f40colsoz<-quantile(colsoz,probs=c(0.4))
    #   f50colsoz<-quantile(colsoz,probs=c(0.5))
    #   f60colsoz<-quantile(colsoz,probs=c(0.6))
    #   f70colsoz<-quantile(colsoz,probs=c(0.7))
    #   f80colsoz<-quantile(colsoz,probs=c(0.8))
    #   f90colsoz<-quantile(colsoz,probs=c(0.9))
    #   f100colsoz<-quantile(colsoz,probs=c(1.0))
    #
    #   f10colsozc<-quantile(colsozc,probs=c(0.1))
    #   f20colsozc<-quantile(colsozc,probs=c(0.2))
    #   f30colsozc<-quantile(colsozc,probs=c(0.3))
    #   f40colsozc<-quantile(colsozc,probs=c(0.4))
    #   f50colsozc<-quantile(colsozc,probs=c(0.5))
    #   f60colsozc<-quantile(colsozc,probs=c(0.6))
    #   f70colsozc<-quantile(colsozc,probs=c(0.7))
    #   f80colsozc<-quantile(colsozc,probs=c(0.8))
    #   f90colsozc<-quantile(colsozc,probs=c(0.9))
    #   f100colsozc<-quantile(colsozc,probs=c(1.0))
    #
    #   quantilematrixsozsozc[1,i]=f10colsoz
    #   quantilematrixsozsozc[2,i]=f20colsoz
    #   quantilematrixsozsozc[3,i]=f30colsoz
    #   quantilematrixsozsozc[4,i]=f40colsoz
    #   quantilematrixsozsozc[5,i]=f50colsoz
    #   quantilematrixsozsozc[6,i]=f60colsoz
    #   quantilematrixsozsozc[7,i]=f70colsoz
    #   quantilematrixsozsozc[8,i]=f80colsoz
    #   quantilematrixsozsozc[9,i]=f90colsoz
    #   quantilematrixsozsozc[10,i]=f100colsoz
    #   quantilematrixsozsozc[11,i]=f10colsozc
    #   quantilematrixsozsozc[12,i]=f20colsozc
    #   quantilematrixsozsozc[13,i]=f30colsozc
    #   quantilematrixsozsozc[14,i]=f40colsozc
    #   quantilematrixsozsozc[15,i]=f50colsozc
    #   quantilematrixsozsozc[16,i]=f60colsozc
    #   quantilematrixsozsozc[17,i]=f70colsozc
    #   quantilematrixsozsozc[18,i]=f80colsozc
    #   quantilematrixsozsozc[19,i]=f90colsozc
    #   quantilematrixsozsozc[20,i]=f100colsozc
    #
    # }
    #
    # quantilesname<-c('SOZ(10th)','SOZ(20th)','SOZ(30th)','SOZ(40th)','SOZ(50th)',
    #                  'SOZ(60th)','SOZ(70th)','SOZ(80th)','SOZ(90th)','SOZ(100th)',
    #                  'SOZc(10th)','SOZc(20th)','SOZc(30th)','SOZc(40th)','SOZc(50th)',
    #                  'SOZc(60th)','SOZc(70th)','SOZc(80th)','SOZc(90th)','SOZc(100th)')
    # quantileplot<- expand.grid(Time = stimes20, Stats=quantilesname)
    # quantileplot$Value <- c(t(quantilematrixsozsozc))
    #
    # if(ploton==TRUE){
    #
    #   ggplot(quantileplot, aes(x = Time, y = Stats, fill = Value)) +
    #     geom_tile() +
    #     ggtitle(titlepng)+
    #     labs(x = "Time (s)", y = "Statistic") +
    #     scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +  #
    #     theme_minimal() +
    #     theme(
    #       axis.text.y = element_text(size = 10),     # Adjust depending on electrodes
    #     )
    #
    #   file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_quantiles.png',sep="")
    #   ggsave(file)
    # }
    # quantilearray<-quantileplot$Value
    # quantilesave<-rbind(quantilesave,quantilearray)
    #
    #
    # sozsdp=cmeansoz+csdsoz
    # sozsdm=cmeansoz-csdsoz
    # sozcsdp=cmeansozc+csdsozc
    # sozcsdm=cmeansozc-csdsozc
    #
    # plotmeanstd<-as.data.frame(stimes20)
    # colnames(plotmeanstd)<-"times"
    # plotmeanstd$meansoz<-cmeansoz
    # plotmeanstd$sozsdp<-sozsdp
    # plotmeanstd$sozsdm<-sozsdm
    # plotmeanstd$meansozc<-cmeansozc
    # plotmeanstd$sozcsdp<-sozcsdp
    # plotmeanstd$sozcsdm<-sozcsdm
    #
    # if(ploton==TRUE){
    #   ggplot(plotmeanstd, aes(x=times, y=cmeansoz)) +
    #     xlab('Time around seizure in s')+
    #     ylab('Fragility')+
    #     ggtitle(titlepng)+
    #     geom_line(aes(y = meansoz),color='red') +
    #     geom_line(aes(y = sozsdp),color='red',linetype="dotted") +
    #     geom_line(aes(y = sozsdm),color='red',linetype="dotted") +
    #     geom_line(aes(y = meansozc),color='black')+
    #     geom_line(aes(y = sozcsdp),color='black',linetype="dotted") +
    #     geom_line(aes(y = sozcsdm),color='black',linetype="dotted")+
    #     geom_ribbon(aes(ymin=sozsdm,ymax=sozsdp), fill="red",alpha=0.5)+
    #     geom_ribbon(aes(ymin=sozcsdm,ymax=sozcsdp), fill="black",alpha=0.5)
    #
    #   file=paste(path,subject_code,'/',subject_code,'_seizure',as.character(j),'_FragilityDistribution.png',sep="")
    #   ggsave(file)
    #}
  }



}

file=paste(path,'FragilityThresholdmean.csv',sep="")
write.csv(resthresholdtot,file)
