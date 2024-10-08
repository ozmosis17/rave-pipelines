fragility_map_plot <- function(repository, adj_frag_info, threshold_elec, display_electrodes, sz_onset, elec_list, sort_fmap = 1, height = 10, threshold) {

  m <- adj_frag_info$frag[as.character(display_electrodes),]
  elecsort <- sort(as.numeric(attr(m, "dimnames")[[1]])) # electrode indices sorted by ascending number
  fsort <- as.numeric(attr(sort(adj_frag_info$avg), "names")) # electrode indices sorted by descending fragility

  if (sort_fmap == 1) {
    elec_order <- elecsort # by electrode (ascending)
  } else if (sort_fmap == 2) {
    elec_order <- rev(elecsort) # by electrode (descending)
  } else if (sort_fmap == 3) {
    elec_order <- fsort # by fragility (ascending)
  } else if (sort_fmap == 4) {
    elec_order <- rev(fsort) # by fragility (descending)
  }

  y <- elecsort # determine what order to display electrodes in
  x <- 1:dim(m)[2]
  m <- t(m[as.character(rev(y)),]) # rev to make display descending from top to bottom

  attr(m, "xlab") = "Time (s)"
  attr(m, "ylab") = "Electrode"
  attr(m, "zlab") = "Fragility"

  tp <- repository$voltage$dimnames$Time

  if (!all(elec_list$Label == "NoLabel")) {
    elec_i <- match(elec_order, elec_list$Electrode)
    y <- paste0(elec_list$Label[elec_i], "(", elec_order, ")")
    elec_i <- match(fsort, elec_list$Electrode)
    f_list <- paste0(elec_list$Label[elec_i], "(", fsort, ")")
  }

  # for electrode label spacing on y axis
  yi = seq_along(y)
  if(length(y) > 10) {
    .seq = seq(1, length(y), length.out=height)
    y = y[.seq]
    yi = .seq
  }

  # map x axis from timewindows (x) to time (for mtext)
  xtime <- round(seq(tp[1], tp[length(tp)], length.out = 9), digits = 2)
  xi <- seq(1, length(x), length.out = 9)

  # map seizure onset from time (from slider input) to timewindows (for abline)
  secs <- seq(tp[1], tp[length(tp)])
  onset <- seq(1, length(x), length.out = length(secs))[match(sz_onset,secs)]

  # convert threshold-identified electrodes from numbers to names
  threshold_elec_i <- as.numeric(threshold_elec$elecnames)
  if (!all(elec_list$Label == "NoLabel")) {
    elec_i <- match(threshold_elec_i, elec_list$Electrode)
    threshold_elec_names <- paste0(elec_list$Label[elec_i], collapse = ", ")
  } else {
    threshold_elec_names <- dipsaus::deparse_svec(threshold_elec_i)
  }

  print(paste0("Electrodes with Fragility > ", threshold, ": ", threshold_elec_names))

  # draw fragility map
  # change color scheme
  ravebuiltins:::draw_many_heat_maps(list(
    list(
      data = m,
      x = x,
      y = seq_along(elecsort),
      has_trials = TRUE,
      range = 0:1
    )
  ), axes = c(FALSE,FALSE), PANEL.LAST = ravebuiltins:::add_decorator(function(...) {
    abline(v = onset, lty = 2, lwd = 2)
    mtext(rev(y), side=2, line=-1.5, at=yi, cex=(ravebuiltins:::rave_cex.lab*0.6), las=1)
    mtext(xtime, side=1, line=0, at=xi, cex=(ravebuiltins:::rave_cex.lab*0.6), las=1)
    mtext(paste0(repository$subject$subject_code, " Electrodes with Fragility > ", threshold, ": \n", threshold_elec_names),side = 3)
  }, ravebuiltins:::spectrogram_heatmap_decorator())
  )
}

voltage_recon_plot <- function(repository, adj_frag_info, t_window, t_step, trial_num, timepoints = 1:250, elec_num = 1, percentile = 0.1, lambda) {

  A <- adj_frag_info$adj
  S <- length(repository$voltage$dimnames$Time) # S is total number of timepoints
  N <- length(repository$voltage$dimnames$Electrode) # N is number of electrodes

  if(S %% t_step != 0) {
    # truncate S to greatest number evenly divisible by timestep
    S <- trunc(S/t_step) * t_step
  }
  n_steps <- S/t_step - (t_window/t_step) + 1 # n_steps is number of time windows

  # generate filearray for original voltage trace
  v_recon <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(S, N),
    # dimnames = c("Timepoint", "Electrode"),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature
  )

  # import voltage trace from repository$voltage$data_list
  raveio::lapply_async(repository$voltage$data_list, function(v) {
    e <- dimnames(v)$Electrode
    idx_e <- which(repository$electrode_list == e)

    v_recon[,idx_e] <- v[1:S, trial_num, 1, drop = TRUE]
    return()
  })

  signalScaling <- 10^floor(log10(max(v_recon[])))
  v_recon[] <- v_recon[]/signalScaling
  v_orig <- v_recon[]

  # populate v_recon with predicted values using adjacency matrix
  raveio::lapply_async(seq_len(n_steps), function(iw) {
    si <- seq_len(t_window-1) + (iw-1)*t_step
    xt <- v_orig[si,]
    pred_xtp1 <- xt %*% A[,,iw]
    v_recon[append(si,max(si)+1),] <- rbind(xt[1,],pred_xtp1)
    return()
  })

  # graph voltage traces for comparison

  y1 <- v_orig[timepoints,elec_num]
  y2 <- v_recon[timepoints,elec_num]
  df <- data.frame(timepoints,y1,y2)

  # check stability of adjacency matrix
  me <- raveio::lapply_async(seq_len(n_steps), function(iw) {
    AEigen <- eigen(A[,,iw])
    e <- Mod(AEigen$values)
    max(e)
  })
  print(paste0("largest eigenvalue norm: ", max(unlist(me))))

  # calculate mean squared error between y1 and y2
  mse <- mean((v_orig - v_recon[])^2)
  print(paste0("MSE: ", mse))

  R2_percentile <- mean(apply(adj_frag_info$R2,2,quantile,percentile))

  g <- ggplot2::ggplot(df, aes(timepoints)) +
    geom_line(aes(y=y1, color = "original")) +
    geom_line(aes(y=y2, color = "reconstructed")) +
    labs(x = "Time (ms)",
         y = paste0("Voltage"),
         color = "Legend",
         caption = paste0("Lambda: ", lambda, "\n",
                          percentile*100,
                          "th percentile of R2 (mean across time windows): ",
                          R2_percentile,
                          "\n Largest eigenvalue norm: ", max(unlist(me)),
                          "\n MSE: ", mse)) +
    scale_color_manual(values = c("original" = "black", "reconstructed" = "red")) +
    ggtitle(paste0("Electrode ", elec_num, " Voltage Reconstruction"))

  g
}

frag_quantile <- function(repository, f, t_window, t_step, soz, sozc){
  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)
  n_steps <- floor((n_tps - t_window) / t_step) + 1
  epoch_time_window <- repository$time_windows[[1]]
  fs <- repository$sample_rate
  if(any(repository$electrode_table$Label == "NoLabel")) {
    elec_names <- repository$electrode_table$Electrode[match(c(soz,sozc), repository$electrode_table$Electrode)]
    elec_names <- as.character(elec_names)
  } else {
    elec_names <- repository$electrode_table$Label[match(c(soz,sozc), repository$electrode_table$Electrode)]
  }

  # create fragility map with soz electrodes separated from sozc electrodes
  fmap <- f[as.character(c(soz,sozc)),]
  stimes <- (seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]

  # raw fragility map
  fplot <- expand.grid(Time = stimes, Electrode = elec_names)
  fplot$Value <- c(t(fmap))

  # create separate heatmaps for soz and sozc for quantile calcs
  hmapsoz <- fmap[as.character(soz),]
  hmapsozc <- fmap[as.character(sozc),]

  #f90soz=quantile(hmapsoz, probs=c(0.9))
  #f90sozc=quantile(hmapsozc,probs=c(0.9))
  #interpretabilityratiosoz=f90soz/f90sozc

  quantilematrixsozsozc=matrix(0,20,length(stimes))
  cmeansoz=c(1:length(stimes))*0
  cmeansozc=c(1:length(stimes))*0
  csdsoz=c(1:length(stimes))*0
  csdsozc=c(1:length(stimes))*0

  for(i in 1:length(stimes)){

    colsoz=hmapsoz[,i]
    colsozc=hmapsozc[,i]

    meansoz=mean(colsoz)
    sdsoz=sd(colsoz)
    meansozc=mean(colsozc)
    sdsozc=sd(colsozc)

    cmeansoz[i]=meansoz
    cmeansozc[i]=meansozc
    csdsoz[i]=sdsoz
    csdsozc[i]=sdsozc

    f10colsoz<-quantile(colsoz,probs=c(0.1))
    f20colsoz<-quantile(colsoz,probs=c(0.2))
    f30colsoz<-quantile(colsoz,probs=c(0.3))
    f40colsoz<-quantile(colsoz,probs=c(0.4))
    f50colsoz<-quantile(colsoz,probs=c(0.5))
    f60colsoz<-quantile(colsoz,probs=c(0.6))
    f70colsoz<-quantile(colsoz,probs=c(0.7))
    f80colsoz<-quantile(colsoz,probs=c(0.8))
    f90colsoz<-quantile(colsoz,probs=c(0.9))
    f100colsoz<-quantile(colsoz,probs=c(1.0))

    f10colsozc<-quantile(colsozc,probs=c(0.1))
    f20colsozc<-quantile(colsozc,probs=c(0.2))
    f30colsozc<-quantile(colsozc,probs=c(0.3))
    f40colsozc<-quantile(colsozc,probs=c(0.4))
    f50colsozc<-quantile(colsozc,probs=c(0.5))
    f60colsozc<-quantile(colsozc,probs=c(0.6))
    f70colsozc<-quantile(colsozc,probs=c(0.7))
    f80colsozc<-quantile(colsozc,probs=c(0.8))
    f90colsozc<-quantile(colsozc,probs=c(0.9))
    f100colsozc<-quantile(colsozc,probs=c(1.0))

    quantilematrixsozsozc[1,i]=f10colsoz
    quantilematrixsozsozc[2,i]=f20colsoz
    quantilematrixsozsozc[3,i]=f30colsoz
    quantilematrixsozsozc[4,i]=f40colsoz
    quantilematrixsozsozc[5,i]=f50colsoz
    quantilematrixsozsozc[6,i]=f60colsoz
    quantilematrixsozsozc[7,i]=f70colsoz
    quantilematrixsozsozc[8,i]=f80colsoz
    quantilematrixsozsozc[9,i]=f90colsoz
    quantilematrixsozsozc[10,i]=f100colsoz
    quantilematrixsozsozc[11,i]=f10colsozc
    quantilematrixsozsozc[12,i]=f20colsozc
    quantilematrixsozsozc[13,i]=f30colsozc
    quantilematrixsozsozc[14,i]=f40colsozc
    quantilematrixsozsozc[15,i]=f50colsozc
    quantilematrixsozsozc[16,i]=f60colsozc
    quantilematrixsozsozc[17,i]=f70colsozc
    quantilematrixsozsozc[18,i]=f80colsozc
    quantilematrixsozsozc[19,i]=f90colsozc
    quantilematrixsozsozc[20,i]=f100colsozc

  }

  quantilesname<-c("SOZ(10th)","SOZ(20th)","SOZ(30th)","SOZ(40th)","SOZ(50th)",
                   "SOZ(60th)","SOZ(70th)","SOZ(80th)","SOZ(90th)","SOZ(100th)",
                   "SOZc(10th)","SOZc(20th)","SOZc(30th)","SOZc(40th)","SOZc(50th)",
                   "SOZc(60th)","SOZc(70th)","SOZc(80th)","SOZc(90th)","SOZc(100th)")
  quantileplot<- expand.grid(Time = stimes, Stats=quantilesname)
  quantileplot$Value <- c(t(quantilematrixsozsozc))

  dimnames(quantilematrixsozsozc) <- list(
    Quantile = quantilesname,
    Time = stimes
  )

  return(list(
    fplot = fplot,
    q_matrix = quantilematrixsozsozc,
    q_plot = quantileplot
  ))
}

mean_f_plot <- function(repository, f, soz, sozc) {
  mean_f_soz <- rep(0,dim(f)[2])
  mean_f_sozc <- rep(0,dim(f)[2])
  se_f_soz <- rep(0,dim(f)[2])
  se_f_sozc <- rep(0,dim(f)[2])
  for (i in seq_len(dim(f)[2])){
    mean_f_soz[i] <- mean(f[as.character(soz),i])
    se_f_soz[i] <- sd(f[as.character(soz),i])/sqrt(length(soz))
    mean_f_sozc[i] <- mean(f[as.character(sozc),i])
    se_f_sozc[i] <- sd(f[as.character(sozc),i]/sqrt(length(sozc)))
  }
  return(list(
    mean_f_soz = mean_f_soz,
    mean_f_sozc = mean_f_sozc,
    se_f_soz = se_f_soz,
    se_f_sozc = se_f_sozc
  ))
}

output_files <- function(repository,f,pipeline_settings,export,note) {

  raveio::dir_create2(paste0(export,"/",note))

  subject_code <- pipeline_settings$subject_code
  sz_num <- pipeline_settings$condition
  t_window <- pipeline_settings$t_window
  t_step <- pipeline_settings$t_step
  soz <- pipeline_settings$soz
  sozc <- pipeline_settings$sozc
  sz_onset <- pipeline_settings$sz_onset
  epoch_time_window <- pipeline_settings$epoch_time_window

  # save fragility matrix results to csv
  raveio::safe_write_csv(
    f,
    file.path(export, paste0(note,"/",subject_code, "_", sz_num,"_fragility_",note,".csv"))
  )

  quantile_results <- frag_quantile(repository, f, t_window, t_step, soz, sozc)

  # fragility heatmap
  colorelec <- rep("black",length(c(soz,sozc)))
  colorelec[1:length(soz)]="blue"

  titlepng=paste(subject_code,as.character(sz_num),note,sep=" ")

  ggplot(quantile_results$fplot, aes(x = Time, y = Electrode, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Electrode") +
    scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5)+  #
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 5,colour=colorelec),     # Adjust depending on electrodes
    )
  img <- paste0(export,"/",note,"/",subject_code,"_",sz_num,"_map_",note,".png")
  ggsave(img)

  # quantile
  raveio::safe_write_csv(
    quantile_results$q_matrix,
    file.path(export, paste0(note,"/",subject_code, "_", sz_num,"_quantile_",note,".csv"))
  )

  # quantile map
  titlepng=paste(subject_code,as.character(sz_num),"Quantiles",note,sep=" ")

  ggplot(quantile_results$q_plot, aes(x = Time, y = Stats, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Statistic") +
    scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +  #
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),     # Adjust depending on electrodes
    )
  q_image <- paste0(export,"/",note,"/",subject_code,"_",sz_num,"_qmap_",note,".png")
  ggsave(q_image)

  # average f over time windows
  mean_f <- mean_f_plot(repository,f,soz,sozc)

  # calculate seizure onset timewindow
  sz_onset_twindow <- (sz_onset - epoch_time_window[1]) * dim(f)[2] / (epoch_time_window[2] - epoch_time_window[1])

  titlepng=paste(subject_code,"Seizure",as.character(sz_num),"Avg Fragility Over Time",note,sep=" ")

  df <- data.frame(1:length(mean_f$mean_f_soz),mean_f)

  ggplot(df, aes(1:length(mean_f_soz))) +
    geom_line(aes(y=mean_f_soz, color = "SOZ +/- sem")) +
    geom_line(aes(y=mean_f_sozc, color = "SOZc +/- sem")) +
    geom_ribbon(aes(ymin = mean_f_soz - se_f_soz, ymax = mean_f_soz + se_f_soz), fill = "indianred3", alpha = 0.7) +
    geom_ribbon(aes(ymin = mean_f_sozc - se_f_sozc, ymax = mean_f_sozc + se_f_sozc), fill = "grey30", alpha = 0.6) +
    geom_vline(xintercept = sz_onset_twindow, color = "blue") +
    scale_color_manual(values = c("SOZ +/- sem" = "red", "SOZc +/- sem" = "black")) +
    labs(x = "Timewindow", y = "Average fragility per timewindow", color = "Legend") +
    ggtitle(titlepng)

  mean_plot <- paste0(export,"/",note,"/",subject_code,"_",sz_num,"_meanplot_",note,".png")
  ggsave(mean_plot)

  line_plot_df <- data.frame(mean_f,attr(quantile_results$q_matrix, "dimnames")$Time)
  names(line_plot_df)[5] <- "time"

  raveio::safe_write_csv(
    line_plot_df,
    file.path(export, paste0(note,"/",subject_code, "_", sz_num,"_meandata_",note,".csv"))
  )
}

# y1 <- repository$voltage$data_list$e_1[]
# y2 <- Fragrepository$voltage$data_list$e_1[]
# y1 <- y1[1:1000,2]
# y2 <- y2[1:1000,1]
# df <- data.frame(timepoints,y1,y2)
#
# ggplot2::ggplot(df, aes(timepoints)) +
#   geom_line(aes(y=y1, color = "retrostudy")) +
#   geom_line(aes(y=y2, color = "frag")) +
#   scale_color_manual(values = c("retrostudy" = "black", "frag" = "red"))

export_pdf <- function(expr, path, env = parent.frame(),
                       quoted = FALSE, width = 12, height = 7, useDingbats = FALSE, ...) {
  force(path)
  if(!quoted) {
    expr <- substitute(expr)
  }
  grDevices::pdf(path, width = width, height = height, useDingbats = useDingbats, ...)
  on.exit({
    grDevices::dev.off()
  }, add = TRUE, after = TRUE)
  eval(expr, envir = env)
}
