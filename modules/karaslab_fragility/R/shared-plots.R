fragility_map_plot <- function(repository, adj_frag_info, threshold_elec, display_electrodes, sz_onset, elec_list, sort_fmap = 1, height = 10) {

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

  attr(m, 'xlab') = 'Time (s)'
  attr(m, 'ylab') = 'Electrode'
  attr(m, 'zlab') = 'Fragility'

  tp <- repository$voltage$dimnames$Time

  if (!all(elec_list$Label == 'NoLabel')) {
    elec_i <- match(elec_order, elec_list$Electrode)
    y <- paste0(elec_list$Label[elec_i], '(', elec_order, ')')
    elec_i <- match(fsort, elec_list$Electrode)
    f_list <- paste0(elec_list$Label[elec_i], '(', fsort, ')')
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
  if (!all(elec_list$Label == 'NoLabel')) {
    elec_i <- match(threshold_elec_i, elec_list$Electrode)
    threshold_elec_names <- paste0(elec_list$Label[elec_i], collapse = ", ")
  } else {
    threshold_elec_names <- dipsaus::deparse_svec(threshold_elec_i)
  }

  print(paste0("Electrodes with fragility above specified threshold: ", threshold_elec_names))

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
    mtext(paste0(repository$subject$subject_code, " High Fragility Electrodes: \n", threshold_elec_names),side = 3)
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
  print(paste0('largest eigenvalue norm: ', max(unlist(me))))

  # calculate mean squared error between y1 and y2
  mse <- mean((v_orig - v_recon[])^2)
  print(paste0('MSE: ', mse))

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
    scale_color_manual(values = c("original" = 'black', "reconstructed" = "red")) +
    ggtitle(paste0("Electrode ", elec_num, " Voltage Reconstruction"))

  g
}


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
