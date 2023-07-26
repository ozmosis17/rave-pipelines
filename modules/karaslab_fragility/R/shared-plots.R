fragility_map_plot <- function(repository, f_info, displayed_elec, sz_onset, elec_list, sort_fmap = 1, height = 10) {

  m <- f_info$norm[as.character(displayed_elec),]
  elecsort <- sort(as.numeric(attr(m, "dimnames")[[1]])) # electrode indices sorted by ascending number
  fsort <- as.numeric(attr(sort(f_info$avg), "names")) # electrode indices sorted by descending fragility

  if (sort_fmap == 1) {
    elec_order <- elecsort # by electrode (ascending)
  } else if (sort_fmap == 2) {
    elec_order <- rev(elecsort) # by electrode (descending)
  } else if (sort_fmap == 3) {
    elec_order <- fsort # by fragility (ascending)
  } else if (sort_fmap == 4) {
    elec_order <- rev(fsort) # by fragility (descending)
  }

  y <- rev(elecsort) # determine what order to display electrodes in
  x <- 1:dim(m)[2]
  m <- t(m[as.character(y),])

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

  # draw fragility map
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
    mtext(y, side=2, line=-1.5, at=yi, cex=(ravebuiltins:::rave_cex.lab*0.6), las=1)
    mtext(xtime, side=1, line=0, at=xi, cex=(ravebuiltins:::rave_cex.lab*0.6), las=1)
  }, ravebuiltins:::spectrogram_heatmap_decorator())
  )
}

voltage_recon_plot <- function(repository, A, t_window, t_step, trial_num, timepoints = 1:200, elec_num = 1) {

  S <- length(repository$voltage$dimnames$Time) # S is total number of timepoints
  N <- length(repository$voltage$dimnames$Electrode) # N is number of electrodes

  if(S %% t_step != 0) {
    # truncate S to greatest number evenly divisible by timestep
    S <- trunc(S/t_step) * t_step
  }
  n_steps <- S/t_step - (t_window/t_step) + 1 # n_steps is number of time windows

  # generate filearray for original voltage trace
  v_orig <- filearray::filearray_load_or_create(
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

    v_orig[,idx_e] <- v[1:S, trial_num, 1, drop = TRUE]
    return()
  })

  # generate filearray for reconstructed voltage trace
  v_recon <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(t_window, n_steps, N),
    # dimnames = c("Time", "Step", "Electrode"),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature
  )

  # split data into timewindows by timestep
  raveio::lapply_async(seq_len(N), function(e) {
    trial_voltage <- v_orig[,e]

    idx_e <- e == 1:N
    idx <- seq_len(t_window)
    lapply(seq_len(n_steps), function(step) {
      t_start <- 1 + (step - 1) * t_step
      v_recon[, step, idx_e] <- trial_voltage[t_start + idx - 1]
    })
    return()
  })

  # use adjacency array A to reconstruct timewindows
  raveio::lapply_async(
    seq_len(n_steps), function(step) {
      slice <- v_recon[, step, , drop = FALSE, dimnames = NULL]
      dm <- dim(slice)
      nr <- nrow(slice)
      dim(slice) <- c(nr, dm[[3]])
      x <- slice[-nr, , drop = FALSE]
      y <- A[,,step] %*% t(x)
      v_recon[,step,] <- t(cbind(slice[1,],y))
      return()
    }
  )

  # convert v_recon to S (timepoints) by N (electrodes)
  v_reconstructed <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(S,N),
    # dimnames = c("Time", "Electrode"),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature
  )

  # populate
  raveio::lapply_async(seq_len(N), function(idx_e) {
    trial_voltage <- v_recon[,,idx_e]
    lapply(seq_len(n_steps), function(step) {
      t_start <- 1 + (step - 1) * t_step
      v_reconstructed[t_start:(t_start+t_window-1),idx_e] <- trial_voltage[1:t_window,step]
    })
    return()
  })

  # graph voltage traces for comparison

  y1 <- v_orig[timepoints,elec_num]
  y2 <- v_reconstructed[timepoints,elec_num]
  df <- data.frame(timepoints,y1,y2)

  # check stability of adjacency matrix
  eigv <- abs(eigen(A[,,1], only.values = TRUE)$values)
  print(paste0('largest eigenvalue norm: ', max(eigv)))

  # calculate mean squared error between y1 and y2
  mse <- mean((v_orig[] - v_reconstructed[])^2)
  print(paste0('MSE: ', mse))

  g <- ggplot(df, aes(timepoints)) +
    geom_line(aes(y=y1, color = "original")) +
    geom_line(aes(y=y2, color = "reconstructed")) +
    labs(x = "Time (ms)", y = paste0("Voltage - Electrode ", elec_num), color = "Legend") +
    scale_color_manual(values = c("original" = 'black', "reconstructed" = "red")) +
    ggtitle(paste0("Mean Squared Error: ", format(mse, scientific = TRUE)))

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
