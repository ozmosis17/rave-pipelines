fragility_map_timestamp_plot <- function(repository, f_info, loading_elec, sz_onset) {
  # TODO: make this a function
  f_info$norm <- f_info$norm[as.character(loading_elec),]
  elecsort <- sort(as.numeric(attr(f_info$norm, "dimnames")[[1]])) # electrode indices sorted by ascending number
  fsort <- as.numeric(attr(sort(f_info$avg), "names")) # electrode indices sorted by descending fragility

  # if (sort_fmap == 'Electrode (ascending)') {
  #   elec_order <- elecsort
  # } else if (sort_fmap == 'Electrode (descending)') {
  #   elec_order <- rev(elecsort)
  # } else if (sort_fmap == 'Fragility (ascending)') {
  #   elec_order <- fsort
  # } else if (sort_fmap == 'Fragility (descending)') {
  #   elec_order <- rev(fsort)
  # }

  y <- rev(elecsort) # determine what order to display electrodes in

  f_info$norm <- f_info$norm[as.character(y),]
  x <- 1:dim(f_info$norm)[2]
  m <- t(f_info$norm)

  attr(m, 'xlab') = 'Time (s)'
  attr(m, 'ylab') = 'Electrode'
  attr(m, 'zlab') = 'Fragility'

  tp <- repository$voltage$dimnames$Time

  # for electrode label spacing on y axis
  yi = seq_along(y)
  if(length(y) > 10) {
    .seq = seq(1, length(y), length.out=10)
    y = y[.seq]
    yi = .seq
  }

  # map x axis from timewindows (x) to time (for mtext)
  xtime <- round(seq(tp[1], tp[length(tp)], length.out = 9), digits = 2)
  xi <- seq(1, length(x), length.out = 9)

  # map seizure onset from time (from slider input) to timewindows (for abline)
  secs <- seq(tp[1], tp[length(tp)])
  onset <- seq(1, length(x), length.out = length(secs))[match(sz_onset,secs)]

  plot_fragility_map_timestamp <- Sys.time()

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
    mtext(y, side=2, line=-1, at=yi, cex=(ravebuiltins:::rave_cex.lab*0.8), las=1)
    mtext(xtime, side=1, line=1, at=xi, cex=(ravebuiltins:::rave_cex.lab*0.8), las=1)
  }, ravebuiltins:::spectrogram_heatmap_decorator())
  )
}
