# # inputs required:
# project_name
# subject_code
# selected electrodes (default to preprocessed electrodes?),
# time window size (t_window, default=100)
# time step (t_step, default=100)
# # of lambdas (nlambda, default=16)
# # of cores (ncores, default=see calculation on lines 17-18)
# trial selection (trial_num)

# outputs: f_info which is a list(
# vals = raw fragility values
# norm = fragility values normalized from -1 to 1
# avg = average f value for an electrode across time
# )

# ---- preparation, will not be in the pipeline --------------------------------
# this script require `glmnet`
if(!dipsaus::package_installed("glmnet")) {
  ravemanager:::install_packages("glmnet")
}
# Enable parallel for debugging
options("raveio.auto.parallel" = TRUE)

# Set inputs for debug. In the following context, you can use these variables
# as if they exist.
# Please only enter simple inputs (e.g. strings, numbers, list, ...). Do NOT
# put functions/objects/data/environments here (there are places for that)

# DIPSAUS DEBUG START
# if(Sys.info()["login"] == "dipterix") {
#   project_name <- "demo"
#   subject_code <- "DemoSubject"
#   epoch_name <- "auditory_onset"
#   epoch_time_window <- c(-1,2)
#   reference_name <- "default"
# } else {
#   project_name = "OnsetZone"
#   subject_code = "PT01"
#   epoch_name = "seizure_onset"
#   epoch_time_window = c(-10, 10)
#   reference_name = "car"
# }
# raveio::save_yaml(
#   list(
#     project_name = project_name,
#     subject_code = subject_code,
#     epoch_name = epoch_name,
#     epoch_time_window = epoch_time_window,
#     reference_name = reference_name,
#     load_electrodes = "1:24,26:36,42:43,46:54,56:70,72:95",
#     t_window = 250,
#     t_step = 125,
#     nlambda = 16,
#     ncores = NA,
#     trial_num = 1,
#     sz_onset = 0
#   ),
#   # Save to module settings.yaml for debug use
#   file = file.path(dipsaus::rs_active_project(), "modules",
#                    "karaslab_fragility", "settings.yaml"))
#
# # Loads inputs & shared functions
# raveio::pipeline_setup_rmd("karaslab_fragility", env = globalenv())

# ---- Functions to be moved to shared-functions.R ----------------------------


find_adj_matrix <- function(x, y, nlambda) {

  # x <- t(state_vectors$x)
  # y <- t(state_vectors$x_n)

  nobs <- nrow(x)
  nvars <- ncol(x)

  # scale x and y
  x <- scale(x)
  y <- scale(y)
  scale_x <- attr(x, "scaled:scale")
  scale_y <- attr(y, "scaled:scale")

  # guess possible lambdas
  lambdas <- rev(lambda_path(
    x = x, y = y, nlambda = nlambda, alpha = 0.0,
    standardize = TRUE, intercept = FALSE))

  XtY <- crossprod(x, y)
  XtX <- crossprod(x)

  # SVD decomposition to speed things up
  # XtX == svd$u %*% diag(svd$d) %*% t(svd$v)
  svd <- svd(XtX)

  # solve(XtX) == svd$v %*% diag(1 / svd$d) %*% t(svd$u)
  V <- svd$v
  Ut <- t(svd$u)

  # sanity check: the following should be identity matrix
  # (XtX + diag(12, nvars)) %*% (V %*% diag(1 / (svd$d + 12)) %*% Ut)
  UtXtY <- Ut %*% XtY

  # for each, lambda, fit ridge regression
  dipsaus::forelse(
    x = lambdas,
    FUN = function(lam) {

      # The following 3 methods generate similar/same results
      # MASS::lm.ridge(y[,2] ~ x - 1, lambda = lam * nobs)
      # glmnet::glmnet(y = y[,1], x = x, intercept = FALSE, lambda = lam, alpha = 0.0)$beta

      # no SVD method
      # ident <- diag(as.double(nobs), ncol(x))
      # adj_matrix <- solve(XtX + lam * ident) %*% XtY

      # SVD method
      adj_matrix <- V %*% diag(1 / (svd$d + lam * nobs)) %*% UtXtY

      # right now scale(x) %*% adj_matrix = scale(y), need to scale back
      adj_matrix <- diag( 1 / scale_x ) %*% adj_matrix %*% diag( scale_y )


      eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
      stable <- max(eigv) < 1

      if( !stable ) { return() }
      # return(list(
      #   adj = adj_matrix,
      #   abs_eigv = eigv,
      #   stable = stable
      # ))

      structure(adj_matrix)
      return(adj_matrix)

    },
    ALT_FUN = function() {
      stop('No lambdas result in a stable adjacency matrix. Increase the number of lambdas, or (more likely) there is something wrong with your data.')
    }
  )

}

# generate_state_vectors <- function(vmat, t_start, t_window) {
#
#   idx <- seq_len(t_window - 1)
#
#   state_vectors <- list(
#     # x(t)
#     x = t(vmat[t_start + idx, ]),
#
#     # x(t+1)
#     x_n = t(vmat[(t_start+1) + idx, ])
#   )
#   return(state_vectors)
# }

generate_adjacency_array <- function(repository, trial_num, t_window, t_step, nlambda) {
  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)

  # Number of steps
  n_steps <- floor((n_tps - t_window) / t_step) + 1

  # slice of data
  arr <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(t_window, n_steps, n_elec),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature,
    t_step = t_step, t_window = t_window,
    trial_num = trial_num,

    on_missing = function(arr) {
      arr$set_header("ready", value = FALSE)
    }
  )

  # check if header `ready` is not TRUE
  if(!isTRUE(arr$get_header("ready", FALSE))) {

    loaded_electrodes <- repository$electrode_list
    raveio::lapply_async(repository$voltage$data_list, function(v) {
      e <- dimnames(v)$Electrode
      idx_e <- loaded_electrodes == e

      trial_voltage <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]

      idx <- seq_len(t_window)
      lapply(seq_len(n_steps), function(step) {
        t_start <- 1 + (step - 1) * t_step
        arr[, step, idx_e] <- trial_voltage[t_start + idx]
      })
      return()
    })

  }

  # calculate adjacency arrays
  A <- raveio::lapply_async(
    seq_len(n_steps), function(step) {
      slice <- arr[, step, , drop = FALSE, dimnames = NULL]
      dm <- dim(slice)
      nr <- nrow(slice)
      dim(slice) <- c(nr, dm[[3]])
      # x is x(t) and y is x(t+1), state vectors
      x <- slice[-nr, , drop = FALSE]
      y <- slice[-1, , drop = FALSE]
      as.vector(find_adj_matrix(x = x, y = y, nlambda = nlambda))
      # as.vector(find_adj_matrix_bij(x = x, y = y))
    }
  )
  A <- do.call("cbind", A)
  dim(A) <- c(n_elec, n_elec, n_steps)
  A
}


find_fragility <- function(node, A_k, N, limit) {

  e_k <- vector(mode = 'numeric', length = N)
  e_k[node] <- 1

  argument <- t(e_k) %*% (solve(A_k - limit*diag(N))) # column perturbation
  # argument <- t(e_k) %*% t(solve(A_k - num*diag(N))) # row perturbation

  B <- rbind(Im(argument),Re(argument))

  perturb_mat <- (t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) %*% t(e_k) # column
  # perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row

  norm(perturb_mat, type = '2')
}

generate_fragility_matrix <- function(A, elec, lim = 1i, ncores) {
  print('Generating fragility matrix')

  dm <- dim(A)
  N <- dm[1]
  J <- dm[3]

  f_vals <- raveio::lapply_async(
    seq_len(J),
    function(k) {
      A_k <- A[,, k]
      f_vals_k <- vapply(seq_len(N), function(i){
        find_fragility(i, A_k = A_k, N = N, limit = lim)
      }, FUN.VALUE = 0.0)
      f_vals_k
    }
  )
  # N x J Fragility Matrix
  f_vals <- do.call("cbind", f_vals)

  dimnames(f_vals) <- list(
    Electrode = elec,
    Step = seq_len(J)
  )

  # # scale fragility values from 0 to 1 with 1 being most fragile
  # for (j in 1:J) {
  #   max_f <- max(f_vals[,j])
  #   f_norm[,j] <- sapply(f_vals[,j], function(x) (max_f - x) / max_f)
  # }

  # scale fragility values from -1 to 1 with 1 being most fragile

  # max_f <- max(f_vals)
  # min_f <- min(f_vals)
  # normalize, for each column (margin=2L)
  f_norm <- apply(f_vals, 2, function(f_col) {
    max_f <- max(f_col)
    min_f <- min(f_col)
    2.0 * (f_col - min_f) / (max_f - min_f) - 1.0 # normalize from -1 to 1
    #(f_col - min_f) / (max_f - min_f) # normalize from 0 to 1
  })

  # find average fragility for each electrode across time
  avg_f <- rowMeans(f_norm)

  return(list(
    vals = f_vals,
    norm = f_norm,
    avg = avg_f
  ))
}

draw_heatmap <- function(f_info,requested_electrodes) {
  f_info$norm <- f_info$norm[as.character(requested_electrodes),]
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
  sz_onset <- 0
  secs <- seq(tp[1], tp[length(tp)])
  onset <- seq(1, length(x), length.out = length(secs))[match(sz_onset,secs)]

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

voltage_plots <- function(repository,A,timepoints,elec_num) {
  require(ggplot2)

  S <- length(repository$voltage$dimnames$Time) # S is total number of timepoints
  N <- length(repository$voltage$dimnames$Electrode) # N is number of electrodes

  if(S %% t_step != 0) {
    # truncate S to greatest number evenly divisible by timestep
    S <- trunc(S/t_step) * t_step
  }
  n_steps <- S/t_step - (t_window/t_step) + 1 # J is number of time windows

  # generate matrix for reconstruction of voltage trace using x(t+1) = Ax(t)
  v_recon <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(t_window, n_steps, N),
    type = "float", mode = "readwrite", partition_size = 1L,

    # if repository has changed, re-calculate
    repository_signature = repository$signature
  )

  raveio::lapply_async(repository$voltage$data_list, function(v) {
    e <- dimnames(v)$Electrode
    idx_e <- repository$electrode_list == e

    trial_voltage <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]

    idx <- seq_len(t_window)
    lapply(seq_len(n_steps), function(step) {
      t_start <- 1 + (step - 1) * t_step
      v_recon[, step, idx_e] <- trial_voltage[t_start + idx - 1]
    })
    return()
  })

  v_trace <- v_recon[]
  dim(v_trace) <- c(t_window*n_steps,N)

  # populate v_recon using adjacency matrix A
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

  v_reconstructed <- v_recon[]
  dim(v_reconstructed) <- c(t_window*n_steps,N)

  # graph voltage traces for comparison

  y1 <- v_trace[timepoints,elec_num]
  y2 <- v_reconstructed[timepoints,elec_num]
  df <- data.frame(timepoints,y1,y2)

  g <- ggplot(df, aes(timepoints)) +
    geom_line(aes(y=y1, color = "original")) +
    geom_line(aes(y=y2, color = "reconstructed")) +
    labs(x = "Time (ms)", y = "Voltage", color = "Legend") +
    scale_color_manual(values = c("original" = 'black', "reconstructed" = "red"))


  #plot(x = timepoints, y = y1, type = 'l', main = 'original)
  #plot(x = timepoints, y = y2, type = 'l', main = 'reconstructed')

  # check stability of adjacency matrix
  eigv <- abs(eigen(A[,,1], only.values = TRUE)$values)
  print(paste0('largest eigenvalue norm: ', max(eigv)))
  g
}

# ---- Analysis script --------------------------------------------------------

# Load subject instance
subject <- raveio::RAVESubject$new(project_name = project_name,
                                   subject_code = subject_code,
                                   strict = TRUE)

# Check data to load
loading_elec <- dipsaus::parse_svec(load_electrodes)
loading_elec <- subject$electrodes[subject$electrodes %in% loading_elec]
if(!length(loading_elec)) {
  stop("No valid electrode to load!")
}

# load voltage repository.
repository <- raveio::prepare_subject_voltage_with_epoch(
  subject = subject,
  electrodes = loading_elec,
  epoch_name = epoch_name,
  reference_name = reference_name,
  time_windows = epoch_time_window
)

# voltage data are stored at `repository$voltage`

# use voltage data to calculate adjacency array
A <- generate_adjacency_array(
  repository = repository,
  trial_num = trial_num,
  t_window = t_window,
  t_step = t_step,
  nlambda = nlambda
)

# use adjacency array to find fragility
f_info <- generate_fragility_matrix(
  A = A,
  elec = repository$electrode_list,
  ncores = 4
)

# verify fragility function: load previously generated adjacency matrix
# adj_info_og <- readRDS("/Volumes/OFZ1_T7/karaslab/rave_data/data_dir/OnsetZone/PT01/rave/module_data/PT01_adj_info_trial_1")

# f_info_check <- generate_fragility_matrix(
#   A = adj_info_og$A,
#   elec = repository$electrode_list,
#   ncores = 4
# )

# f_info_og <- readRDS("/Volumes/OFZ1_T7/karaslab/rave_data/data_dir/OnsetZone/PT01/rave/module_data/PT01_f_info_trial_1")

# TEST FRAGILITY BY CREATING HEATMAP
# loading_elec <- c(33,34,62:69,88:91)

draw_heatmap(f_info,loading_elec)
# draw_heatmap(f_info_check,loading_elec)
# draw_heatmap(f_info_og,loading_elec)

as.numeric(attr(sort(f_info$avg), "names"))
# as.numeric(attr(sort(f_info_og$avg), "names"))

# TEST RECONSTRUCTION USING ADJACENCY MATRIX
voltage_plots(repository,A,timepoints = 1:500, elec_num = 1)
