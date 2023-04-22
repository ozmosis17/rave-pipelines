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
# raveio::save_yaml(
#   list(
#     project_name = "OnsetZone",
#     subject_code = "PT01",
#     epoch_name = "PT01_sz",
#     epoch_time_window = c(-20, 20),
#     reference_name = "car",
#     load_electrodes = "1:24,26:36,42:43,46:54,56:70,72:95",
#     analyze_electrodes = "14",
#     t_window = 100,
#     t_step = 100,
#     nlambda = 16,
#     ncores = NA,
#     trial_num = NULL
#   ),
#   # Save to module settings.yaml for debug use
#   file = file.path(dipsaus::rs_active_project(), "modules",
#                    "karaslab_fragility", "settings.yaml"))
#
# # Loads inputs & shared functions
# raveio::pipeline_setup_rmd("karaslab_fragility", env = globalenv())

# ---- Global, will not be in the pipeline --------------------------------

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

v <- repository$voltage

# obtain voltage data
# volt <- module_tools$get_voltage()
# v <- volt$get_data()

# voltage data are stored at `repository$voltage`

# use voltage data to calculate adjacency array
A <- generate_adj_array(t_window, t_step, v, trial_num, nlambda, ncores)

# use adjacency array to find f_info
f_info <- generate_fragility_matrix(
  A = A,
  elec = attr(v, "dimnames")$Electrode,
  ncores = ncores
)

generate_state_vectors <- function(t_start,vmat,t_window) {

  state_vectors <- list(
    # x(t)
    x = vmat[,t_start:(t_start+t_window-2)],

    # x(t+1)
    x_n = vmat[,(t_start+1):(t_start+t_window-1)]
  )
  return(state_vectors)
}

generate_adj_array <- function(t_window, t_step, v, trial_num, nlambda, ncores) {
  print('Generating adjacency array')

  S <- length(v$dimnames$Time) # S is total number of timepoints
  N <- length(v$dimnames$Electrode) # N is number of electrodes

  if(S %% t_step != 0) {
    # truncate S to greatest number evenly divisible by timestep
    S <- trunc(S/t_step) * t_step
  }

  J <- S/t_step - (t_window/t_step) + 1 # J is number of time windows

  # A will be the adjacency array, contains J adjacency matrices (one per time window)
  A <- array(dim = c(N,N,J))

  # populate adjacency array
  raveio::lapply_async(1:J, function(k) {
    # k is timewindow index, from 1 to J
    t_start <- 1+(k-1)*t_step # calc timepoints within kth time window

    # create N by timepoints matrix of voltage values for trial_num trial
    vmat <- matrix(nrow = N, ncol = length(v$dimnames$Time))
    for (n in 1:v$dim[3]) {
      e_name <- paste0('e_',v$dimnames$Electrode[n])
      vmat[n,] <- v$data_list[[e_name]][,trial_num,1]
    }

    # find state vectors x(t) and x(t+1) for kth timewindow
    svec <- generate_state_vectors(t_start, vmat, t_window)

    # calculate adjacency matrix for kth timewindow
    A_list <- find_adj_matrix(svec,N, t_window, nlambda)

    # convert into matrix and insert into kth slot of adjacency array
    A[,,k] <- matrix(unlist(A_list), nrow = N, ncol = N)
  }, callback = function(k) {
    sprintf("Calculating...|Processing timewindow %s", k)
  })

  return(A)
}

find_adj_matrix <- function(state_vectors, N, t_window, nlambda) {
  # vectorize x(t+1)
  b <- c(state_vectors$x_n)

  # initialize big H matrix for system of linear equations
  H <- matrix(0, nrow = N*(t_window-1), ncol = N^2)

  # populate H matrix
  r <- 1
  for (ii in 1:(t_window-1)) {
    c <- 1
    for (jj in 1:N) {
      H[r,c:(c+N-1)] <- state_vectors$x[,ii]
      c <- c + N
      r <- r + 1
    }
  }

  # solve system using glmnet package least squares, with L2-norm regularization
  # aka ridge filtering

  # find optimal lambda
  cv.ridge <- glmnet::cv.glmnet(H, b, alpha = 0, nfolds = 3, parallel = FALSE, nlambda = nlambda)
  lambdas <- rev(cv.ridge$lambda)

  test_lambda <- function(l, H, b) {
    ridge <- glmnet::glmnet(H, b, alpha = 0, lambda = l)
    N <- sqrt(dim(H)[2])
    adj_matrix <-  matrix(ridge$beta, nrow = N, ncol = N, byrow = TRUE)
    eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
    stable <- max(eigv) < 1
    list(
      adj = adj_matrix,
      abs_eigv = eigv,
      stable = stable
    )
  }

  l <- 1
  stable_i <- FALSE

  while (!stable_i) {
    results <- test_lambda(lambdas[l], H = H, b = b)
    stable_i <- results$stable

    l <- l + 1

    if (l > length(lambdas)) {
      break
    }
  }

  if (!stable_i) {
    stop('No lambdas result in a stable adjacency matrix. Increase the number of lambdas, or (more likely) there is something wrong with your data.')
  }

  adj_matrix <- results$adj

  return(adj_matrix)
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

  N <- dim(A)[1]
  J <- dim(A)[3]
  f_vals <- matrix(nrow = N, ncol = J)
  fprogress = rave::progress(title = 'Generating Fragility Matrix (Step 2 of 2)', max = J)
  shiny::showNotification('Calculating estimated time remaining...', id = 'first_est', duration = NULL)

  for (k in 1:J) {
    start_time <- Sys.time()
    print(paste0('Current timewindow: ', k, ' out of ', J))
    fprogress$inc(paste0('Current timewindow: ', k, ' out of ', J))

    for (i in seq(1,N,ncores)) {
      if (i+ncores-1 <= N) {
        is <- i:(i+ncores-1)
        f_vals_list <- rave::lapply_async3(is,find_fragility,A_k = A[,,k], N = N, limit = lim)
        f_vals[is,k] <- unlist(f_vals_list)
      } else {
        is <- i:N
        f_vals_list <- rave::lapply_async3(is,find_fragility,A_k = A[,,k], N = N, limit = lim)
        f_vals[is,k] <- unlist(f_vals_list)
      }
    }

    end_time <- Sys.time()
    print(end_time - start_time)

    if (k == 1) {
      shiny::removeNotification(id = 'first_est')
      t_avg <- 0
    }

    t_avg <- (t_avg*(k-1) + as.numeric(difftime(end_time, start_time, units='sec')))/k
    shiny::showNotification(paste0('Estimated time remaining: ', (t_avg*(J-k))%/%60, ' minutes'), id = 'est_time', duration = NULL)
  }

  shiny::removeNotification(id = 'est_time')
  fprogress$close()
  rownames(f_vals) <- elec
  colnames(f_vals) <- 1:J

  f_norm <- f_vals

  # # scale fragility values from 0 to 1 with 1 being most fragile
  # for (j in 1:J) {
  #   max_f <- max(f_vals[,j])
  #   f_norm[,j] <- sapply(f_vals[,j], function(x) (max_f - x) / max_f)
  # }

  # scale fragility values from -1 to 1 with 1 being most fragile
  for (j in 1:J) {
    max_f <- max(f_vals[,j])
    f_norm[,j] <- sapply(f_vals[,j], function(x) 2*(max_f - x)/max_f - 1)
  }

  # find average fragility for each electrode across time
  avg_f <- rowMeans(f_norm)

  f_info <- list(
    vals = f_vals,
    norm = f_norm,
    avg = avg_f
  )
}
