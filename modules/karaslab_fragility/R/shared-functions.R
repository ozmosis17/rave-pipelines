# Put shared functions in `shared-*.R` so the pipeline is clean

get_default_cores <- function(round = TRUE) {
  re <- (raveio::raveio_getopt("max_worker") + 1) / 2
  if( round ) {
    re <- ceiling(re)
  }
  re
}

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

  # scale fragility values from -1 to 1 with 1 being most fragile

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
