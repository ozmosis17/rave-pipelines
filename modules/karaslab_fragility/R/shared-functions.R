# Put shared functions in `shared-*.R` so the pipeline is clean

get_default_cores <- function(round = TRUE) {
  re <- (raveio::raveio_getopt("max_worker") + 1) / 2
  if( round ) {
    re <- ceiling(re)
  }
  re
}

ridge <- function(xt, xtp1, lambda, intercept = FALSE, iw) {
  if (!identical(dim(xt), dim(xtp1))) {
    stop("Unmatched dimension")
  }
  nel <- ncol(xt)
  ## Coefficient matrix A
  ## each column is coefficients from a linear regression
  ## formula: xtp1 = xt*A + E
  A <- matrix(0, nel + intercept, nel)
  ## for each electrode
  for (i in seq_len(nel)) {
    y <- xtp1[, i]
    fit <- glmnet::glmnet(xt, y,
                  alpha = 0, lambda = lambda,
                  standardize = FALSE, intercept = intercept
    )
    # fit <- glmnet::cv.glmnet(xt, y,
    #               alpha = 0,
    #               standardize = FALSE, intercept = intercept
    # )

    if (intercept) {
      A[, i] <- as.numeric(coef(fit))
    } else {
      A[, i] <- coef(fit)[-1]
    }
  }
  AEigen <- eigen(A)
  e <- Mod(AEigen$values)
  me <- max(e)
  if(me >= 1){
    print(paste0("Solution in timewindow ",iw," is not stable (",me,")"))
  } else {
    print(paste0("Solution in timewindow ",iw," is stable (",me,")"))
  }
  A
}

ridgeR2 <- function(xt, xtp1, A) {
  nel <- ncol(xt)
  ypredMat <- predictRidge(xt, A)

  R2 <- rep(0, nel)
  for (i in seq_len(nel)) {
    y <- xtp1[, i]
    ypred <- ypredMat[, i]
    sst <- sum((y - mean(y))^2)
    sse <- sum((ypred - y)^2)
    rsq <- 1 - sse / sst
    R2[i] <- rsq
  }
  R2
}

fragilityRow <- function(A, nSearch = 100) {
  ## The adjacency matrix A here is a transpose of the
  ## adjacency matrix in the original paper
  nel <- ncol(A)
  e <- Mod(eigen(A)$values)
  me <- max(e)

  if (me >= 1) {
    #return(0)
  }


  fragcol <- matrix(0, nel, nel)
  fragNorm <- rep(0, nel)
  omvec <- seq(0, 1, length.out = nSearch + 1)[-1]

  b <- c(0, -1)
  ## for each electrode
  for (i in 1:nel) {
    ## indicate which electrode is disturbed (ith)
    ek <- rep(0, nel)
    ek[i] <- 1
    tek <- t(ek)
    minNorm <- 100000
    minPerturbColumn <- NA
    for (k in seq_len(nSearch)) {
      ## imaginary part
      om <- omvec[k]
      ## real part
      sigma <- sqrt(1 - om^2)
      ## target eigenvalue
      lambda <- complex(real = sigma, imaginary = om)
      ## A - (sigma + j* omega)*I
      mat <- A - lambda * diag(nel)
      imat <- t(solve(mat))

      argument <- tek %*% imat
      B <- rbind(Im(argument), Re(argument))
      ## B^T*(B*B^T)^-1*b
      invBtB <- solve(B %*% t(B))
      prov <- t(B) %*% invBtB %*% b

      sigma_hat <- ek %*% t(prov)

      ## validation
      if (FALSE) {
        A2 <- A + sigma_hat
        e2 <- eigen(A2)$values
        closestIndex <- which.min(abs(e2 - lambda))
        e2[closestIndex]
      }

      norm_sigma_hat <- norm(sigma_hat, type = "2")
      if (norm_sigma_hat < minNorm) {
        minPerturbColumn <- prov
        minNorm <- norm(prov, type = "2")
      }
    }

    fragcol[, i] <- minPerturbColumn
    fragNorm[i] <- minNorm
  }

  maxf <- max(fragNorm)
  fragNorm2 <- (maxf - fragNorm) / maxf

  return(fragNorm2)
}

calc_adj_frag <- function(repository, trial_num, t_window, t_step, lambda = 0.0001) {

  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)

  # Number of steps
  n_steps <- floor((n_tps - t_window) / t_step) + 1

  # slice of data
  arr <- filearray::filearray_load_or_create(
    filebase = tempfile(),
    dimension = c(n_tps, n_elec),
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

      arr[,idx_e] <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]

      return()
    })
  }

  signalScaling <- 10^floor(log10(max(arr[])))
  arr[] <- arr[]/signalScaling

  ## create adjacency array (array of adj matrices for each time window)
  ## iw: The index of the window we are going to calculate fragility

  res <- raveio::lapply_async(seq_len(n_steps), function(iw) {
    ## Sample indices for the selected window
    si <- seq_len(t_window-1) + (iw-1)*t_step
    ## measurements at time point t
    xt <- arr[si,]
    ## measurements at time point t plus 1
    xtp1 <- arr[si + 1,]

    ## Coefficient matrix A (adjacency matrix)
    ## each column is coefficients from a linear regression
    ## formula: xtp1 = xt*A + E
    Ai <- ridge(xt, xtp1, intercept = F, lambda = lambda, iw = iw)

    R2 <- ridgeR2(xt,xtp1,Ai)

    return(list(Ai = Ai, R2 = R2))
  }, callback = function(iw) {
    sprintf("Running XXX|Step %s", iw)
  })

  A <- unlist(raveio::lapply_async(res, function(w){
    w$Ai
  }))
  dim(A) <- c(n_elec, n_elec, n_steps)
  dimnames(A) <- list(
    Electrode1 = repository$electrode_list,
    Electrode2 = repository$electrode_list,
    Step = seq_len(n_steps)
  )

  R2 <- unlist(raveio::lapply_async(res, function(w){
    w$R2
  }))
  dim(R2) <- c(n_elec, n_steps)
  dimnames(R2) <- list(
    Electrode = repository$electrode_list,
    Step = seq_len(n_steps)
  )

  # calculate fragility
  f <- unlist(raveio::lapply_async(seq_len(n_steps), function(iw){
    fragilityRow(A[,,iw])
  }))
  dim(f) <- c(n_elec, n_steps)
  f_naked=f
  dimnames(f) <- list(
    Electrode = repository$electrode_list,
    Step = seq_len(n_steps)
  )

  ## use ranking to increase contrast
  f_contrast <- matrix(rank(f), nrow(f), ncol(f))
  attributes(f_contrast) <- attributes(f)

  # scale fragility values from -1 to 1 with 1 being most fragile

  # normalize, for each column (margin=2L)
  f_norm <- apply(f_contrast, 2, function(f_col) {
    max_f <- max(f_col)
    min_f <- min(f_col)
    2.0 * (f_col - min_f) / (max_f - min_f) - 1.0 # normalize from -1 to 1
    #(f_col - min_f) / (max_f - min_f) # normalize from 0 to 1
  })

  return(list(
    adj = A,
    frag = f_norm,
    R2 = R2
  ))
}

# calc_error <- function(repository, trial_num, t_window, t_step, index) {
#   n_tps <- length(repository$voltage$dimnames$Time)
#   n_elec <- length(repository$voltage$dimnames$Electrode)
#
#   # Number of steps
#   n_steps <- floor((n_tps - t_window) / t_step) + 1
#
#   # slice of data
#   arr <- filearray::filearray_load_or_create(
#     filebase = tempfile(),
#     dimension = c(n_tps, n_elec),
#     type = "float", mode = "readwrite", partition_size = 1L,
#
#     # if repository has changed, re-calculate
#     repository_signature = repository$signature,
#     t_step = t_step, t_window = t_window,
#     trial_num = trial_num,
#
#     on_missing = function(arr) {
#       arr$set_header("ready", value = FALSE)
#     }
#   )
#
#   # check if header `ready` is not TRUE
#   if(!isTRUE(arr$get_header("ready", FALSE))) {
#
#     loaded_electrodes <- repository$electrode_list
#     raveio::lapply_async(repository$voltage$data_list, function(v) {
#       e <- dimnames(v)$Electrode
#       idx_e <- loaded_electrodes == e
#
#       arr[,idx_e] <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]
#
#       return()
#     })
#   }
#
#   signalScaling <- 10^floor(log10(max(arr[])))
#   arr[] <- arr[]/signalScaling
#
#   library(foreach)
#   library(doParallel)
#   registerDoParallel(parallel::detectCores())
#
#   si <- seq_len(ntw) + (index-1)*ntsl
#
#
#   ## measurements at time point t
#   xt <- timeSeries[si,]
#   ## measurements at time point t plus 1
#   xtp1 <- timeSeries[si + 1,]
#
#   ## Train/validation split
#   train_prop <- 0.8
#   train_size <- floor(nrow(xt)*train_prop)
#   train_ind <- sample(seq_len(nrow(xt)), train_size)
#   train_xt <- xt[train_ind,]
#   train_xtp1 <- xtp1[train_ind,]
#   test_xt <- xt[-train_ind,]
#   test_xtp1 <- xtp1[-train_ind,]
#
#
#
#   #####################
#   ## Verify the ridge regression cv works
#   #####################
#   ## Coefficient matrix A (adjacency matrix)
#   ## each column is coefficients from a linear regression
#   ## formula: xtp1 = xt*A
#   A <- ridge(train_xt, train_xtp1, intercept = F, lambda = 0.0001)
#   Acv <- ridgecv(train_xt, train_xtp1, parallel = TRUE)
#
#   ## two matrix should have the same dimension
#   stopifnot(all.equal(dim(A), dim(Acv)))
#
#   ## prediction
#   pred <- predictRidge(test_xt, A)
#   predcv <- predictRidge(test_xt, Acv)
#
#   ## calculate error
#   err <- pred - test_xtp1
#   errcv <- predcv - test_xtp1
#
#   ## total error squared
#   err2 <- mean(as.matrix(err^2))
#   err2cv <- mean(as.matrix(errcv^2))
#
#   err2
#   err2cv
#
#   ## Check Goodness of Fit
#   R2 <- ridgeR2(xt,xtp1,A)
#   hist(R2)
#
#   R2cv <- ridgeR2(xt,xtp1,Acv)
#   hist(R2cv)
#
#   ## Check heatmap
#   ## set diagonal to 0 to see if there is any outlier
#   ## (diagonal is always outlier)
#   heatA <- A
#   diag(heatA) <- 0
#   heatmap(heatA, Colv = NA, Rowv = NA, scale="none")
#   range(heatA)
#
#
#   heatAcv <- Acv
#   diag(heatAcv) <- 0
#   heatmap(heatAcv, Colv = NA, Rowv = NA, scale="none")
#   range(heatAcv)
# }

threshold_fragility <- function(repository, adj_frag_info, t_start, t_end, threshold = 0.5) {
  n_windows <- dim(adj_frag_info$adj)[3]
  t_step <- floor(length(repository$voltage$dimnames$Time)/n_windows)

  # convert from input t_start and t_end to timewindow indices
  tw_start <- floor(which(t_start==repository$voltage$dimnames$Time)/t_step)
  tw_end <- floor(which(t_end==repository$voltage$dimnames$Time)/t_step)
  if (tw_end > n_windows) { tw_end <- n_windows }

  # subset fragility matrix to specified timewindows
  mat <- adj_frag_info$frag[,tw_start:tw_end]

  avg_f <- rowMeans(mat)
  elec <- which(avg_f > threshold)

  return(list(
    avg_f = avg_f,
    elecnames = attr(elec, "names")
  ))
}

ridgecv <- function(xt, xtp1, parallel=FALSE) {
    if (!identical(dim(xt), dim(xtp1))) {
        stop("Unmatched dimension")
    }
    nel <- ncol(xt)
    x <- as.matrix(xt)
    ## parallel computing backend
    if (parallel){
        library(doParallel)
        library(foreach)
        if(!getDoParRegistered()){
            registerDoParallel(cores=parallel::detectCores())
        }
        A <- foreach(i = seq_len(nel), .packages=c("glmnet"), .export = "ridgecvTwoPass",.combine = cbind) %dopar%{
            y <- xtp1[, i]
            fitcoef <- ridgecvTwoPass(x,y)
            fitcoef <- as.numeric(fitcoef)
            fitcoef[-1]
        }
    }else{
        ## Coefficient matrix A
        ## each column is coefficients from a linear regression
        ## formula: xtp1 = xt*A + E
        A <- matrix(0, nel, nel)
        ## for each electrode
        for (i in seq_len(nel)) {
            y <- xtp1[, i]
            fitcoef <- tryCatch({
                ridgecvTwoPass(x,y)
            }, error = function(e) {
                fit <- glmnet::glmnet(x, y,
                    alpha = 0, lambda = 0,
                    standardize = FALSE,
                    intercept = FALSE
                )
                coef(fit)
            })

            fitcoef <- as.numeric(fitcoef)
            A[, i] <- fitcoef[-1]
        }
        A
    }
}

ridgecvTwoPass <- function(x, y, lambdaRange = 10^-rev(1:10)){
    set.seed(1)
    ## first pass: determine the scale of lambda
    fitcv1 <- glmnet::cv.glmnet(x, y,
        alpha = 0,
        standardize = FALSE,
        intercept = FALSE,
        type.measure="mse",
        lambda = lambdaRange,
        thresh = 1e-10
    )
    # plot(fitcv1)
    bestlambda <- fitcv1$lambda.min

    ## if the best lambda is the smallest one, then we do not use ridge regression
    if (bestlambda==min(lambdaRange)) {
        fit <- glmnet::glmnet(x, y,
            alpha = 0, lambda = 0,
            standardize = FALSE,
            intercept = FALSE
        )
        return(coef(fit))
    }

    ## second pass: determine the optimal lambda
    lambdaRange2 <- seq(bestlambda/10, bestlambda*10, length.out = 1000)
    fitcv2 <- glmnet::cv.glmnet(x, y,
        alpha = 0,
        standardize = FALSE,
        intercept = FALSE,
        type.measure="mse",
        lambda = lambdaRange2,
        thresh = 1e-10
    )

    coef(fitcv2, s = fitcv2$lambda.min)
}

predictRidge <- function(xt, A) {
    ## the data matrix
    if (nrow(A) == ncol(A) + 1) {
        x <- cbind(1, as.matrix(xt))
    } else {
        x <- as.matrix(xt)
    }
    x %*% A
}


##### Deprecated functions

# generate_adjacency_array <- function(repository, trial_num, t_window, t_step, nlambda, signalScaling) {
#   n_tps <- length(repository$voltage$dimnames$Time)
#   n_elec <- length(repository$voltage$dimnames$Electrode)
#
#   # Number of steps
#   n_steps <- floor((n_tps - t_window) / t_step) + 1
#
#   # slice of data
#   arr <- filearray::filearray_load_or_create(
#     filebase = tempfile(),
#     dimension = c(t_window, n_steps, n_elec),
#     type = "float", mode = "readwrite", partition_size = 1L,
#
#     # if repository has changed, re-calculate
#     repository_signature = repository$signature,
#     t_step = t_step, t_window = t_window,
#     trial_num = trial_num,
#
#     on_missing = function(arr) {
#       arr$set_header("ready", value = FALSE)
#     }
#   )
#
#   # check if header `ready` is not TRUE
#   if(!isTRUE(arr$get_header("ready", FALSE))) {
#
#     loaded_electrodes <- repository$electrode_list
#     raveio::lapply_async(repository$voltage$data_list, function(v) {
#       e <- dimnames(v)$Electrode
#       idx_e <- loaded_electrodes == e
#
#       trial_voltage <- v[, trial_num, 1, drop = TRUE, dimnames = NULL]/signalScaling
#
#       idx <- seq_len(t_window)
#       lapply(seq_len(n_steps), function(step) {
#         t_start <- 1 + (step - 1) * t_step
#         #t_start <- (step - 1) * t_step
#         arr[, step, idx_e] <- trial_voltage[t_start + idx]
#       })
#       return()
#     })
#
#   }
#
#   # calculate adjacency arrays
#   A <- raveio::lapply_async(
#     seq_len(n_steps), function(step) {
#       slice <- arr[, step, , drop = FALSE, dimnames = NULL]
#       dm <- dim(slice)
#       nr <- nrow(slice)
#       dim(slice) <- c(nr, dm[[3]])
#       # x is x(t) and y is x(t+1), state vectors
#       x <- slice[-nr, , drop = FALSE]
#       y <- slice[-1, , drop = FALSE]
#       Ai <- ridge(xt,xtp1, intercept = FALSE, lambda)
#       # AEigen <- eigen(Ai)
#       # e <- Mod(AEigen$values)
#       # me <- max(e)
#       # if(me <1){
#       #   message("Stable solution")
#       # }else{
#       #   message("Solution is not stable")
#       #   message(me)
#       # }
#       as.vector(Ai)
#       # as.vector(find_adj_matrix(x = x, y = y, nlambda = nlambda))
#       # as.vector(find_adj_matrix_bij(x = x, y = y))
#     }
#   )
#   A <- do.call("cbind", A)
#   dim(A) <- c(n_elec, n_elec, n_steps)
#   A
# }
#
# find_adj_matrix <- function(x, y, nlambda) {
#
#   # x <- t(state_vectors$x)
#   # y <- t(state_vectors$x_n)
#
#   nobs <- nrow(x)
#   nvars <- ncol(x)
#
#   # scale x and y
#   x <- scale(x)
#   y <- scale(y)
#   scale_x <- attr(x, "scaled:scale")
#   scale_y <- attr(y, "scaled:scale")
#
#   # guess possible lambdas
#   lambdas <- rev(lambda_path(
#     x = x, y = y, nlambda = nlambda, alpha = 0.0,
#     standardize = TRUE, intercept = FALSE))
#
#   XtY <- crossprod(x, y)
#   XtX <- crossprod(x)
#
#   # SVD decomposition to speed things up
#   # XtX == svd$u %*% diag(svd$d) %*% t(svd$v)
#   svd <- svd(XtX)
#
#   # solve(XtX) == svd$v %*% diag(1 / svd$d) %*% t(svd$u)
#   V <- svd$v
#   Ut <- t(svd$u)
#
#   # sanity check: the following should be identity matrix
#   # (XtX + diag(12, nvars)) %*% (V %*% diag(1 / (svd$d + 12)) %*% Ut)
#   UtXtY <- Ut %*% XtY
#
#   # for each, lambda, fit ridge regression
#   dipsaus::forelse(
#     x = lambdas,
#     FUN = function(lam) {
#
#       # The following 3 methods generate similar/same results
#       # MASS::lm.ridge(y[,2] ~ x - 1, lambda = lam * nobs)
#       # glmnet::glmnet(y = y[,1], x = x, intercept = FALSE, lambda = lam, alpha = 0.0)$beta
#
#       # no SVD method
#       # ident <- diag(as.double(nobs), ncol(x))
#       # adj_matrix <- solve(XtX + lam * ident) %*% XtY
#
#       # SVD method
#       adj_matrix <- V %*% diag(1 / (svd$d + lam * nobs)) %*% UtXtY
#
#       # right now scale(x) %*% adj_matrix = scale(y), need to scale back
#       adj_matrix <- diag( 1 / scale_x ) %*% adj_matrix %*% diag( scale_y )
#
#
#       eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
#       stable <- max(eigv) < 1
#
#       if( !stable ) { return() }
#       # return(list(
#       #   adj = adj_matrix,
#       #   abs_eigv = eigv,
#       #   stable = stable
#       # ))
#
#       structure(adj_matrix)
#       return(adj_matrix)
#
#     },
#     ALT_FUN = function() {
#       stop('No lambdas result in a stable adjacency matrix. Increase the number of lambdas, or (more likely) there is something wrong with your data.')
#     }
#   )
#
# }
#
# generate_fragility_matrix <- function(A, elec, lim = 1i, ncores) {
#   print('Generating fragility matrix')
#
#   dm <- dim(A)
#   N <- dm[1]
#   J <- dm[3]
#
#   f_vals <- raveio::lapply_async(
#     seq_len(J),
#     function(k) {
#       A_k <- A[,, k]
#       f_vals_k <- vapply(seq_len(N), function(i){
#         find_fragility(i, A_k = A_k, N = N, limit = lim)
#       }, FUN.VALUE = 0.0)
#       f_vals_k
#       #fragilityRow(A[,,k])
#     }
#   )
#   # N x J Fragility Matrix
#   f_vals <- do.call("cbind", f_vals)
#
#   dimnames(f_vals) <- list(
#     Electrode = elec,
#     Step = seq_len(J)
#   )
#
#   # scale fragility values from -1 to 1 with 1 being most fragile
#
#   # normalize, for each column (margin=2L)
#   f_norm <- apply(f_vals, 2, function(f_col) {
#     max_f <- max(f_col)
#     min_f <- min(f_col)
#     2.0 * (f_col - min_f) / (max_f - min_f) - 1.0 # normalize from -1 to 1
#     #(f_col - min_f) / (max_f - min_f) # normalize from 0 to 1
#   })
#
#   # find average fragility for each electrode across time
#   avg_f <- rowMeans(f_norm)
#
#   return(list(
#     vals = f_vals,
#     norm = f_norm,
#     avg = avg_f
#   ))
# }
#
# find_fragility <- function(node, A_k, N, limit) {
#
#   e_k <- vector(mode = 'numeric', length = N)
#   e_k[node] <- 1
#
#   argument <- t(e_k) %*% (solve(A_k - limit*diag(N))) # column perturbation
#   # argument <- t(e_k) %*% t(solve(A_k - num*diag(N))) # row perturbation
#
#   B <- rbind(Im(argument),Re(argument))
#
#   perturb_mat <- (t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) %*% t(e_k) # column
#   # perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row
#
#   norm(perturb_mat, type = '2')
# }
