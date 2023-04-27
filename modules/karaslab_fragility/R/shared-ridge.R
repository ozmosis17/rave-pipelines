lambda_path <- function(
    x, y, weights = NULL, nlambda = 100,
    lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04),
    alpha = 0.0, offset = NULL, family = stats::gaussian(),
    standardize = TRUE, intercept = TRUE) {

  # DIPSAUS DEBUG START
  # a = readRDS("~/Desktop/junk")
  # x = a$x
  # y = a$y
  # nlambda = 16
  # np = dim(x)
  # if (is.null(np) | (np[2] <= 1)) {
  #   stop("x should be a matrix with 2 or more columns")
  # }
  # nobs = as.integer(np[1])
  # nvars = as.integer(np[2])
  # list2env(list(weights = NULL,
  #               lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04),
  #               alpha = 0, offset = NULL, family = stats::gaussian(),
  #               standardize = TRUE, intercept = TRUE), envir=.GlobalEnv)


  np = dim(x)
  if (is.null(np) | (np[2] <= 1)) {
    stop("x should be a matrix with 2 or more columns")
  }
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  if (any(is.na(x)))  {
    stop("x has missing values; consider using makeX() to impute them")
  }
  if (is.null(weights)) {
    weights = rep(1, nobs)
  } else if (length(weights) != nobs) {
    stop(paste("number of elements in weights (", length(weights),
               ") not equal to the number of rows of x (", nobs,
               ")", sep = ""))
  }

  if( lambda.min.ratio >= 1 ) {
    stop("lambda.min.ratio must be less than 1")
  } else if ( lambda.min.ratio <= 0 ) {
    lambda.min.ratio <- ifelse(nobs < nvars, 0.01, 1e-04)
  }

  nlambda <- round(nlambda)
  if( nlambda < 1 ) {
    stop("nlambda must be positive integer")
  }

  weights <- as.double(weights)
  etastart = 0
  mustart = NULL
  start = NULL
  eval(family$initialize)
  y <- drop(y)
  is.offset <- !(is.null(offset))
  if (is.offset == FALSE) {
    offset <- as.double(y * 0)
  }

  # calculate weighted mean and sd for x
  weights <- weights/sum(weights)
  xm <- drop(t(weights) %*% x)
  xv <- drop(t(weights) %*% scale(x, xm, FALSE)^2)
  xv[xv < 10 * .Machine$double.eps] <- 0
  xs <- sqrt(xv)

  if (!intercept) {
    xm <- rep(0, times = nvars)
  }
  if (!standardize) {
    xs <- rep(1, times = nvars)
  }
  x <- scale(x, xm, xs)

  if (intercept) {
    if (is.offset) {
      suppressWarnings(tempfit <- glm(y ~ 1, family = family,
                                      weights = weights, offset = offset))
      mu <- tempfit$fitted.values
    }
    else {
      mu <- drop(t(weights) %*% y)
      # dim(mu) <- c(nobs, length(mu) / nobs)
      mu <- rep(mean(mu), times = nobs)
    }
  } else {
    mu <- family$linkinv(offset)
  }
  ju <- rep(1, nvars)
  r <- y - mu
  eta <- family$linkfun(mu)
  v <- family$variance(mu)
  m.e <- family$mu.eta(eta)
  weights <- weights/sum(weights)
  rv <- r/v * m.e * weights
  g <- abs(drop(t(rv) %*% x))
  g <- g * ju
  lambda_max <- max(g)/max(alpha, 0.001)
  lambda_max * lambda.min.ratio
  exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                  length.out = nlambda))


}
