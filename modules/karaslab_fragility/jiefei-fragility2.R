#' Calculate the fragility of an electrode
#'
#' @param A The adjacency matrix
#' @param nSearch The search number for the grid search of the optimal complex number
#'
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




ridge <- function(xt, xtp1, lambda = 0.0001, intercept = FALSE) {
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
        fit <- glmnet(xt, y,
            alpha = 0, lambda = lambda,
            standardize = FALSE, intercept = intercept
        )

        if (intercept) {
            A[, i] <- as.numeric(coef(fit))
        } else {
            A[, i] <- coef(fit)[-1]
        }
    }
    A
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
                fit <- glmnet(x, y,
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
    fitcv1 <- cv.glmnet(x, y,
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
        fit <- glmnet(x, y,
            alpha = 0, lambda = 0,
            standardize = FALSE, 
            intercept = FALSE
        )
        return(coef(fit))
    }

    ## second pass: determine the optimal lambda
    lambdaRange2 <- seq(bestlambda/10, bestlambda*10, length.out = 1000)
    fitcv2 <- cv.glmnet(x, y,
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

ridgeR2 <- function(xt, xtp1, A) {
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
