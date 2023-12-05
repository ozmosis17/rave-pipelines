#'Calculate the fragility of an electrode
#'
#'@param A The adjacency matrix
#'@param nSearch The search number for the grid search of the optimal complex number
#'
fragilityRow <- function(A, nSearch = 100) {
    ## The adjacency matrix A here is a transpose of the
    ## adjacency matrix in the original paper
    nel <- ncol(A)
    e <- Mod(eigen(A)$values)
    me <- max(e)

    if (me >= 1) {
        # return(0)
    }


    fragcol <- matrix(0, nel,nel)
    fragNorm <- rep(0,nel)
    omvec <- seq(0, 1, length.out = nSearch+1)[-1]

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
            if(FALSE){
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

        fragcol[,i] <- minPerturbColumn
        fragNorm[i] <- minNorm
    }

    maxf <- max(fragNorm)
    fragNorm2 <- (maxf - fragNorm) / maxf

    return(fragNorm2)
}




ridge <- function(xt, xtp1, lambda = 0.0001, intercept = FALSE){
    if(!identical(dim(xt),dim(xtp1)))
        stop("Unmatched dimension")
    nel <- ncol(xt)
    ## Coefficient matrix A
    ## each column is coefficients from a linear regression
    ## formula: xtp1 = xt*A + E
    A <- matrix(0, nel + intercept,nel)
    ## for each electrode
    for(i in seq_len(nel)){
        y=xtp1[,i]
        fit <- glmnet(xt, y, alpha = 0, lambda  = lambda,
                      standardize =FALSE,intercept =intercept)
        if(intercept)
            A[,i] <- as.numeric(coef(fit))
        else
            A[,i] <- coef(fit)[-1]
    }
    A
}


ridgeR2 <- function(xt, xtp1, A){
    if(!identical(dim(xt),dim(xtp1)))
        stop("Unmatched dimension")
    nel <- ncol(xt)

    ## the data matrix
    if(nrow(A)==ncol(A)+1){
        x <- cbind(1, as.matrix(xt))
    }else{
        x <- as.matrix(xt)
    }

    ypredMat <- x %*% A

    R2 <- rep(0, nel)
    for(i in seq_len(nel)) {
        y <- xtp1[,i]
        ypred <- ypredMat[,i]
        sst <- sum((y - mean(y))^2)
        sse <- sum((ypred - y)^2)
        rsq <- 1 - sse/sst
        R2[i] <- rsq
    }
    R2
}





