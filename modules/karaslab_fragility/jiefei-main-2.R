library(glmnet)
source("fragility2.R")

# library(doRedis)
# registerDoRedis("RJOBS")
# getDoParWorkers()
# setProgress(value = TRUE)
# doRedis::removeQueue("RJOBS")

## parallel computing backend
library(foreach)
library(doParallel)
registerDoParallel(cores = parallel::detectCores())


## Determine the signal strength (larger scaling -> weaker signal)
signalScaling <- 1000000

library(readxl)
channelDefs <- read_xls("Pt01ictalRun01EcoGChannels.xls")
goodChannels <- c(1:4,7:36,42:43,46:69,72:95)
channelDefs <- channelDefs[goodChannels,]


## Read signal data, each column will be one electrode
library(R.matlab)
x <- readMat("timeSeriesSubpt01Run1.mat")
timeSeries <- data.frame(x$timeSeries/signalScaling)
colnames(timeSeries) <- channelDefs$name
## Number of electrodes
nel <- ncol(timeSeries)
## Number of time points
nt <- nrow(timeSeries)



#####################
# Signal visualization
#####################
## select channels and standardize signal
gaps <- 2
displayChannels <- c(31,32,53:56,57:60, 77:80)
displayNames <- channelDefs$name[displayChannels]
displayNames
plotData <- timeSeries[,displayChannels]
for(i in seq_along(plotData)){
    plotData[, i] <- plotData[, i]- mean(plotData[, i]) +
        (ncol(plotData)-i)*gaps
}

plot(plotData[, 1],type="l" ,cex=0.1,
     ylim = range(plotData), yaxt = "n")
for(i in 2:ncol(plotData)){
    lines(plotData[, i])
}
axis(2, at = rev(seq_along(displayChannels) - 1)*gaps,
     labels = displayNames,las=1)


#####################
# Compute fragility at a time window
#####################
## The number of samples in each window
ntw <- 250
## sliding time increment
ntsl <- 125
## The index of the window we are going to calculate fragility
iw <- 80
## Sample indices for the selected window
si <- seq_len(ntw) + (iw-1)*ntsl


## measurements at time point t
xt <- timeSeries[si,]
## measurements at time point t plus 1
xtp1 <- timeSeries[si + 1,]

## Train/validation split
train_prop <- 0.8
train_size <- floor(nrow(xt)*train_prop)
train_ind <- sample(seq_len(nrow(xt)), train_size)
train_xt <- xt[train_ind,]
train_xtp1 <- xtp1[train_ind,]
test_xt <- xt[-train_ind,]
test_xtp1 <- xtp1[-train_ind,]



#####################
## Verify the ridge regression cv works
#####################
## Coefficient matrix A (adjacency matrix)
## each column is coefficients from a linear regression
## formula: xtp1 = xt*A
A <- ridge(train_xt, train_xtp1, intercept = F, lambda = 0.0001)
Acv <- ridgecv(train_xt, train_xtp1, parallel = TRUE)

## two matrix should have the same dimension
stopifnot(all.equal(dim(A), dim(Acv)))

## prediction
pred <- predictRidge(test_xt, A)
predcv <- predictRidge(test_xt, Acv)

## calculate error
err <- pred - test_xtp1
errcv <- predcv - test_xtp1

## total error squared
err2 <- mean(as.matrix(err^2))
err2cv <- mean(as.matrix(errcv^2))

err2
err2cv

## Check Goodness of Fit
R2 <- ridgeR2(xt,xtp1,A)
hist(R2)

R2cv <- ridgeR2(xt,xtp1,Acv)
hist(R2cv)

## Check heatmap
## set diagonal to 0 to see if there is any outlier
## (diagonal is always outlier)
heatA <- A
diag(heatA) <- 0
heatmap(heatA, Colv = NA, Rowv = NA, scale="none")
range(heatA)


heatAcv <- Acv
diag(heatAcv) <- 0
heatmap(heatAcv, Colv = NA, Rowv = NA, scale="none")
range(heatAcv)


#####################
## fit the ridge regression with all data
#####################
A <- ridgecv(xt, xtp1, parallel = TRUE)


## checking adjacency matrix stability
AEigen <- eigen(A)
e <- Mod(AEigen$values)
me <- max(e)
if(me <1){
    message("Stable solution")
}else{
    message("Solution is not stable")
}


## fragility for each electrode
fragrow <- fragilityRow(A)

fragrow



#####################
# Compute Fragility matrix
#####################
## parallel computing
## The number of samples in each window
ntw <- 250
## sliding time increment
ntsl <- 125
## total number of possible windows
nWins <- floor(nrow(timeSeries)/ntsl-1)
## iw: The index of the window we are going to calculate fragility
start_time <- Sys.time()
res <- foreach(iw = seq_len(nWins),
               .packages = "glmnet",
               .errorhandling = "pass")%dopar%{
    ## Sample indices for the selected window
    si <- seq_len(ntw) + (iw-1)*ntsl
    ## measurements at time point t
    xt <- timeSeries[si,]
    ## measurements at time point t plus 1
    xtp1 <- timeSeries[si + 1,]

    ## Coefficient matrix A (adjacency matrix)
    ## each column is coefficients from a linear regression
    ## formula: xtp1 = xt*A + E
    A <- ridgecv(xt,xtp1)
    # A <- ridge(xt,xtp1)

    ## fragility vector
    fragilityRow(A)
}
end_time <- Sys.time()
end_time - start_time

# saveRDS(res, "fragility_matrix.rds")

## combine list by column
res2 <- do.call(cbind, res)
rownames(res2) <- channelDefs$name
colnames(res2) <- as.character(round((seq_len(nWins)-1)*ntsl/1000))

## use ranking to increase contrast
res3 <- matrix(rank(res2), nrow(res2), ncol(res2))
attributes(res3) <- attributes(res2)
library(RColorBrewer)
coul <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(256)
heatmap(res3, Colv = NA, Rowv = NA,
        scale="none", col = coul, xlab = "time")
