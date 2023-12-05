library(glmnet)
source("fragility.R")

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


data("iris")


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


## Coefficient matrix A (adjacency matrix)
## each column is coefficients from a linear regression
## formula: xtp1 = xt*A
A <- ridge(xt,xtp1, intercept = F, lambda = 0.0001)
R2 <- ridgeR2(xt,xtp1,A)

## Check heatmap
## set diagonal to 0 to see if there is any outlier
## (diagonal is always outlier)
heatA <- A
diag(heatA) <- 0
heatmap(heatA, Colv = NA, Rowv = NA, scale="none")

## Check Goodness of Fit
hist(R2)

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
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores())

## The number of samples in each window
ntw <- 250
## sliding time increment
ntsl <- 125
## total number of possible windows
nWins <- floor(nrow(timeSeries)/ntsl-1)
## iw: The index of the window we are going to calculate fragility
res <- foreach(iw = seq_len(nWins),
               .combine = cbind,
               .packages = "glmnet")%dopar%{
    ## Sample indices for the selected window
    si <- seq_len(ntw) + (iw-1)*ntsl
    ## measurements at time point t
    xt <- timeSeries[si,]
    ## measurements at time point t plus 1
    xtp1 <- timeSeries[si + 1,]

    ## Coefficient matrix A (adjacency matrix)
    ## each column is coefficients from a linear regression
    ## formula: xtp1 = xt*A + E
    A <- ridge(xt,xtp1, intercept = F)

    ## fragility vector
    fragilityRow(A)
}
rownames(res) <- channelDefs$name
colnames(res) <- as.character(round((seq_len(nWins)-1)*ntsl/1000))

## use ranking to increase contrast
res2 <- matrix(rank(res), nrow(res), ncol(res))
attributes(res2) <- attributes(res)
library(RColorBrewer)
coul <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(256)
heatmap(res2, Colv = NA, Rowv = NA,
        scale="none", col = coul, xlab = "time")
