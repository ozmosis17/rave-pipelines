#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
# library to read Matlab data formats into R
library(R.matlab)
# R library to be able to import python libraries
library(reticulate)


# importing python libraries
use_python("/Users/ozhou/anaconda3/bin/python3")
np <- import("numpy")
plt <- import("matplotlib.pyplot")

# script for testing one time window-----
repository <- vplot_data$repository
trial_num <- vplot_data$trial_num
t_window <- vplot_data$t_window
t_step <- vplot_data$t_step
nlambda <- 16
step <- 1

# from generate_adjacency_array function
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

slice <- arr[, step, , drop = FALSE, dimnames = NULL]
dm <- dim(slice)
nr <- nrow(slice)
dim(slice) <- c(nr, dm[[3]])
# x is x(t) and y is x(t+1), state vectors
x <- slice[-nr, , drop = FALSE]
y <- slice[-1, , drop = FALSE]

h <-matrix(0,n_elec*(t_window-1),n_elec*n_elec)
b   <- vector(mode="numeric", length=(t_window-1)*n_elec)

# from ACL find_adjacency_matrix function
# from Li 2017
# n_elec - number of electrodes
# t_window - size of timewindow
# xt - state vector 1
# xtp - state vector 2
# h -
# b -
# return sol  adjacency matrix as a column

nch <- n_elec*n_elec
nr <-  n_elec*(t_window-1)

h[1:nr,1:nch]<-0.0

nb  <- (t_window-1)*n_elec

b[1:nb]<-0.0

# set of undetermined linear equations to be solved Hx=b

for(ie in 1:n_elec){
  for(it in 1:(t_window-1)){
    b[ie+n_elec*(it-1)]<-y[it,ie];
  }
}

for(it in 1:(t_window-1)){
  for(je in 1:n_elec){
    for(ie in 1:n_elec){
      h[je+(it-1)*n_elec,ie+(je-1)*n_elec]=x[it,ie];
    }
  }
}

# get H pseudo inverse

start <- Sys.time()
hinv <- np$linalg$pinv(h)
print(Sys.time()-start)

sol <- hinv %*% b

A1 <- as.vector(sol)
dim(A1) <- c(n_elec, n_elec)

eigv <- abs(eigen(A1, only.values = TRUE)$values)
print(paste0('largest eigenvalue norm: ', max(eigv)))

# graph voltage traces for comparison
y1 <- x
y2 <- t(vplot_data$A[,,1] %*% t(x))
y2 <- t(A1 %*% t(x))
df <- data.frame(seq_len(t_window-1),y1,y2)

# calculate mean squared error between y1 and y2
mse <- mean((y1 - y2)^2)
print(paste0('MSE: ', mse))

g <- ggplot(df, aes(seq_len(t_window-1))) +
  geom_line(aes(y=y1[seq_len(t_window-1),elec_num], color = "original")) +
  geom_line(aes(y=y2[seq_len(t_window-1),elec_num], color = "reconstructed")) +
  labs(x = "Time (ms)", y = paste0("Voltage - Electrode ", elec_num), color = "Legend") +
  scale_color_manual(values = c("original" = 'black', "reconstructed" = "red")) +
  ggtitle(paste0("Mean Squared Error: ", format(mse, scientific = TRUE)))

g

#################
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

  start <- Sys.time()
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
      #as.vector(find_adj_matrix(x = x, y = y, nlambda = nlambda))


      h <-matrix(0,n_elec*(t_window-1),n_elec*n_elec)
      b   <- vector(mode="numeric", length=(t_window-1)*n_elec)
      as.vector(find_adjacency_matrix(n_elec = n_elec, t_window = t_window, x = x, y = y, h = h, b = b))
    }
  )
  A <- do.call("cbind", A)
  dim(A) <- c(n_elec, n_elec, n_steps)
  A
  print(Sys.time()-start)
}

# from Li 2017
# n_elec - number of electrodes
# t_window - size of timewindow
# xt - state vector 1
# xtp - state vector 2
# h -
# b -
# return sol  adjacency matrix as a column

find_adjacency_matrix<-function(n_elec, t_window, x, y, h, b){

  nch <- n_elec*n_elec
  nr <-  n_elec*(t_window-1)

  h[1:nr,1:nch]<-0.0

  nb  <- (t_window-1)*n_elec

  b[1:nb]<-0.0

  # set of undetermined linear equations to be solved Hx=b

  for(ie in 1:n_elec){
    for(it in 1:(t_window-1)){
      b_old[ie+n_elec*(it-1)]<-y[it,ie];
    }
  }

  for(it in 1:(t_window-1)){
    for(je in 1:n_elec){
      for(ie in 1:n_elec){
        h[je+(it-1)*n_elec,ie+(je-1)*n_elec]=x[it,ie];
      }
    }
  }

  # get H pseudo inverse

  start <- Sys.time()
  hinv <- np$linalg$pinv(h)
  print(Sys.time()-start)

  sol <- hinv %*% b

  return(sol)

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

#################

n_elec <- 83
nom <- 100
adjm <- vplot_data$A[,,1]

testrow <- generate_fragilityrowwindow(n_elec=83,nom=100,adjm = vplot_data$A[,,1])
fmatrix <- matrix(nrow=83,ncol=159)
for(window in 1:159) {
  print(window)
  start <- Sys.time()
  fmatrix[,window] <- generate_fragilitycolwindow(n_elec=83,nom=100,adjm = vplot_data$A[,,window])
  print(Sys.time()-start)
}
# normalize fmatrix
f_normalized <- apply(f_matrix, 2, function(f_col) {
  max_f <- max(f_col)
  min_f <- min(f_col)
  2.0 * (f_col - min_f) / (max_f - min_f) - 1.0 # normalize from -1 to 1
  #(f_col - min_f) / (max_f - min_f) # normalize from 0 to 1
})

testcol <- generate_fragilitycolwindow(n_elec=83,nom=100,adjm = vplot_data$A[,,1])
# n_elec - number of electrodes
# nom - number of discretized omegas (100)
# adjm - adjacency matrix
generate_fragilityrowwindow<-function(n_elec, nom, adjm){

  # compute row perturbation
  e_k <- vector(mode = 'numeric', length = n_elec)
  frag_row <- vector(mode = 'numeric', length = n_elec)

  omvec <- seq(0.0, 1.0, length.out=nom)

  minindnorm <-100000
  iomm <-1

  ie<-1
  e_k[ie]<-1

  for(iom in 2:nom){

    #iom<-2

    om<- omvec[iom]

    sigma = sqrt(1-om*om)

    lambda <- sigma+ 1i*om

    argument <- t(e_k) %*% t(solve(adjm - lambda*diag(n_elec))) # row perturbation
    B <- rbind(Im(argument),Re(argument))

    perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row

    # compute 2-induced matrix norm of perturbation matrix delta : sqrt(lamba_max (t(delta)delta))
    tpertpert<-t(perturb_mat)%*%perturb_mat
    ev <- eigen(tpertpert)
    lammax <- ev$values[1]

    twoindnorm<-sqrt(lammax)
    #indnormv[iom]<-twoindnorm


    if(twoindnorm<minindnorm){
      minindnorm<-twoindnorm
      iomm<-iom
      #print(iomm)
    }


  }
  #print(iomm)
  #print(minindnorm)
  frag_row[ie]<-minindnorm

  for(ie in 2:n_elec){

    e_k[ie-1]<-0
    e_k[ie]<-1

    minindnorm <-100000
    iomm <-1


    for(iom in 2:nom){

      #iom<-2

      om<- omvec[iom]

      sigma = sqrt(1-om*om)

      lambda <- sigma+ 1i*om

      argument <- t(e_k) %*% t(solve(adjm - lambda*diag(n_elec))) # row perturbation
      B <- rbind(Im(argument),Re(argument))

      perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row

      # compute 2-induced matrix norm of perturbation matrix delta : sqrt(lamba_max (t(delta)delta))
      tpertpert<-t(perturb_mat)%*%perturb_mat
      ev <- eigen(tpertpert)
      lammax <- ev$values[1]

      twoindnorm<-sqrt(lammax)
      #indnormv[iom]<-twoindnorm


      if(twoindnorm<minindnorm){
        minindnorm<-twoindnorm
        iomm<-iom
        #print(iomm)
      }


    }
    #print(iomm)
    #print(minindnorm)
    frag_row[ie]<-minindnorm

  }

  # maxf<-max(frag_row)
  #
  # frag_row[1:n_elec]<-(maxf-frag_row[1:n_elec])/maxf

  return(frag_row)

}


generate_fragilitycolwindow<-function(n_elec, nom, adjm){

  # compute row perturbation
  e_k <- vector(mode = 'numeric', length = n_elec)
  frag_col <- vector(mode = 'numeric', length = n_elec)

  omvec <- seq(0.0, 1.0, length.out=nom)

  minindnorm <-100000
  iomm <-1

  ie<-1
  e_k[ie]<-1

  for(iom in 2:nom){

    #iom<-2

    om<- omvec[iom]

    sigma = sqrt(1-om*om)

    lambda <- sigma+ 1i*om

    #argument <- t(e_k) %*% t(solve(adjm - lambda*diag(n_elec))) # row perturbation
    argument <- t(e_k) %*% (solve(adjm - lambda*diag(n_elec))) # column perturbation
    B <- rbind(Im(argument),Re(argument))

    #perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row
    perturb_mat <- (t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) %*% t(e_k) # column

    # compute 2-induced matrix norm of perturbation matrix delta : sqrt(lamba_max (t(delta)delta))
    tpertpert<-t(perturb_mat)%*%perturb_mat
    ev <- eigen(tpertpert)
    lammax <- ev$values[1]

    twoindnorm<-sqrt(lammax)
    #indnormv[iom]<-twoindnorm


    if(twoindnorm<minindnorm){
      minindnorm<-twoindnorm
      iomm<-iom
      #print(iomm)
    }


  }
  #print(iomm)
  #print(minindnorm)
  frag_col[ie]<-minindnorm

  for(ie in 2:n_elec){

    e_k[ie-1]<-0
    e_k[ie]<-1

    minindnorm <-100000
    iomm <-1


    for(iom in 2:nom){

      #iom<-2

      om<- omvec[iom]

      sigma = sqrt(1-om*om)

      lambda <- sigma+ 1i*om


      #argument <- t(e_k) %*% t(solve(adjm - lambda*diag(n_elec))) # row perturbation
      argument <- t(e_k) %*% (solve(adjm - lambda*diag(n_elec))) # column perturbation
      B <- rbind(Im(argument),Re(argument))

      #perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row
      perturb_mat <- (t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) %*% t(e_k) # column

      # compute 2-induced matrix norm of perturbation matrix delta : sqrt(lamba_max (t(delta)delta))
      tpertpert<-t(perturb_mat)%*%perturb_mat
      ev <- eigen(tpertpert)
      lammax <- ev$values[1]

      twoindnorm<-sqrt(lammax)
      #indnormv[iom]<-twoindnorm


      if(twoindnorm<minindnorm){
        minindnorm<-twoindnorm
        iomm<-iom
        #print(iomm)
      }


    }
    #print(iomm)
    #print(minindnorm)
    frag_col[ie]<-minindnorm

  }

  # maxf<-max(frag_row)
  #
  # frag_row[1:n_elec]<-(maxf-frag_row[1:n_elec])/maxf

  return(frag_col)

}

leaddata<-read.csv('/Users/aclesage/Documents/RAVEProjects/fragility/fragilityDiamondData/leadloc.csv')
leadpos<-as.matrix(leaddata)
# scatterplot3d(leadpos, pch = 16, color="steelblue",
#               main="Lead Plot",
#               xlab = "x",
#               ylab = "y",
#               zlab = "z")

# estimated from epileptogenic index
onset=12500
onsetlead<-c(10,17,31,41,45)

# read in our data
inputtimeseries <- readMat("/Users/aclesage/Documents/RAVEProjects/fragility/fragilityDiamondData/timeseries.mat")
#inputtimeseries[1,1]
mode(inputtimeseries)
length(inputtimeseries)

# inputtimeseries is read like a list
names(inputtimeseries) <- c("data")

# inputtimeseries is read like a list
nt<-144000; # number of time steps
n_elec<-64; # number of electrodes
ns<-1
ne<-144000
t<-1:144000;
f=1000;  # acquisition frequency
dt<-1/f; # time step
t<-t*dt;

#############################
#
#.   Epoching
###########################

nepoch<-4000

xepoch1 <-matrix(0,nepoch,n_elec)

for(ie in 1:n_elec){
  ns<-(ie-1)*nt+10000
  ne<-ns+nepoch-1
  xepoch1[1:nepoch,ie]<-inputtimeseries$"data"[ns:ne]
}

xepoch <- scale(xepoch1)

##############################################################################################################
# Fragility analysis

t_window=250 # window size for fragility analysis
ntsl=125 # sliding step
nw=30 # window number

xt  <- matrix(0,t_window,n_elec)
xtp <- matrix(0,t_window,n_elec)

nch <- n_elec*n_elec
nr <-  n_elec*(t_window-1)


nb  <- t_window*n_elec
b   <- vector(mode="numeric", length=nb)


nch <- n_elec*n_elec
nr <-  n_elec*t_window

h <-matrix(0,nr,nch)

nb  <- t_window*n_elec
b   <- vector(mode="numeric", length=nb)

msev   <- vector(mode="numeric", length=nw)

fragility <-matrix(0,n_elec,nw)

nom<-101

for(iw in 1:nw){
  #for(iw in 1:29){
  #iw=1

  print(iw)

  nts=1+(iw-1)*ntsl
  nte=nts+t_window-1

  ntsp=nts+1
  ntep=nte+1

  for(ie in 1:n_elec){
    xt[1:t_window,ie]=xepoch[nts:nte,ie]
    xtp[1:t_window,ie]=xepoch[ntsp:ntep,ie]
  }

  start <- Sys.time()
  print(start)
  sol<-find_adjacency_matrix(n_elec,t_window,xt,xtp,h,b)
  end <- Sys.time()
  print(end - start)

  # build adjacency matrix/ transition matrix (rearrange linear equation solution)
  adjm <-matrix(sol,n_elec,n_elec,byrow=TRUE)


  checkval=FALSE



  txt<-t(xt)
  txtp<-t(xtp)
  # validate data reconstruction with transition matrix

  xtpv = adjm %*% txt


  # mean square error
  sse=sum((txtp-xtpv)^2)
  mse = sse/(ne*ne)
  print(mse)

  msev[iw]<-mse

  if(checkval==TRUE){
    # basic R plotting with layout
    # see help(plot)
    layout(matrix(1:4,2,2))
    layout.show(4)


    plot(txtp[1,1:t_window],type='l')
    plot(xtpv[1,1:t_window],type='l')


    plot(txtp[10,1:t_window],type='l')
    plot(xtpv[10,1:t_window],type='l')

  }


  frag_row<-generate_fragilityrowwindow(n_elec,nom,adjm)
  fragility[1:n_elec,iw]<-frag_row[1:n_elec]

  # frag_col<-generate_fragilitycolwindow(n_elec,nom,adjm)
  # fragility[1:n_elec,iw]<-frag_col[1:n_elec]

}

fragsave <- fragility

maxf<-max(fragility)
fragility[1:n_elec,1:nw]<-(maxf-fragility[1:n_elec,1:nw])/maxf

fragilmap <-matrix(0,nw,n_elec)
fragilmap[1:nw,1:n_elec]<-fragility[1:n_elec,1:nw]

# Library
library(ggplot2)
library(reshape)

colnames(fragilmap) <- paste("E", 1:n_elec)
rownames(fragilmap) <- paste("W", 1:30)

# Transform the matrix in long format
df <- melt(fragilmap)
colnames(df) <- c("x", "y", "value")
library(viridis)
ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile()+
  #scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  #scale_fill_distiller(palette = "Spectral", direction = -1) +
  #scale_fill_viridis_c(option = "B", direction = -1) +
  scale_fill_viridis(option="turbo")+
  coord_fixed()
