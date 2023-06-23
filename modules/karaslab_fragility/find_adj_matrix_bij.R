checkstability <- function(V,UtXtY,svd,scale_x,scale_y,lam){

  adj_matrix <- V %*% diag(1 / (svd$d + lam)) %*% UtXtY

  # right now scale(x) %*% adj_matrix = scale(y), need to scale back
  adj_matrix <- diag( 1 / scale_x ) %*% adj_matrix %*% diag( scale_y )


  eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
  stable <- max(eigv) < 1
  return(stable)
}

find_adj_matrix_bij <- function(x, y, a=0.00001, b=100, n=20, tol=1e-5) {

  # x <- t(state_vectors$x)
  # y <- t(state_vectors$x_n)

  nobs <- nrow(x)
  nvars <- ncol(x)

  # scale x and y
  x <- scale(x)
  y <- scale(y)
  scale_x <- attr(x, "scaled:scale")
  scale_y <- attr(y, "scaled:scale")

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

  # check lambda range
  stablea <- checkstability(V,UtXtY,svd,scale_x,scale_y,a)
  stableb <- checkstability(V,UtXtY,svd,scale_x,scale_y,b)

  if ((stablea==stableb)|(stableb==FALSE)) {
    print(stablea)
    print(stableb)
    stop(' Range is not appropriate to find optimal lambda')
  } else if ((stablea==FALSE)&(stableb==TRUE)) {
    print(' Range is appropriate to find optimal lambda')
  }

  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint
    stablec <- checkstability(V,UtXtY,svd,scale_x,scale_y,c)

    if(stablec==TRUE){
      b=c
    }else{
      a=c
    }
  }

  print(b)

  adj_matrix <- V %*% diag(1 / (svd$d + b)) %*% UtXtY

  # right now scale(x) %*% adj_matrix = scale(y), need to scale back
  adj_matrix <- diag( 1 / scale_x ) %*% adj_matrix %*% diag( scale_y )


  eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
  stable <- max(eigv) < 1

  structure(adj_matrix)
  return(adj_matrix)

}
