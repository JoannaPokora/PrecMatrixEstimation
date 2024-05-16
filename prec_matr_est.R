estimate_precision_matrix <- function(dat, penalty, treshold){
  # Estimate covariance matrix
  est.cov.matr <- apply(dat, 2, function(variable) variable-mean(variable))
  est.cov.matr <- 1/nrow(dat) *
    Reduce('+', apply(est.cov.matr, 1, function(obs.vec){
      obs.vec %*% t(obs.vec)
    }, simplify = FALSE))
  
  # Calculate penaltized estimated covariance matrix
  est.cov.pen.matr <- est.cov.matr + diag(penalty, nrow = ncol(dat))
  
  # Initialize beta estimator
  beta.hat <- list(old = 0, new = 0)
  # Graphical LASSO
  while(any(abs(beta.hat[["new"]]-beta.hat[["old"]]) > treshold) |
        all(beta.hat[["new"]] == 0)){
    sapply(1:ncol(dat), function(i){
      dif <- est.cov.matr[-i,i]- sapply(1:ncol(dat)-1, function(j){
        sum(est.cov.pen.matr[-i,-i][-j,j]*beta.hat[["new"]][i])
      })
      soft.tresh.op <- sign(dif)*pmax(abs(dif)-penalty, 0)
      beta.hat[["old"]] <<- beta.hat[["new"]]
      beta.hat[["new"]] <<- soft.tresh.op/diag(est.cov.pen.matr[-i,-i])
      
      # Update penaltizied covariance matrix
      est.cov.pen.matr[-i,i] <<- est.cov.pen.matr[-i,-i] %*%
        beta.hat[["new"]][-i]
    })
  }
  
  prec.matr <- diag(est.cov.pen.matr)
  sapply(1:ncol(dat), function(i){
    prec.matr[-i,i] <<- -beta.hat[["new"]][-i]*
      solve(est.cov.pen.matr[i,i]-
              t(est.cov.pen.matr[-i,i])%*%beta.hat[["new"]][-i])
  })
  
  prec.matr
}

dat <- matrix(c(4,5,3,8,2,1,5,6,3,9,4,2), ncol = 3)

estimate_precision_matrix(dat, 0, 0.00001)

          