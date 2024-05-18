estimate_precision_matrix <- function(dat, penalty, treshold){
  # Estimate covariance matrix
  sample.mean <- apply(dat, 2, mean)
  est.cov.matr <- 1/nrow(dat) *
    Reduce('+', apply(dat, 1, function(sample){
      scaled.sample <- sample - sample.mean
      scaled.sample %*% t(scaled.sample)
    }, simplify = FALSE))
  
  # Calculate inverse precision matrix
  inv.prec.matr <- est.cov.matr + diag(penalty, nrow = ncol(dat))
  inv.prec.matr.old <- inv.prec.matr
  
  # Initialize beta estimator
  beta.hat <- list(old = rep(0, ncol(dat)),
                   new = rep(0, ncol(dat)))
  
  # Perform lasso algorithm
  while(any(abs(as.vector(inv.prec.matr-inv.prec.matr.old)) > treshold) |
        all(as.vector(inv.prec.matr == inv.prec.matr.old))){
    inv.prec.matr.old <- inv.prec.matr
    sapply(1:ncol(dat), function(i){
      crop.beta.hat <- beta.hat[["new"]][-i]
      crop.est.cov.matr <- est.cov.matr[-i,i]
      crop.inv.prec.matr <- inv.prec.matr[-i,-i]
      while(any(abs(crop.beta.hat-beta.hat[["old"]][-i]) > treshold) |
            all(beta.hat[["new"]] == beta.hat[["old"]])){
        print(abs(crop.beta.hat-beta.hat[["old"]][-i]))
        beta.hat[["old"]][-i] <<- crop.beta.hat
        sapply(1:(ncol(dat)-1), function(j){
          dif <- crop.est.cov.matr[j]-
            sum(crop.inv.prec.matr[-j,j]*crop.beta.hat[-j])
          soft.tresh.op <- sign(dif)*max(abs(dif)-penalty,0)
          crop.beta.hat[j] <<- soft.tresh.op/crop.inv.prec.matr[j,j]
        })
      }
      beta.hat[["new"]][-i] <<- crop.beta.hat
      
      # Update inverse precision matrix
      inv.prec.matr[-i,i] <<- inv.prec.matr[-i,-i] %*%
        beta.hat[["new"]][-i]
    })
  }
  
  # Calculate precision matrix
  prec.matr <- matrix(0, nrow = ncol(dat), ncol = ncol(dat))
  sapply(1:ncol(dat), function(i){
    prec.matr[i,i] <<- solve(inv.prec.matr[i,i]-
                               t(inv.prec.matr[-i,i])%*%beta.hat[["new"]][-i])
    prec.matr[-i,i] <<- -beta.hat[["new"]][-i]*prec.matr[i,i]
  })
  
  prec.matr
}

dat <- matrix(c(4,5,3,8,2,1,5,6,3,9,4,2), ncol = 3)

estimate_precision_matrix(dat, 0, 0.001)
