set.seed(12)

library(CVglasso)
library(ggplot2)
library(MASS)
library(cvCovEst)
library(dplyr)
library(mvtnorm)
library(Matrix)
library(matrixcalc)

estimate_covariance_matrix <- function(dat){
  sample.mean <- apply(dat, 2, mean)
  est.cov.matr <- 1/nrow(dat) *
    Reduce('+', apply(dat, 1, function(sample){
      scaled.sample <- sample - sample.mean
      scaled.sample %*% t(scaled.sample)
    }, simplify = FALSE))
  est.cov.matr
}

glasso <- function(dat, penalty, treshold = 1e-4){
  est.cov.matr <- estimate_covariance_matrix(dat)
  
  # Calculate inverse precision matrix
  inv.prec.matr <- est.cov.matr + diag(penalty, nrow = ncol(dat))
  inv.prec.matr.old <- inv.prec.matr
  
  # Initialize beta estimator
  beta.hat <- list(old = rep(0, ncol(dat)),
                   new = rep(0, ncol(dat)))
  
  iter <- 0
  
  # Perform lasso algorithm
  while(((any(abs(as.vector(inv.prec.matr-inv.prec.matr.old)) > treshold) |
        iter == 0) & iter < 100)){
    iter <- iter + 1
    inv.prec.matr.old <- inv.prec.matr
    for(i in 1:ncol(dat)){
      crop.beta.hat <- beta.hat[["new"]][-i]
      crop.est.cov.matr <- est.cov.matr[-i,i]
      crop.inv.prec.matr <- inv.prec.matr[-i,-i]
      
      iter.beta <- 0
      while((any(abs(crop.beta.hat-beta.hat[["old"]][-i]) > treshold) |
            iter.beta == 0) & iter.beta < 100){
        iter.beta <- iter.beta + 1
        beta.hat[["old"]][-i] <- crop.beta.hat
        for(j in 1:(ncol(dat)-1)){
          dif <- crop.est.cov.matr[j]-
            sum(crop.inv.prec.matr[-j,j]*crop.beta.hat[-j])
          soft.tresh.op <- sign(dif)*max(abs(dif)-penalty,0)
          crop.beta.hat[j] <- soft.tresh.op/crop.inv.prec.matr[j,j]
        }
      beta.hat[["new"]][-i] <- crop.beta.hat
      
      # Update inverse precision matrix
      inv.prec.matr[-i,i] <- inv.prec.matr[-i,-i] %*%
        beta.hat[["new"]][-i]
     }
    }
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

# Ostateczne symulacje
cases <- expand.grid(Methood = c("glasso", "CVglasso", "shrinkage"),
                     Distribution = c("normal", "t-distribution"),
                     rep = 1:10,
                     obs.num = c(100, 500),
                     var.num = c(100, 200),
                     cor.lvl = c("low", "high"),
                     cor.var.num = c(5, 50, "all"),
                     stringsAsFactors = FALSE) %>%
  mutate(cor.var.num = ifelse(cor.var.num == "all", as.numeric(var.num),
                              as.numeric(cor.var.num)))

prec.matrs <- apply(unique(cases[,4:7]), 1, function(case){
  var.num <- as.numeric(case["var.num"])
  cor.var.num <- as.numeric(case["cor.var.num"])
  if(case["cor.lvl"] == "low") cor.lvl <- c(0.3, 0.5)
  else cor.lvl <- c(0.7, 0.9)
  
  coord.choose <- diag(0, nrow = var.num)
  coord.choose[upper.tri(coord.choose)][sample(1:((var.num^2-var.num)/2),
                                               cor.var.num)] <- 1
  coord <- which(coord.choose == 1, arr.ind = TRUE)

  prec.matr <- as.matrix(sparseMatrix(coord[,1], coord[,2],
                         x = runif(cor.var.num, cor.lvl[1], cor.lvl[2]),
                         dims = c(var.num, var.num),
                         symmetric = TRUE))
  diag(prec.matr) <- var.num
  prec.matr
})

cases <- cbind(cases, Times = 1:nrow(cases),
               MSE = 1:nrow(cases), DifSup = 1:nrow(cases))

cases.dat <- unique(cases[3:7])
dats <- lapply(1:(nrow(cases.dat)), function(case){
  prec.matr <- prec.matrs[[(case-1)%/%20+1]]
  list(mvrnorm(cases.dat[case, 2], rep(0, nrow(prec.matr)), solve(prec.matr)),
       rmvt(cases.dat[case, 2], solve(prec.matr), df = 10, checkSymmetry = FALSE))
})
dats <- unlist(dats, recursive = FALSE)

calc_stat <- function(est.prec.matr, prec.matr){
  dif <- est.prec.matr-prec.matr
  c(mse = mean((dif)^2), dif.sup = max(abs(dif)))
}

for(case in 1:length(dats)){
   prec.matr <- prec.matrs[[(case-1)%/%40+1]]

  CVglasso.start <- Sys.time()
  CVglasso.res <- CVglasso(X = dats[[case]])
  cases[case*3-1, 8] <- Sys.time() - CVglasso.start

  CVglasso.matr <- CVglasso.res[["Omega"]]
  CVglasso.lambda <- CVglasso.res[["Tuning"]][2]

  cases[case*3-1, 9:10] <- calc_stat(CVglasso.matr, prec.matr)

  glasso.start <- Sys.time()
  glasso.matr <- glasso(dats[[case]], CVglasso.lambda)
  cases[case*3-2, 8] <- Sys.time() - glasso.start

  cases[case*3-2, 9:10] <- calc_stat(glasso.matr, prec.matr)
  
  shrink.start <- Sys.time()
  shrink.matr <- solve(linearShrinkLWEst(dats[[case]]))
  cases[case*3, 8] <- Sys.time() - shrink.start
  
  cases[case*3, 9:10] <- calc_stat(shrink.matr, prec.matr)
}

results <- cases %>%
  group_by(Methood, Distribution, obs.num, var.num, cor.var.num, cor.lvl) %>%
  summarise(Mean.time = mean(Times),
            Mean.MSE = mean(MSE),
            Mean.dif.sup = mean(DifSup)) %>%
  ungroup()

ggplot(results, aes(x = cor.var.num, y = Mean.time, color = Distribution)) +
  geom_point() +
  geom_line() +
  facet_grid(obs.num + var.num ~ Methood + cor.lvl)

ggplot(results, aes(x = cor.var.num, y = Mean.MSE, color = Distribution)) +
  geom_point() +
  geom_line() +
  facet_grid(obs.num + var.num ~ Methood + cor.lvl)

ggplot(results, aes(x = cor.var.num, y = Mean.dif.sup, color = Distribution)) +
  geom_point() +
  geom_line() +
  facet_grid(obs.num + var.num ~ Methood + cor.lvl)if.sup, color = Distribution)) +
  geom_point() +
  geom_line() +
  facet_grid(obs.num + var.num ~ Methood + rot)
