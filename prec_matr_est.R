set.seed(12)

library(CVglasso)
library(ggplot2)
library(MASS)
library(cvCovEst)
library(dplyr)
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(hrbrthemes)

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
results <- expand.grid(Method = c("glasso", "CVglasso", "shrinkage"),
                     Distribution = c("normal", "t-distribution"),
                     rep = 1:10,
                     obs.num = c(100, 500),
                     var.num = c(100, 200),
                     cov.val = c("0.3-0.5", "0.7-0.9"),
                     cor.var.num = c(5, 50, "all"),
                     stringsAsFactors = FALSE) %>%
  mutate(cor.var.num = ifelse(cor.var.num == "all", as.numeric(var.num),
                              as.numeric(cor.var.num)))

prec.matrs <- apply(unique(results[,4:7]), 1, function(result){
  var.num <- as.numeric(result["var.num"])
  cor.var.num <- as.numeric(result["cor.var.num"])
  if(result["cov.val"] == "0.3-0.5") cov.val <- c(0.3, 0.5)
  else cov.val <- c(0.7, 0.9)
  
  coord.choose <- diag(0, nrow = var.num)
  coord.choose[upper.tri(coord.choose)][sample(1:((var.num^2-var.num)/2),
                                               cor.var.num)] <- 1
  coord <- which(coord.choose == 1, arr.ind = TRUE)

  prec.matr <- as.matrix(sparseMatrix(coord[,1], coord[,2],
                         x = runif(cor.var.num, cov.val[1], cov.val[2]),
                         dims = c(var.num, var.num),
                         symmetric = TRUE))
  diag(prec.matr) <- var.num
  prec.matr
})

results <- cbind(results, Times = 1:nrow(results),
               MSE = 1:nrow(results), DifSup = 1:nrow(results))

results.dat <- unique(results[3:7])
dats <- lapply(1:(nrow(results.dat)), function(result){
  prec.matr <- prec.matrs[[(result-1)%/%20+1]]
  list(mvrnorm(results.dat[result, 2], rep(0, nrow(prec.matr)), solve(prec.matr)),
       rmvt(results.dat[result, 2], solve(prec.matr), df = 10, checkSymmetry = FALSE))
})
dats <- unlist(dats, recursive = FALSE)

calc_stat <- function(est.prec.matr, prec.matr){
  dif <- est.prec.matr-prec.matr
  c(mse = mean((dif)^2), dif.sup = max(abs(dif)))
}

for(result in 1:length(dats)){
   prec.matr <- prec.matrs[[(result-1)%/%40+1]]

  CVglasso.start <- Sys.time()
  CVglasso.res <- CVglasso(X = dats[[result]])
  results[result*3-1, 8] <- Sys.time() - CVglasso.start

  CVglasso.matr <- CVglasso.res[["Omega"]]
  CVglasso.lambda <- CVglasso.res[["Tuning"]][2]

  results[result*3-1, 9:10] <- calc_stat(CVglasso.matr, prec.matr)

  glasso.start <- Sys.time()
  glasso.matr <- glasso(dats[[result]], CVglasso.lambda)
  results[result*3-2, 8] <- Sys.time() - glasso.start

  results[result*3-2, 9:10] <- calc_stat(glasso.matr, prec.matr)
  
  shrink.start <- Sys.time()
  shrink.matr <- solve(linearShrinkLWEst(dats[[result]]))
  results[result*3, 8] <- Sys.time() - shrink.start
  
  results[result*3, 9:10] <- calc_stat(shrink.matr, prec.matr)
}

results <- results %>%
  mutate(obs.num = paste0(obs.num, " observations"),
         var.num = paste0(var.num, " variables"))

plots_theme <- function(){
  theme
}

results %>%
  group_by(Method, Distribution, obs.num, var.num, cor.var.num, cov.val) %>%
  summarise(Mean.time = mean(Times)) %>%
  ggplot(aes(x = cor.var.num, y = Mean.time, color = Method)) +
  geom_point() +
  geom_line() +
  labs(x = "Number of conditional dependent variables",
       y = "Mean time") +
  scale_color_manual(values = ) +
  facet_grid(cov.val + Distribution ~ var.num + obs.num, scales = "free") +
  theme_light()

ggplot(results, aes(x = cor.var.num, y = MSE, color = Distribution)) +
  geom_boxplot() +
  facet_grid(obs.num + var.num ~ Method + cov.val)

ggplot(results, aes(x = cor.var.num, y = Mean.dif.sup, color = Distribution)) +
  geom_point() +
  geom_line() +
  facet_grid(obs.num + var.num ~ Methood + cov.val)
