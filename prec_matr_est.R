set.seed(12)

library(CVglasso)
library(ggplot2)
library(MASS)
library(cvCovEst)
library(dplyr)
library(mvtnorm)
library(Matrix)
library(matrixcalc)
library(ggpubr)

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
      while((any(abs(crop.beta.hat-beta.hat[["old"]][-i]) > 0.01) |
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
rep.times <- 10

results <- expand.grid(Method = c("glasso", "CVglasso", "shrinkage"),
                     Distribution = c("normal distribution", "t-distribution"),
                     rep = 1:rep.times,
                     obs.num = c(100, 500),
                     var.num = c(100, 200),
                     cor.var.num = c("5", "50", "all"),
                     stringsAsFactors = FALSE)

prec.matrs <- apply(unique(results[,4:6]), 1, function(result){
  var.num <- as.numeric(result["var.num"])
  
  if(result["cor.var.num"] == "all") cor.var.num <- (var.num^2-var.num)/2
  else cor.var.num <- as.numeric(result["cor.var.num"])
  
  coord.choose <- diag(0, nrow = var.num)
  coord.choose[upper.tri(coord.choose)][sample(1:((var.num^2-var.num)/2),
                                               cor.var.num)] <- 1
  coord <- which(coord.choose == 1, arr.ind = TRUE)

  prec.matr <- as.matrix(sparseMatrix(coord[,1], coord[,2],
                         x = rnorm(cor.var.num),
                         dims = c(var.num, var.num),
                         symmetric = TRUE))
  diag(prec.matr) <- var.num
  prec.matr
})

results <- cbind(results, Times = 1:nrow(results),
                 sq.frob.norm = 1:nrow(results),
                 DifSup = 1:nrow(results),
                 power = 1:nrow(results),
                 fdp = 1:nrow(results))

results.dat <- unique(results[3:6])
dats <- lapply(1:(nrow(results.dat)), function(result){
  prec.matr <- prec.matrs[[(result-1)%/%rep.times+1]]
  list(mvrnorm(results.dat[result, 2],
               rep(0, nrow(prec.matr)),
               solve(prec.matr)),
       rmvt(results.dat[result, 2],
            solve(prec.matr), df = 5,
            checkSymmetry = FALSE))
})

dats <- unlist(dats, recursive = FALSE)

calc_stat <- function(est.prec.matr, prec.matr){
  dif <- est.prec.matr-prec.matr
  c(matrix.mse = mean((dif%*%t(dif))^2), dif.sup = max(abs(dif)),
    power = sum(est.prec.matr[which(prec.matr != 0)] != 0)/sum(prec.matr != 0),
    fdp = sum(est.prec.matr[which(prec.matr == 0)] != 0))/
          max(sum(est.prec.matr != 0), 1)
}

est.prec.matrs <- as.list(rep(0, nrow(results)))

for(result in 1:length(dats)){
  print(result)
  prec.matr <- prec.matrs[[(result-1)%/%(2*rep.times)+1]]

  est.prec.matrs[[result*3-1]] <- CVglasso(X = dats[[result]])
  results[result*3-1, 7] <- est.prec.matrs[[result*3-1]][["Time"]]

  CVglasso.matr <- est.prec.matrs[[result*3-1]][["Omega"]]
  CVglasso.lambda <- est.prec.matrs[[result*3-1]][["Tuning"]][2]

  results[result*3-1, 8:11] <- calc_stat(CVglasso.matr, prec.matr)

  glasso.res <- tryCatch({
    glasso.start <- Sys.time()
    glasso.res <- glasso(dats[[result]], CVglasso.lambda)
    glasso.time <- Sys.time() - glasso.start
    list(glasso.res, glasso.time)
  }, error = function(e){
    print("error")
    NaN
  })
  
  if(any(is.na(glasso.res))){
    results[result*3-2, 7:11] <- NaN
    est.prec.matrs[[result*3-2]] <- NaN
  }else{
    results[result*3-2, 7] <- glasso.res[[2]]
    est.prec.matrs[[result*3-2]] <- glasso.res[[1]]
    results[result*3-2, 8:11] <- calc_stat(est.prec.matrs[[result*3-2]],
                                          prec.matr)
  }
  
  shrink.start <- Sys.time()
  est.prec.matrs[[result*3]] <- solve(linearShrinkLWEst(dats[[result]]))
  results[result*3, 7] <- Sys.time() - shrink.start
  
  results[result*3, 8:11] <- calc_stat(est.prec.matrs[[result*3]], prec.matr)
}

results <- results %>%
  mutate(obs.num = paste0(obs.num, " observations"),
         var.num = paste0(var.num, " variables"))

results <- results %>%
  filter(!(is.na(Times)))

results %>%
  group_by(Method, Distribution, obs.num, var.num, cor.var.num) %>%
  summarise(Mean.time = mean(Times)) %>%
  ggplot(aes(x = cor.var.num, y = Mean.time, color = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(x = "Number of nonzero elements in a precision matrix",
       y = "Mean time") +
  scale_color_manual(values = c("CVglasso" = "#FFB909",
                                "glasso" = "#30CB14",
                                "shrinkage" = "#0b2206")) +
  facet_wrap(Distribution ~ var.num + obs.num, scales = "free_y", nrow = 2) +
  theme_light()

results %>%
  group_by(Method, Distribution, obs.num, var.num, cor.var.num) %>%
  summarise(MSE = mean(sq.frob.norm)) %>%
  ggplot(aes(x = cor.var.num, y = MSE, color = Method, fill = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(x = "Number of nonzero elements in a precision matrix",
       y = "Mean squared error") +
  scale_color_manual(values = c("CVglasso" = "#FFB909",
                                "glasso" = "#30CB14",
                                "shrinkage" = "#0b2206")) +
  scale_fill_manual(values = c("CVglasso" = "#FFEAB5",
                               "glasso" = "#BDFCB1",
                               "shrinkage" = "#6A8663")) +
  facet_wrap(Distribution ~ var.num + obs.num, scales = "free_y", nrow = 2) +
  theme_light()

ggplot(results, aes(x = cor.var.num, y = DifSup, color = Method, fill = Method)) +
  geom_boxplot() +
  labs(x = "Number of nonzero elements in a precision matrix",
       y = "Maximum difference") + 
  scale_color_manual(values = c("CVglasso" = "#FFB909",
                                "glasso" = "#30CB14",
                                "shrinkage" = "#0b2206")) +
  scale_fill_manual(values = c("CVglasso" = "#FFEAB5",
                               "glasso" = "#BDFCB1",
                               "shrinkage" = "#6A8663")) +
  facet_wrap(Distribution ~ var.num + obs.num, scales = "free_y", nrow = 2) +
  theme_light()

results %>%
  group_by(Method, Distribution, obs.num, var.num, cor.var.num) %>%
  summarise(Mean.power = mean(power)) %>%
  ggplot(aes(x = cor.var.num, y = Mean.power, color = Method, fill = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(x = "Number of nonzero elements in a precision matrix",
       y = "Mean power") +
  scale_color_manual(values = c("CVglasso" = "#FFB909",
                                "glasso" = "#30CB14",
                                "shrinkage" = "#0b2206")) +
  scale_fill_manual(values = c("CVglasso" = "#FFEAB5",
                               "glasso" = "#BDFCB1",
                               "shrinkage" = "#6A8663")) +
  facet_wrap(Distribution ~ var.num + obs.num, nrow = 2) +
  theme_light()

results %>%
  group_by(Method, Distribution, obs.num, var.num, cor.var.num) %>%
  summarise(FDR = mean(fdp)) %>%
  ggplot(aes(x = cor.var.num, y = FDR, color = Method, fill = Method)) +
  geom_point() +
  geom_line(aes(group = Method)) +
  labs(x = "Number of nonzero elements in a precision matrix",
       y = "False discovery rate") +
  scale_color_manual(values = c("CVglasso" = "#FFB909",
                                "glasso" = "#30CB14",
                                "shrinkage" = "#0b2206")) +
  scale_fill_manual(values = c("CVglasso" = "#FFEAB5",
                               "glasso" = "#BDFCB1",
                               "shrinkage" = "#6A8663")) +
  facet_wrap(Distribution ~ var.num + obs.num, nrow = 2) +
  theme_light()
