set.seed(12)

library(glasso)
library(ggplot2)
library(MASS)
library(pracma)
library(cvCovEst)
library(dplyr)

estimate_covariance_matrix <- function(dat){
  sample.mean <- apply(dat, 2, mean)
  est.cov.matr <- 1/(nrow(dat) - 1) *
    Reduce('+', apply(dat, 1, function(sample){
      scaled.sample <- sample - sample.mean
      scaled.sample %*% t(scaled.sample)
    }, simplify = FALSE))
  est.cov.matr
}

estimate_precision_matrix <- function(dat, penalty, treshold = 1e-4){
  est.cov.matr <- estimate_covariance_matrix(dat)
  
  # Calculate inverse precision matrix
  inv.prec.matr <- est.cov.matr + diag(penalty, nrow = ncol(dat))
  inv.prec.matr.old <- inv.prec.matr
  
  # Initialize beta estimator
  beta.hat <- list(old = rep(0, ncol(dat)),
                   new = rep(0, ncol(dat)))
  
  iter <- 0
  
  # Perform lasso algorithm
  while((any(abs(as.vector(inv.prec.matr-inv.prec.matr.old)) > treshold) |
        iter == 0) & iter < 100){
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
      }
      beta.hat[["new"]][-i] <- crop.beta.hat
      
      # Update inverse precision matrix
      inv.prec.matr[-i,i] <- inv.prec.matr[-i,-i] %*%
        beta.hat[["new"]][-i]
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

# penalty <- c(1, 2, 5, 8, 10)
# treshold <- 1e-4
# 
# variables.num <- c(5, 100, 200)
# 
# diag.prec.matr <- lapply(variables.num, function(k){
#   diag(runif(k, 1, 5))
# })
# 
# low.prec.matr <- lapply(diag.prec.matr, function(diagonal){
#   diagonal[diagonal == 0] <- 0.1
#   diagonal
# })
# 
# high.prec.matr <- lapply(diag.prec.matr, function(diagonal){
#   diagonal[diagonal == 0] <- 0.8
#   diagonal
# })
# 
# sparse.prec.matr <- lapply(diag.prec.matr, function(diagonal){
#   rand.cov <- sample(c(0, 0.4),
#                      (ncol(diagonal)^2-ncol(diagonal))/2,
#                      replace = TRUE,
#                      prob = c(0.8, 0.2))
#   diagonal[lower.tri(diagonal)] <- rand.cov
#   diagonal %*% t(diagonal)
# })
# 
# # rand.prec.matr <- lapply(diag.prec.matr, function(diagonal){
# #   diagonal[lower.tri(diagonal)] <- rnorm((nrow(diagonal)^2-nrow(diagonal))/2)
# #   diagonal %*% t(diagonal)
# # })
# 
# all.prec.matr <- list(diag.prec.matr = diag.prec.matr,
#                       low.prec.matr = low.prec.matr,
#                       high.prec.matr = high.prec.matr,
#                       sparse.prec.matr = sparse.prec.matr)
# 
# dat <- lapply(all.prec.matr, function(case){
#   lapply(case, function(prec.matr){
#     MASS::mvrnorm(100, mu = rep(0, ncol(prec.matr)), Sigma = solve(prec.matr))
#   })
# })
# 
# results <- setNames(lapply(penalty, function(lambda){
#   setNames(lapply(dat, function(case){
#     setNames(lapply(case, function(dat){
#         est.prec.matr.start <- Sys.time()
#         est.prec.matr <- estimate_precision_matrix(dat, lambda, treshold)
#         est.prec.matr.time <- Sys.time() - est.prec.matr.start
#         
#         glasso.start <- Sys.time()
#         est.cov.matr <- estimate_covariance_matrix(dat)
#         glasso.matr <- glasso(est.cov.matr, lambda)
#         glasso.time <- Sys.time() - glasso.start
#       
#         list(est.prec.matr = list(matr = est.prec.matr,
#                                   time = est.prec.matr.time),
#              glasso.matr = list(matr = glasso.matr[["wi"]],
#                                 time = glasso.time))
#     }), paste0("variables=", variables.num))
#   }), paste0("case=", names(all.prec.matr)))
# }), paste0("penalty=", penalty))
# 
# MSE <- lapply(results, function(penalty){
#   mapply(function(res.case, prec.matr.case){
#     mapply(function(res.meth, prec.matr){
#       lapply(res.meth, function(res){
#         mean((res[[1]]-prec.matr)^2)
#       })
#     }, res.case, prec.matr.case, SIMPLIFY = FALSE)
#   }, penalty, all.prec.matr, SIMPLIFY = FALSE)
# })
# 
# MSE.df <- data.frame(MSE = unlist(MSE),
#                       Penalty = rep(penalty, each = 24),
#                       Case = rep(names(all.prec.matr), each = 6),
#                       Variables.num = rep(variables.num, each = 2),
#                       Method = c("est.cov.matr", "glasso"))
# 
# ggplot(MSE.df, aes(x = Penalty, y = MSE, color = Method)) +
#   geom_point() +
#   geom_line() +
#   facet_grid(Case ~ Variables.num, scales = "free_y")
# 
# false.disc <- lapply(results, function(penalty){
#   mapply(function(res.case, prec.matr.case){
#     mapply(function(res.meth, prec.matr){
#       lapply(res.meth, function(res){
#         res[[1]][which(prec.matr == 0)]
#       })
#     }, res.case, prec.matr.case, SIMPLIFY = FALSE)
#   }, penalty, all.prec.matr, SIMPLIFY = FALSE)
# })
# 
# times <- lapply(results, function(penalty){
#   lapply(penalty, function(case){
#     lapply(case, function(variables.num){
#       lapply(variables.num, function(result){
#         result[[2]]
#       })
#     })
#   })
# })
# times.df <- data.frame(Time = unlist(times),
#                        Penalty = rep(penalty, each = 24),
#                        Case = rep(names(all.prec.matr), each = 6),
#                        Variables.num = rep(variables.num, each = 2),
#                        Method = c("est.cov.matr", "glasso"))
# 
# ggplot(times.df, aes(x = Penalty, y = Time, color = Method)) +
#   geom_point() +
#   geom_path() +
#   facet_grid(Case ~ Variables.num)

# Lambda

rep.times <- 10
lambda.vals <- c(1, 2, 5)
var.nums <- c(5, 30, 50, 100)
obs.num <- 50

cov.matrs <- lapply(var.nums, function(var.num){
  rand.matr <- diag(runif(var.num, 5, 10))
  rand.matr[lower.tri(rand.matr)] <- runif((var.num^2-var.num)/2, 1, 10)
  rand.matr %*% t(rand.matr)
})

rep.dat <- lapply(1:rep.times, function(i){
  lapply(cov.matrs, function(cov.matr){
    mvrnorm(obs.num, rep(0, nrow(cov.matr)), cov.matr)
  })
})

res <- lapply(rep.dat, function(i){
  lapply(i, function(dat){
    lapply(lambda.vals, function(lambda){
      # est.prec.matr.start <- Sys.time()
      # est.prec.matr <- estimate_precision_matrix(dat, lambda)
      # est.prec.matr.time <- Sys.time() - est.prec.matr.start
      #            
      glasso.start <- Sys.time()
      est.cov.matr <- estimate_covariance_matrix(dat)
      glasso.matr <- glasso(est.cov.matr, lambda)
      glasso.time <- Sys.time() - glasso.start
               
      list(glasso.matr = list(prec.matr = glasso.matr[["wi"]],
                              time = glasso.time))
      # est.prec.matr = list(prec.matr = est.prec.matr,
      #                      time = est.prec.matr.time),
    })
  })
})

false.disc <- lapply(res, function(i){
     mapply(function(res.val.num, cov.matr){
       lapply(res.val.num, function(res.lambda){
         lapply(res.lambda, function(res){
           sum(res[[1]][which(solve(cov.matr) == 0)] != 0)
         })
       })
     }, i, cov.matrs, SIMPLIFY = FALSE)
   })

true.disc.prop <- lapply(res, function(i){
  mapply(function(res.var.num, cov.matr){
    lapply(res.var.num, function(res.lambda){
      lapply(res.lambda, function(res){
        prec.matr <- solve(cov.matr)
        sum(res[[1]][which(prec.matr != 0)] != 0)/ sum(prec.matr != 0)
      })
    })
  }, i, cov.matrs, SIMPLIFY = FALSE)
})

true.disc.prop.df <- data.frame(TD.prop = unlist(true.disc.prop),
                                Var.num = rep(var.nums, each = rep.times),
                                Lambda =
                                  rep(lambda.vals,
                                      each = rep.times*length(var.nums))) %>%
  group_by(Var.num, Lambda) %>%
  summarise(TD.prop = mean(TD.prop))

ggplot(true.disc.prop.df, aes(x = Var.num, y = TD.prop,
                              color = as.factor(Lambda))) +
  geom_point() +
  geom_line()

shrink.res <- lapply(rep.dat, function(i){
  lapply(i, function(dat){
      shrink.start <- Sys.time()
      shrink.matr <- linearShrinkLWEst(dat)
      shrink.time <- Sys.time() - shrink.start
      
      list(prec.matr = solve(shrink.matr),
           time = shrink.time)
  })
})

shrink.true.disc.prop <- lapply(shrink.res, function(i){
  mapply(function(res, cov.matr){
    prec.matr <- solve(cov.matr)
    sum(res[[1]][which(prec.matr != 0)] != 0)/ sum(prec.matr != 0)
  }, i, cov.matrs, SIMPLIFY = FALSE)
})

shrink.true.disc.prop.df <- data.frame(TD.prop = unlist(shrink.true.disc.prop),
                                Var.num = rep(var.nums, each = rep.times)) %>%
  group_by(Var.num) %>%
  summarise(TD.prop = mean(TD.prop))

ggplot(shrink.true.disc.prop.df, aes(x = Var.num, y = TD.prop)) +
  geom_point() +
  geom_line()
