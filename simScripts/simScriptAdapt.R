library(data.table)
library(hal9001)
library(sl3)
library(causalHAL)
library(doMC)
doMC::registerDoMC(cores = 11)

do_sims <- function(niter, n, pos_const, muIsHard, do_local_alt = FALSE) {
  sim_results <- rbindlist(lapply(1:niter, function(iter) {
    print(paste0("Iteration number: ", iter))
    try({
    data_list <- get_data(n, pos_const, muIsHard)
    return(as.data.table(get_estimates(data_list$W, data_list$A, data_list$Y,iter)))
    })
    return(data.table())
  }))
  key <- paste0("iter=", niter, "_n=", n, "_pos=", pos_const, "_hard=", muIsHard)
  fwrite(sim_results, paste0("~/causalHAL/simResults/sim_results_", key, ".csv"))
}


#' generates dataset of size n.
#' constant in propensity score can be used to vary overlap.
#' two settings for outcome regression: easy form and hard form
get_data <- function(n, pos_const, muIsHard = TRUE) {
  d <- 7
  W <- replicate(d, runif(n, -1, 1))
  colnames(W) <- paste0("W", 1:d)
  pi0 <- plogis(pos_const * ((1+ 1.5*W[,1])*sin(5*W[,1]) + (1 + 1.5*W[,2])* cos(5*W[,2]) + sin(3*W[,3]) + sin(5*W[,7]) + cos(3*W[,6])))
  A <- rbinom(n, 1, pi0)
  if(muIsHard) {
    mu0 <-  W[,1] + sin(5*W[,2]) + W[,3] + (1 + W[,4])*sin(5*W[,4]) + cos(5*W[,5])
  } else {
    mu0 <-  W[,1] + W[,2]  + W[,4] + W[,5]
  }
  tau <- 1 + W[,1] + W[,3] + W[,5]
  Y <- rnorm(n,  mu0 + A * tau, 0.5)
  return(list(W=W, A = A, Y = Y, ATE = 1))
}

get_data_local_alt <- function(n, pos_const, muIsHard = TRUE) {
  d <- 7
  W <- replicate(d, runif(n, -1, 1))
  colnames(W) <- paste0("W", 1:d)
  pi0 <- plogis(pos_const * ((1+ 1.5*W[,1])*sin(5*W[,1]) + (1 + 1.5*W[,2])* cos(5*W[,2]) + sin(3*W[,3]) + sin(5*W[,7]) + cos(3*W[,6])))
  A <- rbinom(n, 1, pi0)
  if(muIsHard) {
    mu0 <-  W[,1] + sin(5*W[,2]) + W[,3] + (1 + W[,4])*sin(5*W[,4]) + cos(5*W[,5])
  } else {
    mu0 <-  W[,1] + W[,2]  + W[,4] + W[,5]
  }
  # Add local alternative fluctuation
  mu0 <- mu0 + (1 + W[,3] + W[,4] + W[,6] + W[,7])/sqrt(n)
  tau <- 1 + W[,1] + W[,3] + W[,5]
  # Add local alternative fluctuation
  tau <- tau + (1 + W[,2] + W[,4] + W[,6] + W[,7])/sqrt(n)
  Y <- rnorm(n,  mu0 + A * tau, 0.5)
  return(list(W=W, A = A, Y = Y, ATE = 1))
}

#' Given simulated data (W,A,Y) and simulation iteration number `iter`,
#' computes ATE estimates, se, and CI for plug-in T-learner HAL, plug-in R-learner HAL, partially linear intercept model, AIPW.
get_estimates <- function(W, A, Y,iter) {
  n <- length(Y)
  if(n <= 500) {
    num_knots <- c(10, 5)
  } else if(n <= 1000) {
    num_knots <- c(25, 10)
  } else if(n <= 3000) {
    num_knots <- c(50, 15)
  } else{
    num_knots <- c(100, 30)
  }
  fit_T <- fit_hal_cate(W, A, Y,  max_degree = 1, num_knots = num_knots, smoothness_orders = 1,max_degree_cate =1, num_knots_cate = num_knots, smoothness_orders_cate = 1,    screen_variables = TRUE, fit_control = list(parallel = TRUE))
  ate_T <- unlist(inference_cate(fit_T))
  ate_T[1] <- "Tlearner"


  lrnr_A <- Lrnr_gam$new(family = "binomial")
  fit_R <- fit_hal_pcate(W, A, Y, lrnr_Y = NULL, lrnr_A = lrnr_A,  max_degree =2, num_knots = num_knots,  max_degree_cate = 1, num_knots_cate = num_knots, smoothness_orders_cate = 1,    screen_variables = TRUE, formula_cate = ~ h(.),fit_control = list(parallel = TRUE))
  ate_R<-  unlist(inference_cate(fit_R))
  ate_R[1] <- "Rlearner"

  lrnr_A <- fit_R$internal$fit_EAW
  #fit_R_intercept <- fit_hal_pcate(W, A, Y, lrnr_Y = NULL, lrnr_A = lrnr_A,  max_degree =2, num_knots = num_knots,  max_degree_cate = 1, num_knots_cate = 2, smoothness_orders_cate = 1,    screen_variables = TRUE, formula_cate = ~ h(W7), fit_control = list(parallel = TRUE))
  #mean(fit_R$internal$data$tau_relaxed)
  #ate_intercept <-  unlist(inference_cate(fit_R_intercept))
  #ate_intercept[1] <- "intercept"

  pi <- fit_R$internal$data$pi
  m <- fit_R$internal$data$mu
  tau_int <- mean((A-pi) * (Y - m)) / mean((A-pi)^2)
  IF <- (A - pi) / mean((A-pi)^2) * (Y - m - (A-pi)*tau_int)
  CI <- tau_int + 1.96*c(-1,1)*sd(IF)/sqrt(n)
  ate_intercept <- c("intercept", tau_int, sd(IF)/sqrt(n), CI)
  names(ate_intercept) <- c("method", "coef","se", "CI_left", "CI_right")




  mu1 <- fit_T$internal$data$mu1
  mu0 <- fit_T$internal$data$mu1
  mu <- ifelse(A==1, mu1, mu0)
  IF <- mu1 - mu0 + (A - pi) / ((1-pi)*pi) * (Y - mu)
  est_AIPW <-  mean(IF)
  CI <- est_AIPW + 1.96*c(-1,1)*sd(IF)/sqrt(n)
  ate_aipw <- c("AIPW", est_AIPW, sd(IF)/sqrt(n), CI)
  names(ate_aipw) <- c("method", "coef","se", "CI_left", "CI_right")

  mat <- cbind(iter,rbind(ate_T, ate_R, ate_intercept, ate_aipw))
  colnames(mat)  <- c("iter", "method", "coef","se", "CI_left", "CI_right")
  return(mat)
}
