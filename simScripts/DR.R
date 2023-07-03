library(data.table)
library(sl3)
library(doFuture)
library(future)

doFuture::registerDoFuture()
future::plan(multisession, workers = 7)


do_sims <- function(n, pos_const, nsims, misp) {
  # true ATE
true <- 0.8082744
# crossfit estimator
lrnr_cv <- Pipeline$new(Lrnr_cv$new(Stack$new(Lrnr_gam$new(), Lrnr_earth$new(degree = 1), Lrnr_ranger$new(), Lrnr_xgboost$new(max_depth = 5))), Lrnr_cv_selector$new(loss_squared_error))
lrnr_pi <- lrnr_mu <- lrnr_cv
lrnr_misp <- Lrnr_mean$new()
# 1 is both correct, 2 is just outcome, 3 is just treatment, 4 is neither
if(misp==2 || misp == 4) {
  lrnr_pi <- lrnr_misp
}
if(misp==3 || misp == 4) {
  lrnr_mu <- lrnr_misp
}
 sim_results <- lapply(1:nsims, function(i){
  try({
  print(paste0("iter: ", i))
  data_list <- get_data(n, pos_const)
  W <- data_list$W
  A <- data_list$A
  Y <- data_list$Y
  initial_estimators <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu, lrnr_pi = lrnr_pi)
  print("initial")
  folds <- initial_estimators$folds
  out_AIPW <- compute_AIPW(A,Y, initial_estimators$mu1, initial_estimators$mu0, initial_estimators$pi1, initial_estimators$pi0)
  out_AuDRIE <- compute_AuDRIE_boot(A,Y, initial_estimators$mu1, initial_estimators$mu0, initial_estimators$pi1, initial_estimators$pi0, nboot = 5000, folds = folds, alpha = 0.05)
  out <- unlist(c(out_AuDRIE, out_AIPW))
  names(out) <- c("estimate_audrie", "CI_left_audrie", "CI_right_audrie", "estimate_AIPW", "CI_left_AIPW", "CI_right_AIPW")
  return(out)
  })
  return(data.table())
})
sim_results <- data.table::rbindlist(sim_results)
key <- paste0("DR_iter=", nsims, "_n=", n, "_pos=", pos_const, "_mode=", misp )
try({fwrite(sim_results, paste0("~/causalHAL/simResultsDR/sim_results_", key, ".csv"))})
return(sim_results)
}



get_data <- function(n, pos_const) {
  d <- 4
  W <- replicate(d, runif(n, -1, 1))
  colnames(W) <- paste0("W", 1:d)
  pi0 <- plogis(pos_const * ( W[,1] + sin(4*W[,1]) +   W[,2] + cos(4*W[,2]) + W[,3] + sin(4*W[,3]) + W[,4] + cos(4*W[,4]) ))
  A <- rbinom(n, 1, pi0)
  mu0 <-  sin(4*W[,1]) + sin(4*W[,2]) + sin(4*W[,3])+  sin(4*W[,4]) + cos(4*W[,2])
  tau <- 1 + W[,1] + sin(4*W[,2]) + cos(4*W[,3]) + W[,4]
  Y <- rnorm(n,  mu0 + A * tau, 0.5)
  return(list(W=W, A = A, Y = Y, ATE = 0.8082744, pi = pi0))
}



compute_AuDRIE_boot <-  function(A,Y, mu1, mu0, pi1, pi0, nboot = 5000, folds, alpha = 0.05) {
  data <- data.table(A, Y, mu1, mu0, pi1, pi0)
  folds <- lapply(folds, `[[`, "validation_set")
  tau_n <- compute_AuDRIE(A,Y, mu1, mu0, pi1, pi0)

  bootstrap_estimates <- sapply(1:nboot, function(iter){
    try({
      bootstrap_indices <- unlist(lapply(folds, function(fold) {
        sample(fold, length(fold), replace = TRUE)
      }))
      data_boot <- data[bootstrap_indices,]
      tau_boot <- compute_AuDRIE(data_boot$A, data_boot$Y, mu1 = data_boot$mu1, mu0 = data_boot$mu0, pi1 = data_boot$pi1, pi0 =data_boot$pi0)
      return(tau_boot)
    })
    return(NULL)
  })




  CI <- tau_n + quantile(bootstrap_estimates - mean(bootstrap_estimates), c(alpha/2, 1-alpha/2), na.rm = TRUE)
  return(list(estimate = tau_n, CI = CI))
}



compute_AuDRIE <- function(A,Y, mu1, mu0, pi1, pi0) {
  calibrated_estimators  <- calibrate_nuisances(A,Y, mu1 = mu1, mu0 = mu0, pi1 = pi1, pi0 = pi0)
  mu1_star <- calibrated_estimators$mu1_star
  mu0_star <- calibrated_estimators$mu0_star
  mu_star <- ifelse(A==1, mu1_star, mu0_star)
  pi1_star <- calibrated_estimators$pi1_star
  pi0_star <- calibrated_estimators$pi0_star

  alpha_n <- ifelse(A==1, 1/pi1_star, - 1/pi0_star)
  tau_n <-  mean(mu1_star - mu0_star + alpha_n * (Y - mu_star))
  return(tau_n)
}

compute_AIPW <- function(A,Y, mu1, mu0, pi1, pi0) {
  mu <- ifelse(A==1, mu1, mu0)
  alpha_n <- ifelse(A==1, 1/pi1, - 1/pi0)
  tau_n <-  mean(mu1 - mu0 + alpha_n * (Y - mu))
  CI <- tau_n + c(-1,1) * qnorm(1-0.025) * sd(mu1 - mu0 + alpha_n * (Y - mu))/sqrt(n)
  return(list(estimate = tau_n, CI = CI))
}




calibrate_nuisances <- function(A, Y,mu1, mu0, pi1, pi0) {
  calibrator_mu1 <- as.stepfun(isoreg(mu1[A==1], Y[A==1]))
  mu1_star <- calibrator_mu1(mu1)
  calibrator_mu0 <- as.stepfun(isoreg(mu0[A==0], Y[A==0]))
  mu0_star <- calibrator_mu0(mu0)


  calibrator_pi1 <- as.stepfun(isoreg(pi1, A))
  pi1_star <- calibrator_pi1(pi1)
  calibrator_pi0 <- as.stepfun(isoreg(pi0, 1-A))
  pi0_star <- calibrator_pi0(pi0)
  return(list(mu1_star= mu1_star, mu0_star=mu0_star, pi1_star = pi1_star, pi0_star = pi0_star))
}




compute_initial <- function(W,A,Y, lrnr_mu, lrnr_pi) {
  data <- data.table(W,A,Y)
  taskY <- sl3_Task$new(data, covariates = colnames(W), outcome  = "Y", outcome_type = "continuous")
  folds <- taskY$folds
  taskA <- sl3_Task$new(data, covariates = colnames(W), outcome  = "A", folds = folds, outcome_type = "binomial")
  mu1 <- lrnr_mu$train(taskY[A==1])$predict(taskY)
  mu0 <- lrnr_mu$train(taskY[A==0])$predict(taskY)
  pi1 <- lrnr_pi$train(taskA)$predict(taskA)
  pi1 <- truncate_pscore_adaptive(A, pi1)
  return(list(mu1 = mu1, mu0 = mu0, pi1 = pi1, pi0 = 1-pi1, folds = folds))
}

truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  risk_function <- function(cutoff, level) {
    pi <- pmax(pi, cutoff)
    pi <- pmin(pi, 1 - cutoff)
    alpha <- A/pi - (1-A)/(1-pi) #Riesz-representor
    alpha1 <- 1/pi
    alpha0 <- - 1/(1-pi)
    mean(alpha^2 - 2*(alpha1 - alpha0))
  }
  cutoff <- optim(1e-5, fn = risk_function, method = "Brent", lower = min_trunc_level, upper = 0.5, level = 1)$par
  pi <- pmin(pi, 1 - cutoff)
  pi <- pmax(pi, cutoff)
  pi
}





