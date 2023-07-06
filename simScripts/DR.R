library(data.table)
library(sl3)
library(doFuture)
library(future)


#Lrnr_gam$new(), Lrnr_earth$new(degree = 1),
#



do_sims <- function(n, pos_const, nsims) {


  stack <- Lrnr_hal9001$new(smoothness_orders = 0, num_knots = 100, max_degree = 1) #Stack$new(Lrnr_hal9001$new(smoothness_orders = 0, num_knots = 50, max_degree = 1))
  lrnr_mu <- Pipeline$new(Lrnr_cv$new(stack), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_pi <- Pipeline$new(Lrnr_cv$new(stack), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_misp_pi <- Lrnr_cv$new(Lrnr_glm$new())
  lrnr_misp_mu <- Lrnr_cv$new(Lrnr_glm$new())


  sim_results <- lapply(1:nsims, function(i){
    try({
      print(paste0("iter: ", i))
      data_list <- get_data(n, pos_const)
      W <- data_list$W
      A <- data_list$A
      Y <- data_list$Y

      initial_estimators <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu, lrnr_pi = lrnr_pi, folds = 5, invert = FALSE)
      folds <- initial_estimators$folds
      initial_estimators_misp <- compute_initial(W,A,Y, lrnr_mu = lrnr_misp_mu, lrnr_pi = lrnr_misp_pi, folds = folds)
      out_list <- list()
      for(misp in c("1", "2", "3", "4")) {
        mu1 <- initial_estimators$mu1
        mu0 <- initial_estimators$mu0
        pi1 <- initial_estimators$pi1
        pi0 <- initial_estimators$pi0
        if(misp == "2" || misp == "4") {
          mu1 <- initial_estimators_misp$mu1
          mu0 <- initial_estimators_misp$mu0
        }
        if(misp == "3" || misp == "4") {
          pi1 <- initial_estimators_misp$pi1
          pi0 <- initial_estimators_misp$pi0
        }


        out_AIPW <- compute_AIPW(A,Y, mu1=mu1, mu0 =mu0, pi1 = pi1, pi0 = pi0)
        out_AuDRIE <- compute_AuDRIE_boot(A,Y,  mu1=mu1, mu0 =mu0, pi1 = pi1, pi0 = pi0, nboot = 1000, folds = folds, alpha = 0.05)
        out <- matrix(unlist(c(misp, out_AuDRIE, out_AIPW)), nrow=1)
        out <- as.data.table(rbind(unlist(out_AuDRIE), unlist(out_AIPW)))
        colnames(out) <- c("estimate", "CI_left", "CI_right")
        out$misp <- misp
        out$estimator <- c("auDRI", "AIPW")
        out_list[[misp]] <- out
      }
      out <- rbindlist(out_list)
      out$iter <- i
      return(as.data.table(out))
    })
    return(data.table())
  })
  sim_results <- data.table::rbindlist(sim_results)
  key <- paste0("DR_iter=", nsims, "_n=", n, "_pos=", pos_const )
  try({fwrite(sim_results, paste0("~/causalHAL/simResultsDR/sim_results_", key, ".csv"))})
  return(sim_results)
}



get_data <- function(n, pos_const) {
  d <- 3
  W <- replicate(d, runif(n, -1, 1))
  colnames(W) <- paste0("W", 1:d)
  pi0 <- plogis(pos_const * ( -0.5+  W[,1]*sin(2*W[,1]) + W[,2]*cos(3*W[,2]) + W[,3]*cos(2*W[,2])    ))
  A <- rbinom(n, 1, pi0)
  mu0 <-  2*(W[,1]*sin(2*W[,1]) + W[,2]*cos(3*W[,2]) + W[,3]*cos(2*W[,2]))
  #mu1 <-  exp(W[,1]) + W[, 2] + sin(4*W[,3]) + abs(W[,3])

  tau <- 0.5 + 2*(W[,1]*sin(2*W[,1]) + W[,2]*cos(3*W[,2]) + W[,3]*cos(2*W[,2]))

  Y <- rnorm(n,  ifelse(A==1, mu0 + tau, mu0), 0.5)
  return(list(W=W, A = A, Y = Y, ATE = 1.371726, pi = pi0, mu0 = mu0, mu1 = mu0 + tau))
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
  pi1 <- truncate_pscore_adaptive(A, pi1)
  pi0 <- truncate_pscore_adaptive(1-A, pi0)
  n <- length(A)
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




compute_initial <- function(W,A,Y, lrnr_mu, lrnr_pi, folds,   invert = TRUE) {
  data <- data.table(W,A,Y)
  taskY <- sl3_Task$new(data, covariates = colnames(W), outcome  = "Y", outcome_type = "continuous", folds = folds)
  folds <- taskY$folds
  taskA <- sl3_Task$new(data, covariates = colnames(W), outcome  = "A", folds = folds, outcome_type = "binomial")
  print("mu")
  mu1 <- lrnr_mu$train(taskY[A==1])$predict(taskY)
  mu0 <- lrnr_mu$train(taskY[A==0])$predict(taskY)
  print("done_mu")
  print("pi")
  pi1 <- lrnr_pi$train(taskA)$predict(taskA)
  if(invert) {
    data0 <- data
    data0$A <= 1-A
    taskA0 <- sl3_Task$new(data0, covariates = colnames(W), outcome  = "A", folds = folds, outcome_type = "continuous")
    pi0 <- lrnr_pi$train(taskA0)$predict(taskA0)
  } else {
    pi0 <- 1 - pi1
  }

  print("done_pi")
  if(invert) {
    pi1 <- 1/pi1
    pi0 <- 1/pi0
  }
  return(list(mu1 = mu1, mu0 = mu0, pi1 = pi1, pi0 = pi0, folds = folds))
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




