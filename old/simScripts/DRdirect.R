library(data.table)
library(sl3)
library(doFuture)
library(future)


#Lrnr_gam$new(), Lrnr_earth$new(degree = 1),
#



out <- do_sims(3000, 0.5, 20)

do_sims <- function(n, pos_const, nsims) {
  loss_inv <- function (pred, observed) {
    A <- observed
    pred <- truncate_pscore_adaptive(A, pred)
    out <- A/pred^2 - 2 * 1/pred + (1-A)/(1-pred)^2 - 2 * 1/(1-pred)
    attributes(out)$name <- "MSE"
    return(out)
  }

  stack <- Stack$new(Lrnr_xgboost$new(max_depth = 2), Lrnr_xgboost$new(max_depth = 3), Lrnr_xgboost$new(max_depth = 4))
  lrnr_hal <- Lrnr_hal9001$new(smoothness_orders = 0, num_knots = 100, max_degree = 1, screen_variables = FALSE, family = "gaussian") #Stack$new(Lrnr_hal9001$new(smoothness_orders = 0, num_knots = 50, max_degree = 1))
  lrnr_gam_mu <- Lrnr_gam$new(family = "binomial")
  lrnr_hal_inv <- Lrnr_hal9001$new(smoothness_orders = 0, num_knots = c(100), max_degree = 1, screen_variables = FALSE, family = "gaussian") # Stack$new(Lrnr_gam$new(family = "gaussian") , Lrnr_gam$new(family = "binomial")) #

  lrnr_mu <- Pipeline$new(Lrnr_cv$new(stack), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_pi <- Pipeline$new(Lrnr_cv$new(stack), Lrnr_cv_selector$new(loss_inv))
  lrnr_misp_pi <- Lrnr_cv$new(Lrnr_mean$new())
  lrnr_misp_mu <- Lrnr_cv$new(Lrnr_mean$new())

  lrnr_list <- get_learners()
  lrnr_mu <- lrnr_list$mu
  lrnr_pi <- lrnr_list$pi
  lrnr_misp_pi <- lrnr_list$misp_pi
  lrnr_misp_mu <- lrnr_list$misp_mu
  sim_results <- lapply(1:nsims, function(i){
    try({
      print(paste0("iter: ", i))
      data_list <- get_data(n, pos_const)
      W <- data_list$W
      A <- data_list$A
      Y <- data_list$Y

      initial_estimators <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu, lrnr_pi = lrnr_pi, folds = 5, invert = TRUE)
      folds <- initial_estimators$folds
      initial_estimators_misp <- compute_initial(W,A,Y, lrnr_mu =  Lrnr_cv$new(Lrnr_mean$new()), lrnr_pi =  Lrnr_cv$new(Lrnr_mean$new()), folds = folds)
      initial_estimators_misp2 <- compute_initial(W,A,Y, lrnr_mu = lrnr_misp_mu, lrnr_pi = lrnr_misp_pi, folds = folds)

      out_list <- list()
      #print(mean(initial_estimators$mu1 - initial_estimators$mu0))
      #print(mean(A*Y/initial_estimators$pi1) - mean((1-A)*Y/(initial_estimators$pi0)))


      calibrated_estimators  <- calibrate_nuisances(A,Y, mu1 = initial_estimators$mu1, mu0 = initial_estimators$mu0, pi1 = initial_estimators$pi1, pi0 = initial_estimators$pi0)

      plot(1/initial_estimators$pi1[A==1],1/data_list$pi[A==1] )
      plot(1/initial_estimators$pi0[A==0],1/(1-data_list$pi[A==0] ))
      plot(1/calibrated_estimators$pi1[A==1],1/data_list$pi[A==1] )
      plot(1/calibrated_estimators$pi0[A==0],1/(1-data_list$pi[A==0] ))

      sqrt(mean((1/initial_estimators$pi1[A==1] - 1/data_list$pi[A==1])^2))
      sqrt(mean((1/initial_estimators$pi0[A==0] - 1/(1-data_list$pi[A==0]))^2))
      sqrt(mean((1/calibrated_estimators$pi1[A==1] - 1/data_list$pi[A==1])^2))
      sqrt(mean((1/calibrated_estimators$pi0[A==0] - 1/(1-data_list$pi[A==0]))^2))



      for(misp in c("1", "2", "3", "4")) {
        mu1 <- initial_estimators$mu1
        mu0 <- initial_estimators$mu0
        pi1 <- initial_estimators$pi1
        pi0 <- initial_estimators$pi0
        if(misp == "2") {
          mu1 <- initial_estimators_misp$mu1
          mu0 <- initial_estimators_misp$mu0
        } else if(misp == "3" ) {
          pi1 <- initial_estimators_misp$pi1
          pi0 <- initial_estimators_misp$pi0
        } else if(misp == "4") {
          mu1 <- initial_estimators_misp2$mu1
          mu0 <- initial_estimators_misp2$mu0
          pi1 <- initial_estimators_misp2$pi1
          pi0 <- initial_estimators_misp2$pi0
        }


        out_AIPW <- compute_AIPW(A,Y, mu1=mu1, mu0 =mu0, pi1 = pi1, pi0 = pi0)
        out_AuDRIE <- compute_AuDRIE_boot(A,Y,  mu1=mu1, mu0 =mu0, pi1 = pi1, pi0 = pi0, nboot = 1000, folds = folds, alpha = 0.05)
        #out <- matrix(unlist(c(misp, out_AuDRIE, out_AIPW)), nrow=1)
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
  d <- 10
  W <- replicate(d, runif(n, -1, 1))
  colnames(W) <- paste0("W", 1:d)
  #pos_consts <- c(0.25, 0.75,  1.25, 1.75)
  #f <- sapply(pos_consts, function(pos){
  # pi0 <- plogis(pos * ( W[,1] +  cos(3*W[,2]) +  sin(3*W[,3])    ))
  #  print(range(pi0))
  #})
  pi0 <- plogis(pos_const * ( W[,1] + W[,2] + W[,3]   ))

  A <- rbinom(n, 1, pi0)
  mu0 <-  W[,1] - W[,2] + W[,3]  #(sin(3*W[,1]) + W[,2] + cos(3*W[,3]))
  #mu1 <-  exp(W[,1]) + W[, 2] + sin(4*W[,3]) + abs(W[,3])

  tau <- 1 - W[,1] + W[,2] + W[,3]

  Y <- rnorm(n,  mu0 + A * tau, 0.5)


  return(list(W=W, A = A, Y = Y, ATE = 1.092115, pi = pi0, mu0 = mu0, mu1 = mu0 + tau))
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




  CI <- tau_n + quantile(bootstrap_estimates - median(bootstrap_estimates), c(alpha/2, 1-alpha/2), na.rm = TRUE)
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
  #pi1 <- truncate_pscore_adaptive(A, pi1)
  #pi0 <- truncate_pscore_adaptive(1-A, pi0)
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




compute_initial <- function(W,A,Y, lrnr_mu, lrnr_pi, folds,   invert = FALSE) {
  data <- data.table(W,A,Y)
  taskY <- sl3_Task$new(data, covariates = colnames(W), outcome  = "Y", outcome_type = "continuous", folds = folds)
  folds <- taskY$folds
  taskA <- sl3_Task$new(data, covariates = colnames(W), outcome  = "A", folds = folds, outcome_type = "continuous")
  print("mu")
  mu1 <- lrnr_mu$train(taskY[A==1])$predict(taskY)
  mu0 <- lrnr_mu$train(taskY[A==0])$predict(taskY)
  print("done_mu")
  print("pi")
  fit1 <- lrnr_pi$train(taskA)
  pi1 <- fit1$predict(taskA)
  print(fit1$fit_object$learner_fits$Lrnr_cv_selector_NULL$fit_object$cv_risk)
  if(invert) {
    print("invert")
    data0 <- data
    data0$A <- 1-A
    taskA0 <- sl3_Task$new(data0, covariates = colnames(W), outcome  = "A", folds = folds, outcome_type = "continuous")
    fit0 <- lrnr_pi$train(taskA0)
    pi0 <- fit0$predict(taskA0)
    print(fit0$fit_object$learner_fits$Lrnr_cv_selector_NULL$fit_object$cv_risk)


    pi1 <- 1/pi1
    pi0 <- 1/pi0
  } else {
    pi0 <- 1 - pi1
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




