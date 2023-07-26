library(data.table)
library(sl3)

library(future)
d <- 3

out <- do_sims(5000, 2, 20)


 do_sims <- function(n, pos_const, nsims) {
  loss_inv <- function (pred, observed) {
    A <- observed
    pred <- truncate_pscore_adaptive(A, pred)
    out <- A/pred^2 - 2 * 1/pred + (1-A)/(1-pred)^2 - 2 * 1/(1-pred)
    attributes(out)$name <- "MSE"
    return(out)
  }

  #stack <-  Lrnr_earth$new(degree =1, family = binomial(), pmethod = "cv", nk = 100, nfold = 10)
  #lrnr_hal_inv <- Lrnr_hal9001$new(smoothness_orders = 0, num_knots = c(100,100), max_degree = 2, screen_variables = FALSE, family = "gaussian") # Stack$new(Lrnr_gam$new(family = "gaussian") , Lrnr_gam$new(family = "binomial")) #
  cols <- paste0("W", 1:d)
   #formula_A <- paste0("A~", paste0("s(", cols, ", k = 20, bs='bs',m=c(1,0))", collapse = " + "))
  #formula_Y <- paste0("Y~", paste0("s(", cols, ", k = 20, bs='bs',m=c(1,0))", collapse = " + "))
  formula_A_quad <- paste0("A~", paste0("s(", cols, ", k = 25, bs='bs',m=c(1,1))", collapse = " + "))
  formula_Y_quad <- paste0("Y~", paste0("s(", cols, ", k = 25, bs='bs',m=c(1,1))", collapse = " + "))

  stack_earth_Y <-  Lrnr_earth$new(family = "binomial",   degree=1, pmethod = "cv", nfold = 5,    nk = 75)
  stack_earth_A <-  Lrnr_earth$new(family = "binomial",  degree=1, pmethod = "cv", nfold = 5,   nk = 75)

  stack_gam_Y_quad <-  Lrnr_gam$new(family = "binomial",
                               formula = formula_Y_quad)
  stack_gam_A_quad <-  Lrnr_gam$new(family = "binomial",
                               formula = formula_A_quad)

  stack_rf <- Lrnr_ranger$new()
  stack_xg <- Stack$new(
    list(
      Lrnr_xgboost$new(min_child_weight = max(10, (n)^(1/3)), max_depth = 2, nrounds = 30, eta = 0.15 ),
      Lrnr_xgboost$new(min_child_weight = max(10, (n)^(1/3)), max_depth = 3, nrounds = 30, eta = 0.15 ),
      Lrnr_xgboost$new(min_child_weight = max(10, (n)^(1/3)), max_depth = 4, nrounds = 30, eta = 0.15 ),
      Lrnr_xgboost$new(min_child_weight = max(10, (n)^(1/3)), max_depth = 5, nrounds = 30, eta = 0.15 )

    )
  )


  lrnr_mu_earth <-  Pipeline$new(Lrnr_cv$new(stack_earth_Y), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_pi_earth <- Pipeline$new(Lrnr_cv$new(stack_earth_A), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_mu_gam_quad <-  Pipeline$new(Lrnr_cv$new(stack_gam_Y_quad), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_pi_gam_quad <- Pipeline$new(Lrnr_cv$new(stack_gam_A_quad), Lrnr_cv_selector$new(loss_squared_error))
    lrnr_mu_xg <-  Pipeline$new(Lrnr_cv$new(stack_xg), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_pi_xg <- Pipeline$new(Lrnr_cv$new(stack_xg), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_mu_rf <-  Pipeline$new(Lrnr_cv$new(stack_rf), Lrnr_cv_selector$new(loss_squared_error))
  lrnr_pi_rf <- Pipeline$new(Lrnr_cv$new(stack_rf), Lrnr_cv_selector$new(loss_squared_error))


  sim_results <- lapply(1:nsims, function(i){
    try({
      print(paste0("iter: ", i))
      data_list <- get_data(n, pos_const)
      W <- data_list$W
      A <- data_list$A
      Y <- data_list$Y
      ATE <- data_list$ATE
      n <- length(A)
      initial_estimators_earth <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu_earth, lrnr_pi = lrnr_pi_earth, folds = 5, invert = FALSE)
      folds <- initial_estimators_earth$folds
      initial_estimators_gam_quad <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu_gam_quad, lrnr_pi = lrnr_pi_gam_quad, folds = 5, invert = FALSE)
       initial_estimators_xg <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu_xg, lrnr_pi = lrnr_pi_xg, folds = folds, invert = FALSE)
       initial_estimators_rf <- compute_initial(W,A,Y, lrnr_mu = lrnr_mu_rf, lrnr_pi = lrnr_pi_rf, folds = folds, invert = FALSE)

       initial_estimators_misp <- compute_initial(W,A,Y, lrnr_mu =  Lrnr_cv$new(Lrnr_mean$new()), lrnr_pi =  Lrnr_cv$new(Lrnr_mean$new()), folds = folds)

      out_list <- list()

      #
      # plot(1/initial_estimators$pi1[A==1],1/data_list$pi[A==1] )
      # plot(1/initial_estimators$pi0[A==0],1/(1-data_list$pi[A==0] ))
      # plot(1/calibrated_estimators$pi1[A==1],1/data_list$pi[A==1] )
      # plot(1/calibrated_estimators$pi0[A==0],1/(1-data_list$pi[A==0] ))
      #
      # sqrt(mean((1/initial_estimators$pi1[A==1] - 1/data_list$pi[A==1])^2))
      # sqrt(mean((1/initial_estimators$pi0[A==0] - 1/(1-data_list$pi[A==0]))^2))
      # sqrt(mean((1/calibrated_estimators$pi1[A==1] - 1/data_list$pi[A==1])^2))
      # sqrt(mean((1/calibrated_estimators$pi0[A==0] - 1/(1-data_list$pi[A==0]))^2))


      for(lrnr in c("earth", "gam_1",   "xgboost", "rf")) {
        print(lrnr)
        if(lrnr=="earth") {
          initial_estimators <- initial_estimators_earth
        } else if(lrnr=="gam_1") {
          initial_estimators <- initial_estimators_gam_quad
        } else if(lrnr=="rf") {
          initial_estimators <- initial_estimators_rf
        } else
        {
          initial_estimators <- initial_estimators_xg
        }
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
            mu1 <- initial_estimators_misp$mu1
            mu0 <- initial_estimators_misp$mu0
            pi1 <- initial_estimators_misp$pi1
            pi0 <- initial_estimators_misp$pi0
          }


          out_AIPW <- compute_AIPW(A,Y, mu1=mu1, mu0 =mu0, pi1 = pi1, pi0 = pi0)
          out_AuDRIE <- compute_AuDRIE_boot(A,Y,  mu1=mu1, mu0 =mu0, pi1 = pi1, pi0 = pi0, nboot = 1000, folds = folds, alpha = 0.05)
          #out <- matrix(unlist(c(misp, out_AuDRIE, out_AIPW)), nrow=1)
          out <- as.data.table(rbind(unlist(out_AuDRIE), unlist(out_AIPW)))
          colnames(out) <- c("estimate", "CI_left", "CI_right")
          out$misp <- misp
          out$estimator <- c("auDRI", "AIPW")
          out$lrnr <- lrnr
          out_list[[paste0(misp, lrnr)]] <- out
        }
      }
      out <- rbindlist(out_list)
      out$pos_const <- pos_const
      out$n <- n
      out$ATE <- ATE
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

  W <- replicate(d, runif(n, -1, 1))
  colnames(W) <- paste0("W", 1:d)


  link <- 0.75*(sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]) - 0.5)
  pi0 <- plogis(pos_const * link)
  A <- rbinom(n, 1, pi0)
  mu <- plogis(-1 +
    (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3]))) +
      A * (1 + sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]))
  )
  mu0 <- plogis( -1 +
    (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3]))) +
      0 * (1 + sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]))
  )
  mu1 <- plogis(-1 +
    (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3]))) +
      1* (1 + sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]))
  )
  #mu0 <-  plogis(-link - 0.75)
  #mu1 <- plogis(link + 0.75 + (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3])))/2)
  Y <- rbinom(n, 1, mu)


  out <- list(W=W, A = A, Y = Y,   pi = pi0, mu0 = mu0, mu1 = mu1)

  W <- replicate(d, runif(1000000, -1, 1))
  link <- 0.75*(sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]) - 0.5)
  pi0 <- plogis(pos_const * link)
  A <- rbinom(1000000, 1, pi0)
  mu <- plogis(-1 +
                 (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3]))) +
                 A * (1 + sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]))
  )
  mu0 <- plogis( -1 +
                   (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3]))) +
                   0 * (1 + sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]))
  )
  mu1 <- plogis(-1 +
                  (cos(3.14*W[,1]) + W[,2]*sin(W[,2]) + sqrt(abs(W[,3]))) +
                  1* (1 + sign(W[,1]) * sqrt(abs(W[,1])) + sin(3.14*W[,2]) + W[,3]*sin(W[,3]))
  )
  Y <- rbinom(1000000, 1, A*mu1 + (1-A)*mu0)
  ATE <- mean(mu1 - mu0)
  ATE
  mean(Y[A==1]) - mean(Y[A==0])
  out$ATE <- ATE



  return(out)
}


#
# get_data <- function(n, pos_const) {
#
#   W <- replicate(d, runif(n, -1, 1))
#   colnames(W) <- paste0("W", 1:d)
#
#   beta_1 <- (rnorm(d, 0, 1))
#   beta_1 <- beta_1 / sum(abs(beta_1))
#   beta_2 <- (rnorm(d, 0, 1))
#   beta_2 <- beta_2 / sum(abs(beta_2))
#   pi0 <- plogis(pos_const * (W%*%beta_1 + sin(3.14*W)%*%beta_2) )
#   A <- rbinom(n, 1, pi0)
#   beta_3 <- abs(rnorm(d, 0, 1)) * sign(beta_1)
#   beta_3 <- beta_3 / sum(abs(beta_3))
#   mu0 <-  plogis(qlogis(pi0) - 0.5)
#   mu1 <- plogis(qlogis(mu0) + (1 + 2*cos(3.14*W)%*%beta_3))
#   mu <- ifelse(A==1, mu1, mu0)
#   Y <- rbinom(n, 1, mu)
#
#   Wbig <- replicate(d, runif(1000000, -1, 1))
#   pi0big <- plogis(pos_const * (Wbig%*%beta_1 + sin(3.14*Wbig)%*%beta_2) )
#   mu0big <-  plogis(qlogis(pi0big) - 0.5)
#   mu1big <- plogis(qlogis(mu0big) + (1 + 2*cos(3.14*Wbig)%*%beta_3))
#   ATE <- mean(mu1big - mu0big)
#   print(ATE)
#
#
#   return(list(W=W, A = A, Y = Y, ATE = ATE, pi = pi0, mu0 = mu0, mu1 = mu1))
# }



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
  pi1_star <- calibrated_estimators$pi1_star
  pi0_star <- calibrated_estimators$pi0_star

  mu_star <- ifelse(A==1, mu1_star, mu0_star)
  alpha_n <- ifelse(A==1, 1/pi1_star, - 1/pi0_star)
  tau_n <-  mean(mu1_star - mu0_star + alpha_n * (Y - mu_star))
  return(tau_n)
}

compute_AIPW <- function(A,Y, mu1, mu0, pi1, pi0) {
  pi1 <- pmax(pi1,  25/(sqrt(n)*log(n)))
  pi0 <- pmax(pi0,  25/(sqrt(n)*log(n)))

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
  print("mu")
  taskY0 <- sl3_Task$new(data, covariates = colnames(W), outcome  = "Y", outcome_type = "binomial", folds = folds)
  folds <- taskY0$folds
  fit0 <- lrnr_mu$train(taskY0[A==0])
  mu0 <- fit0$predict(taskY0)
  data$mu0 <- qlogis(mu0)
  taskY1 <- sl3_Task$new(data, covariates = c(colnames(W)), outcome  = "Y", outcome_type = "binomial", folds = folds)
  fit1 <- lrnr_mu$train(taskY1[A==1])
  mu1 <-  fit1$predict(taskY1)


  print("done_mu")
  print("pi")
  taskA <- sl3_Task$new(data, covariates = colnames(W), outcome  = "A", folds = folds, outcome_type = "binomial")

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




