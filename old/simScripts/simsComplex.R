
library(autocausalML)
library(data.table)
library(sl3)
library(doMC)
library(hal9001)

doMC::registerDoMC(cores = 5)

run_sims <- function(const, n, nsims, fit_control = list(), formula_hal = ~ h(.) + h(.,A), num_knots = c(1,1), smoothness_orders = 1, max_degree = 2,screen_basis = F, gen_fun, lrnr_pi = Lrnr_glmnet$new(), lrnr_g = Lrnr_glmnet$new(formula = ~ . + A * .),nboots = 500, relaxed_fit = TRUE, weight_screen_by_alpha = FALSE) {
  print(const)
  print(n)
  fit_hal_g_params = list(   smoothness_orders = smoothness_orders, max_degree =max_degree, num_knots = num_knots)

  out_list <- list()




  fit_control$parallel = TRUE
   fit_hal_g_params$fit_control <- fit_control
  fit_hal_g_params$num_knots <- num_knots
  fit_hal_g_params$smoothness_orders <- smoothness_orders
  fit_hal_g_params$formula <- formula_hal
  fit_hal_g_params$max_degree <- max_degree
  fit_hal_g_params$fit_control$relax <- FALSE
  fit_hal_g_params$fit_control$gamma <- 0
  #fit_hal_g_params$fit_control$weights <- weights


  fit_hal_g_params$family = "gaussian"
  fit_hal_g_params$reduce_basis = 25/n

  # fit_hal_g_params_relaxed <- fit_hal_g_params

  # hal_fit <- sl3:::call_with_args( hal9001::fit_hal, fit_hal_g_params)
  # lambda <- hal_fit$lambda_star
  # fit_hal_g_params$lambda <- lambda
  fit_hal_g_params$fit_control$cv_select = TRUE
  fit_hal_g_params$fit_control$parallel = TRUE
  fit_hal_g_params_relaxed <- fit_hal_g_params


  fit_hal_g_params_relaxed$lambda <- NULL
  fit_hal_g_params_relaxed$fit_control$relax <- TRUE


  out <- lapply(1:nsims, function(iter) {
    try({
      print("Current settings:")
      print(n)
      print(const)
      print(iter)
      print(weight_screen_by_alpha)
      datam_list <- gen_fun(const)(n)

      X <- datam_list$X
      A <- datam_list$A
      Y <- datam_list$Y
      # fit_control$parallel = TRUE
      #fit_hal_g_params$fit_control$foldid <- (sample(1:n,n, replace= FALSE) %% 10) + 1
      # fit_hal_g_params$fit_control <- fit_control
      # fit_hal_g_params$num_knots <- num_knots
      # fit_hal_g_params$smoothness_orders <- smoothness_orders
      # fit_hal_g_params$formula <- formula_hal
      # fit_hal_g_params$max_degree <- max_degree
      #fit_hal_g_params$fit_control$foldid <-  fit_hal_g_params$fit_control$foldid

      g_basis_gen <-make_g_basis_generator_HAL(X,A,Y,  fit_hal_g_params = fit_hal_g_params,  screen_basis = screen_basis, relaxed_fit = FALSE, weight_screen_by_alpha = weight_screen_by_alpha)
    print("done")
      # weights <- g_basis_gen$weights
      #g_basis_gen <- g_basis_gen$g_basis
     # g_basis_gen_relaxed <-make_g_basis_generator_HAL(X,A,Y,  fit_hal_g_params = fit_hal_g_params_relaxed, screen_basis = screen_basis, relaxed_fit = TRUE, weight_screen_by_alpha = weight_screen_by_alpha)
      #g_basis_gen_relaxed <- g_basis_gen_relaxed$g_basis
      #print(quantile(weights))
      ATE <- datam_list$ATE

      causal_sieve <- causalsieve$new(X, A, Y, g_basis_gen, nboots = nboots, weights = NULL)
      causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1, name = "ATE")
      causal_sieve$estimate()
      name <- unlist(sapply(causal_sieve$estimates, `[[`, "name"))

      estimates <- unlist(sapply(causal_sieve$estimates, `[[`, "estimate"))
      CI_IF_df <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      causal_sieve$confint(include_se_df_correction = FALSE)
      CI_IF <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))

      CI_boot <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI_boot"))


      out <- cbind(t(as.data.table(c(iter, name))), t(as.data.table(as.numeric(c(estimates, CI_IF, CI_IF_df, CI_boot)))))
      colnames(out) <- c("iter", "name", "estimate", "CI_left", "CI_right", "CI_df_left", "CI_df_right", "CI_boot_left", "CI_boot_right")


      # RELAXED
      # causal_sieve <- causalsieve$new(X, A, Y, g_basis_gen_relaxed, nboots = nboots, weights = NULL)
      # causal_sieve$add_target_parameter(g(A=1,X=X) - g(A=0,X=X) ~ 1, name = "ATE")
      # causal_sieve$estimate()
      # name <- unlist(sapply(causal_sieve$estimates, `[[`, "name"))
      #
      # estimates <- unlist(sapply(causal_sieve$estimates, `[[`, "estimate"))
      # CI_IF_df <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      # causal_sieve$confint(include_se_df_correction = FALSE)
      # CI_IF <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI"))
      # CI_boot <- do.call(rbind, lapply(causal_sieve$estimates, `[[`, "CI_boot"))
      # out2 <- cbind(t(as.data.table(c(iter, name))), t(as.data.table(as.numeric(c(estimates, CI_IF, CI_IF_df, CI_boot)))))
      # colnames(out2) <- c("iter", "name", "estimate_relaxed", "CI_left_relaxed", "CI_right_relaxed", "CI_df_left_relaxed", "CI_df_right_relaxed", "CI_boot_left_relaxed", "CI_boot_right_relaxed")
      # print("HERE")

      g1 <- causal_sieve$regression_fit$g_n1
      g0 <- causal_sieve$regression_fit$g_n0

      data <- as.data.frame(cbind(X,A,Y))
      # g_ests <- compute_g(data, lrnr_g = lrnr_g)
      # g1 <- g_ests$g1
      #g0 <- g_ests$g0
      pi <-compute_pi(as.data.frame(cbind(X,A,Y)), lrnr_pi = lrnr_pi)
      tmle <- compute_TMLE (data, pi, g1, g0,level = 0.05)
      aipw <- compute_AIPW (data, pi, g1, g0,level = 0.05)
      lm <- compute_glm (data,level = 0.05)

      comp <- as.numeric(unlist(c(tmle[-2], aipw[-2], lm[-2])))
      names(comp) <- c(paste0(c("estimate", "CI_left", "CI_right"), "_tmle"),
                       paste0(c("estimate", "CI_left", "CI_right"), "_aipw"),
                       paste0(c("estimate", "CI_left", "CI_right"), "_lm"))
      #
     # out <- cbind(out,out2[,-c(1,2), drop = F], t(as.data.table(comp)))
      out <- cbind(out, t(as.data.table(comp)))

      colnames(out) <- c(
        c("iter", "name", "estimate", "CI_left", "CI_right", "CI_df_left", "CI_df_right", "CI_boot_left", "CI_boot_right"),
      #  c("estimate_relaxed", "CI_left_relaxed", "CI_right_relaxed", "CI_df_left_relaxed", "CI_df_right_relaxed", "CI_boot_left_relaxed", "CI_boot_right_relaxed"),
        c(paste0(c("estimate", "CI_left", "CI_right"), "_tmle"),
          paste0(c("estimate", "CI_left", "CI_right"), "_aipw"),
          paste0(c("estimate", "CI_left", "CI_right"), "_lm"))
      )
      out_list[[iter]] <<- out
      out_full <- as.data.table(do.call(rbind, out_list))









      return(out)
    })
  })

  out <- as.data.frame(do.call(rbind, out_list))
  out$const <- const
  out$n <- n
  return(out)

}




### COMPLEX
const <- as.numeric(const)
n <- as.numeric(n)
nsims <- 1000
library(sl3)
fit_control <- list()
if(n == 5000){
  nknots <- 100
  fit_control$lambda.min.ratio <- 1e-5
} else if(n == 4000){
  nknots <- 80
  fit_control$lambda.min.ratio <- 1e-5
} else if(n == 3000){
  nknots <- 65
  fit_control$lambda.min.ratio <- 1e-4
} else if(n == 2000){
  nknots <- 50
  fit_control$lambda.min.ratio <- 1e-4
} else if(n == 1000){
  nknots <- 35
  fit_control$lambda.min.ratio <- 1e-4
} else if(n == 500){
  nknots <- 20
  fit_control$lambda.min.ratio <- 1e-4
}



fit_control$parallel = TRUE

out <- run_sims(const,n,nsims, fit_control = fit_control, formula_hal = ~ h(.) + h(.,A), num_knots = c(nknots,nknots), screen_basis = TRUE, gen_fun = get_data_generator_nonlinear, lrnr_pi = Lrnr_gam$new(),
                lrnr_g = Lrnr_hal9001$new(formula = ~h(.)  , fit_control = fit_control, smoothness_orders = 1, max_degree =1, num_knots = c(nknots, 1)), nboots=2)

#out_list[[as.character(const)]][[as.character(n)]] <<- out
#out2 <- rbindlist(unlist(out_list, recursive = F))
fwrite(out, file = paste0("ComplexParametricHAL_", const,"_" ,n, ".csv"))
#
# outs <- lapply(c( 3,5,8), function(const) {
#    lapply(rev (c(   500, 1000,  2500 ,5000 )) ,function(n) {
#     fit_control <- list()
#     if(n >= 4000){
#       nknots <- 100
      fit_control$lambda.min.ratio <- 1e-5
#     } else if(n == 2500){
#       nknots <- 75
#       fit_control$lambda.min.ratio <- 1e-4
#     } else if(n == 1000){
#       nknots <- 50
#       fit_control$lambda.min.ratio <- 1e-4
#     } else if(n == 500){
#       nknots <- 30
#       fit_control$lambda.min.ratio <- 1e-4
#     }
#
#
#     fit_control$parallel = TRUE
#
#     out <- run_sims(const,n,2500, fit_control = fit_control, formula_hal = ~ h(.) + h(.,A), num_knots = c(nknots,nknots), screen_basis = TRUE, gen_fun = get_data_generator_nonlinear, lrnr_pi = Lrnr_gam$new(),
#                     lrnr_g = Lrnr_hal9001$new(formula = ~h(.)  , fit_control = fit_control, smoothness_orders = 1, max_degree =1, num_knots = c(nknots, 1)), nboots=2)
#
#     #out_list[[as.character(const)]][[as.character(n)]] <<- out
#     #out2 <- rbindlist(unlist(out_list, recursive = F))
#     fwrite(out, file = paste0("ComplexParametricHAL_", const,"_" ,n, ".csv"))
#     return(NULL)
#
#   })
#
# })
#
