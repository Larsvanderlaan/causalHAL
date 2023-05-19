

#' Doubly-robust nonparametric superefficient estimation of the partial conditional average treatment effect
#' using the highly adaptive lasso R-learner
#' with adaptive MARS-based variable and interaction screening.
#'
#' This method estimates the conditional average treatment effect function `w - > tau(w) := E[Y | A=1, W =w] - E[Y | A=0, W = w]`
#' under the regression model `E[Y | A, W] = E[Y | A=0, W] + A * tau(W)`.
#' @param W A \code{numeric} \code{matrix} of covariate values.
#' @param A A \code{numeric} vector of treatment values. May be binary or continuous.
#' @param Y A \code{numeric} vector of outcome values.
#' @param weights (Optional) A \code{numeric} vector of observation weights.
#' @param max_degree For estimation of nuisance functions `E[Y|W]` and `E[X|W]`.
#' The maximm interaction degree of basis functions generated.
#' Passed to \code{\link[hal9001]{fit_hal}} function of \code{hal9001} package.
#' @param num_knots For estimation of nuisance functions `E[Y|W]` and `E[X|W]`.
#' Passed to \code{\link[hal9001]{fit_hal}} function of \code{hal9001} package.
#' A \code{numeric} vector of length \code{max_degree} where
#' the `d`-th entry specifies the number of univariable spline knot points to use
#' when generating the tensor-product basis functions of interaction degree `d`.
#' @param smoothness_orders For estimation of nuisance functions `E[Y|W]` and `E[X|W]`.
#' An integer taking values in (0,1,2,...)
#' specifying the smoothness order of the basis functions. See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param family_Y A \code{\link[stats]{family}} object specifying the outcome type of the outcome \code{Y}.
#' This is passed internally to \code{\link[hal9001]{fit_hal}} when estimating `E[Y | W]`.
#' @param family_A A \code{\link[stats]{family}} object specifying the outcome type of the treatment \code{A}.
#' This is passed internally to \code{\link[hal9001]{fit_hal}} when estimating `E[A | W]`.
#' @param formula_cate (Optional) A \code{hal9001}-formatted \code{formula} object for the CATE/tau to be passed to \code{\link[hal9001]{formula_hal}}.
#' By default the CATE model is learned data-adaptivelly using MARS-based screening and HAL.
#' @param max_degree_cate (Optional) Same as \code{max_degree} but for CATE model.
#' @param num_knots_cate (Optional) Same as \code{num_knots} but for CATE model.
#' @param smoothness_orders_cate (Optional) Same as \code{smoothness_orders} but for CATE model.
#' @param screen_variables Highly recommended. See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param screen_interactions Highly recommended. See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param ... Other arguments to be passed to \code{\link[hal9001]{fit_hal}}.
#' @import hal9001
#' @export
#'

 fit_hal_cate_Rlearner <- function(W, A, Y, weights = NULL,  pA1W = NULL, EYW = NULL, formula_cate = NULL, max_degree_cate = max(max_degree - 1, 1), num_knots_cate = num_knots, smoothness_orders_cate = smoothness_orders, screen_variables_cate = FALSE, lrnr_A= NULL, lrnr_Y = NULL,   params_hal_EYW = list(smoothness_orders =1, max_degree = 3, num_knots = c(sqrt(length(Y)), length(Y)^(1/3), length(Y)^(1/5)), screen_variables = TRUE), params_hal_pA1W =  params_hal_EYW   , verbose = TRUE,...) {
  if(!is.matrix(W)) W <- as.matrix(W)
  if(is.null(weights)) weights <- rep(1, length(Y))

  X <- cbind(W,A)

  ###############################################
  ### Fit outcome regression m(W) := E[Y | W]
  ###############################################
  if(verbose) print("Fitting m = E[Y|W]")
  # If given use it.
  if(!is.null(EYW)) {
    m <- EYW
    fit_m <- NULL
  } else if(is.null(lrnr_Y)) {
    # if no sl3 learner provides use HAL
    params_hal_EYW$X <- W
    params_hal_EYW$Y <- Y
    params_hal_EYW$weights <- weights
    params_hal_EYW$family <- "gaussian"
    fit_m <- sl3:::call_with_args(fit_hal, params_hal_EYW)
    #fit_m <- fit_hal(W[subset, , drop = F], Y[subset], weights = weights[subset], max_degree = max_degree, num_knots = num_knots, smoothness_orders = smoothness_orders, family = family_Y, screen_variables = screen_variables, screen_control = screen_control,... )
    cols_m <- unique(unlist(lapply(fit_m$basis_list, function(basis) {basis$cols})))
    m <- predict(fit_m, new_data = W)
  if(verbose) print(fit_m$formula)
  } else {
    # use sl3 learner if provided
    task_Y <- sl3_Task$new(data.table(W, Y= Y), covariates = c(colnames(W)), outcome = "Y", outcome_type = "continuous")
    fit_m <- lrnr_Y$train(task_Y)
    m <- fit_m$predict(task_Y)
  }

  ###############################################
  ### Fit propensity score pi(W) := E[A | W]
  ###############################################
  if(verbose) print("Fitting pi = E[A|W]")
  if(!is.null(pA1W)) {
    pi <- pA1W
    fit_pi <- NULL
  } else if(is.null(lrnr_A)) {
    params_hal_pA1W$X <- W
    params_hal_pA1W$Y <- A
    params_hal_pA1W$weights <- weights
    params_hal_pA1W$family <- "gaussian"
    fit_pi <- sl3:::call_with_args(fit_hal, params_hal_pA1W)
    cols_pi <- unique(unlist(lapply(fit_pi$basis_list, function(basis) {basis$cols})))
    pi <- predict(fit_pi, new_data = W)
  } else {
    task_A <- sl3_Task$new(data.table(W, A = A), covariates = colnames(W), outcome = "A", outcome_type = "continuous")
    fit_pi <- lrnr_A$train(task_A)
    pi <-fit_pi$predict(task_A)
  }
  #### truncate propensity score to lie in (c_n,1-c_n) with data-adaptive truncation level c_n.
  pi <- truncate_pscore_adaptive(A, pi)

  ###############################################
  ### Fit CATE tau(W) := E[Y_1 - Y_0 | W] using HAL-based R-learner
  ###############################################
  # R learner is implemented as least-squares regression with the following pseudo-outcomes and pseudo-weights:
  pseudo_outcome <- ifelse(abs(A-pi)<1e-10, 0, (Y - m)/(A-pi))
  pseudo_weights <- (A - pi)^2 * weights


  if(verbose) print("Fitting tau = E[Y|A=1,W] - E[Y|A=0,W]")


  # Observations with near zero weights are dropped.
  keep <- which( abs(A-pi) > 1e-10)
  fit_cate <- fit_hal(W[keep,,drop = F], pseudo_outcome[keep], formula = formula_cate, weights = pseudo_weights[keep], max_degree = max_degree_cate, num_knots = num_knots_cate, smoothness_orders = smoothness_orders_cate, family = "gaussian", screen_variables = screen_variables_cate, screen_control = list(screen_interactions = screen_variables_cate),... )
  tau <- predict(fit_cate, new_data = W)
  # fit relaxed HAL to allow for plug-in inference
  basis_list_reduced_tau <- fit_cate$basis_list[fit_cate$coefs[-1] != 0]
  x_basis_tau <- cbind(1,as.matrix(hal9001::make_design_matrix(W, basis_list_reduced_tau)))
  print(mean(tau))
  if(verbose) print(fit_cate$formula)

  fit_cate_relaxed <- glm.fit(x_basis_tau[keep,,drop = F], pseudo_outcome[keep], family = gaussian(), weights = pseudo_weights[keep], intercept = FALSE)
  beta <- coef(fit_cate_relaxed)
  beta[is.na(beta)] <- 0
  tau_relaxed <- x_basis_tau %*% beta




  # IF_conditional <- (A - pi) * (Y - m - (A-pi)*tau)

  # Compute IF map for projection of CATE
  IF_Y_map <- function(x_proj) {
    n <- length(Y)
    scale <- solve(t(x_proj) %*% diag(pseudo_weights) %*% x_proj / n)
    IF_cate <- (x_proj) *  pseudo_weights* as.vector(pseudo_outcome - tau_relaxed)
    IF_cate <- IF_cate %*% scale + tau_relaxed - mean(tau_relaxed)
  }

  fit_cate$internal <- list(IF_Y_map = IF_Y_map, fit_EAW = fit_pi, fit_EYW = fit_m, data  = list(tau_relaxed = tau_relaxed, W = W, A = A, Y = Y, pi = pi, m = m, tau = tau, pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
  class(fit_cate) <- c("hal9001", "hal_cate")
  return(fit_cate)
}
