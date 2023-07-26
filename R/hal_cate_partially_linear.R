

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
#' @param pi.hat A \code{numeric} vector containing estimates of the propensity score `pi(W) := P(A=1 | W)`.
#' @param m.hat A \code{numeric} vector containing estimates of the treatment-marginalized outcome regression `m(W) := E[Y | W]`.
#' @param sl3_Lrnr_pi.hat If \code{pi.hat} is not provided, a \code{\link[sl3]{Lrnr_base} object for estimation of `pi(W) := P(A=1 | W)`.
#' @param sl3_Lrnr_m.hat If \code{m.hat} is not provided, a \code{\link[sl3]{Lrnr_base} object for estimation of `m(W) := E[Y | W]`.
#' @param formula_cate (Optional) A \code{hal9001}-formatted \code{formula} object for the CATE/tau to be passed to \code{\link[hal9001]{formula_hal}}.
#' By default the CATE model is learned data-adaptively using MARS-based screening and HAL. See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param max_degree_cate (Optional) Same as \code{max_degree} but for CATE model.  See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param num_knots_cate (Optional) Same as \code{num_knots} but for CATE model.  See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param smoothness_orders_cate (Optional) Same as \code{smoothness_orders} but for CATE model.  See documentation for \code{\link[hal9001]{fit_hal}}.
#' @param ... Other arguments to be passed to \code{\link[hal9001]{fit_hal}}.
#' @import hal9001
#' @export
#'

 fit_hal_cate_partially_linear <- function(W, A, Y,   weights = NULL, pi.hat = NULL, m.hat = NULL, formula_cate = NULL, max_degree_cate = 1, smoothness_orders_cate = 1, num_knots_cate = c(50),    sl3_Lrnr_pi.hat= NULL, sl3_Lrnr_m.hat = NULL, verbose = TRUE,...) {
  if(!is.matrix(W)) W <- as.matrix(W)
  if(is.null(weights)) weights <- rep(1, length(Y))

  X <- cbind(W,A)

  ###############################################
  ### Fit outcome regression m(W) := E[Y | W]
  ###############################################
  if(verbose) print("Fitting m = E[Y|W]")
  # If given use it.
  if(!is.null(m.hat)) {
    m <- m.hat
    fit_m <- NULL
  }  else {
    # use sl3 learner if provided
    task_Y <- sl3_Task$new(data.table(W, Y= Y), covariates = c(colnames(W)), outcome = "Y", outcome_type = "continuous")
    fit_m <- sl3_Lrnr_m.hat$train(task_Y)
    m <- fit_m$predict(task_Y)
  }

  ###############################################
  ### Fit propensity score pi(W) := E[A | W]
  ###############################################
  if(verbose) print("Fitting pi = E[A|W]")
  if(!is.null(pi.hat)) {
    pi <- pi.hat
    fit_pi <- NULL
  } else {
    task_A <- sl3_Task$new(data.table(W, A = A), covariates = colnames(W), outcome = "A", outcome_type = "continuous")
    fit_pi <- sl3_Lrnr_pi.hat$train(task_A)
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
  fit_cate <- fit_hal(W[keep,,drop = F], pseudo_outcome[keep], formula = formula_cate, weights = pseudo_weights[keep], max_degree = max_degree_cate, num_knots = num_knots_cate, smoothness_orders = smoothness_orders_cate, family = "gaussian", screen_variables = FALSE,  ... )
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
    gamma <- x_basis_tau %*% solve(t(x_basis_tau) %*% diag((A-pi)^2) %*% x_basis_tau / n) %*% colMeans(x_basis_tau)
    IF_cate <- gamma *  (A-pi) * (Y - m - (A-pi)*tau)
    IF_cate <- IF_cate   + tau - mean(tau)
  }

  fit_cate$internal <- list(IF_Y_map = IF_Y_map, fit_EAW = fit_pi, fit_m.hat = fit_m, data  = list(tau_relaxed = tau_relaxed, W = W, A = A, Y = Y, pi = pi, m = m, tau = tau, pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
  class(fit_cate) <- c("hal9001", "hal_cate")
  return(fit_cate)
}
