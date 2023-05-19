
# turn to plug-in


#' Doubly-robust nonparametric superefficient estimation of the conditional average treatment effect
#' using the highly adaptive lasso plug-in estimator.
#'
#' This method estimates the conditional average treatment effect function `w - > tau(w)`
#' under the regression model `E[Y | A, W] = E[Y | A=0, W] + A * tau(W)`.
#' @param W A \code{matrix} of covariate values.
#' @param A A \code{numeric} binary vector of treatment values.
#' @param Y A \code{numeric} vector of outcome values.
#' @param Delta (Not used)
#' @param weights (Optional) A \code{numeric} vector of observation weights.
#' @param max_degree For estimation of nuisance functions `E[Y|W]` and `E[X|W]`.
#' The maximum interaction degree of basis functions generated.
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

fit_hal_cate <- function(W, A, Y,weights = NULL, formula_cate = NULL, max_degree_cate = 3, num_knots_cate =  c(sqrt(length(Y)), length(Y)^(1/3), length(Y)^(1/5)), smoothness_orders_cate = 1, screen_variable_cate = TRUE,   params_EY0W =  list(max_degree = 3, num_knots =  c(sqrt(length(Y)), length(Y)^(1/3), length(Y)^(1/5)), smoothness_orders = 1, screen_variables = TRUE), include_propensity_score = FALSE,   verbose = TRUE,...) {
    if(!all(A %in% c(0,1))) {
    stop("The treatment `A` should be binary. For continuous treatments consider using `fit_hal_pcate()`.")
  }

  if(!is.matrix(W)) W <- as.matrix(W)

  if(is.null(weights)) weights <- rep(1, length(Y))

  X <- cbind(W,A)

  # mu = E[Y | W]
  if(include_propensity_score) {
    if(verbose) print("Fitting pi = E[A|W]")
    fit_pi <- fit_hal(W, A, weights = weights, max_degree = max_degree, num_knots = num_knots, smoothness_orders = smoothness_orders, family = family_A, screen_variables = screen_variables,... )
    cols_pi <- unique(unlist(lapply(fit_pi$basis_list, function(basis) {basis$cols})))
    pi <- predict(fit_pi, new_data = W)
    if(verbose) print(fit_pi$formula)
  } else {
    pi <- NULL
    fit_pi <- NULL
  }


  if(verbose) print("Fitting mu0 = E[Y|A=0,W]")
  if(include_propensity_score){
    W_pi <- cbind(W,pi)
  } else {
    W_pi <- W
  }
  subset <- A==0
  params_EY0W$X <- W_pi[subset, , drop = F]
  params_EY0W$Y <- Y[subset]
  params_EY0W$weights <- weights[subset]
  params_EY0W$family <- "gaussian"

  fit_mu0 <- sl3:::call_with_args(fit_hal, params_EY0W)
  #fit_mu0 <- fit_hal(W_pi[subset, , drop = F], Y[subset], weights = weights[subset], max_degree = max_degree, num_knots = num_knots, smoothness_orders = smoothness_orders, family = "gaussian", screen_variables = screen_variables,... )
  mu0 <- predict(fit_mu0, new_data = W_pi)
  basis_list_reduced_mu0 <- fit_mu0$basis_list[fit_mu0$coefs[-1] != 0]
  x_basis_mu0 <- cbind(1,as.matrix(hal9001::make_design_matrix(W_pi, basis_list_reduced_mu0)))
  if(verbose) print(fit_mu0$formula)


  if(verbose) print("Fitting tau = E[Y|A=1,W] - E[Y|A=0,W]")
  subset <- A==1
  # Offset not supported by MARS screener so added to outcome since we use least-squares
  fit_cate <- fit_hal(W[subset, , drop = F], Y[subset] - mu0[subset], weights = weights[subset], max_degree = max_degree_cate, formula = formula_cate, num_knots = num_knots_cate, smoothness_orders = smoothness_orders_cate, family = "gaussian", screen_variables = screen_variable_cate,... )
  tau <- predict(fit_cate, new_data = W)
  basis_list_reduced_tau <- fit_mu0$basis_list[fit_cate$coefs[-1] != 0]
  x_basis_tau <- cbind(1,as.matrix(hal9001::make_design_matrix(W, basis_list_reduced_tau)))
  if(verbose) print(fit_cate$formula)


  mu <- mu0 + A * tau
  mu1 <-mu0 + tau

  x_basis <- cbind(x_basis_mu0, A * x_basis_tau)
  x_basis0 <- cbind(x_basis_mu0, 0 * x_basis_tau)
  x_basis1 <- cbind(x_basis_mu0, 1 * x_basis_tau)


  fit_mu_relaxed <- glm.fit(x_basis, Y, family = gaussian(), weights = weights, intercept = FALSE)
  beta_mu <- coef(fit_mu_relaxed)
  mu1_relaxed <- x_basis1 %*% beta_mu
  mu0_relaxed <- x_basis0 %*% beta_mu
  tau_relaxed <- mu1_relaxed - mu0_relaxed


  # Computes Y component of IF for projection of CATE onto x_proj
  IF_Y_map <- function(x_proj) {
    n <- length(A)
   #M <- x_proj %*% solve((t(x_proj) %*% x_proj) / n, t(x_proj))
    alpha <- x_basis %*% solve(t(x_basis) %*% x_basis/n, colMeans(x_basis1 - x_basis0)  )
    IF <- alpha * (Y - x_basis %*% beta_mu) +  (x_basis1 - x_basis0) %*% beta_mu
  }

  # For inference function
  pseudo_outcome <- tau_relaxed
  pseudo_weights <- weights

  fit_cate$internal <- list(IF_Y_map = IF_Y_map, fit_EAW = fit_pi, fit_EY0W = fit_mu0, fit_cate = fit_cate, data  = list(tau_relaxed = tau_relaxed, W = W, A = A, Y = Y, pi = pi, mu1 = mu1, mu0 = mu0, tau = tau, pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
  class(fit_cate) <- c("hal9001", "hal_cate")
  return(fit_cate)
}




