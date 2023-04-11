

#' Doubly-robust nonparametric superefficient estimation of the partial conditional average treatment effect
#' using the highly adaptive lasso R-learner
#' with adaptive MARS-based variable and interaction screening.
#'
#' This method estimates the conditional average treatment effect function `w - > tau(w)`
#' under the regression model `E[Y | A, W] = E[Y | A=0, W] + A * tau(W)`.
#' @param W A \code{matrix} of covariate values.
#' @param A A \code{numeric} vector of treatment values.
#' If binary then the CATE is estimated,
#' and otherwise the partial conditional average treatment effect is estimated.
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
fit_hal_pcate <- function(W, A, Y, Delta = NULL, weights = NULL,  lrnr_A= NULL, lrnr_Y = NULL,max_degree = 3, max_degree_pi = max_degree, num_knots = c(sqrt(length(Y)), length(Y)^(1/3), length(Y)^(1/5)), smoothness_orders = 1, family_Y = c("auto", "gaussian", "binomial", "poisson"),  family_A = c("auto", "gaussian", "binomial", "poisson"), formula_cate = NULL, max_degree_cate = max(max_degree - 1, 1), num_knots_cate = num_knots, smoothness_orders_cate = smoothness_orders,    screen_variables = TRUE, screen_control = list(), verbose = TRUE,...) {
  if (!inherits(family_A, "family")) family_A <- match.arg(family_A)
  if (!inherits(family_Y, "family")) family_Y <- match.arg(family_Y)
  family_A <- detect_family(A, family_A); family_Y <- detect_family(Y, family_Y)

  if(!is.matrix(W)) W <- as.matrix(W)
  if(is.null(Delta)) Delta <- rep(1, length(Y))
  if(is.null(weights)) weights <- rep(1, length(Y))

  X <- cbind(W,A)
  subset <- which(Delta == 1)
  # mu = E[Y | W]
  if(verbose) print("Fitting mu = E[Y|W]")
  if(is.null(lrnr_Y)) {
  fit_mu <- fit_hal(W[subset, , drop = F], Y[subset], weights = weights[subset], max_degree = max_degree, num_knots = num_knots, smoothness_orders = smoothness_orders, family = family_Y, screen_variables = screen_variables, screen_control = screen_control,... )
  cols_mu <- unique(unlist(lapply(fit_mu$basis_list, function(basis) {basis$cols})))
  mu <- predict(fit_mu, new_data = W)
  if(verbose) print(fit_mu$formula)
  } else {
    task_Y <- sl3_Task$new(data.table(W, A = A, Y= Y), covariates = c(colnames(W), "A"), outcome = "Y", outcome_type = "continuous")
    fit_mu <- lrnr_Y$train(task_Y)
    mu <-fit_mu$predict(task_Y)
  }

  # pi <- E[A]
  if(verbose) print("Fitting pi = E[A|W]")
  if(is.null(lrnr_A)) {
  fit_pi <- fit_hal(W, A, weights = weights, max_degree = max_degree_pi, num_knots = num_knots, smoothness_orders = smoothness_orders, family = family_A, screen_variables = screen_variables, screen_control = screen_control,... )
  cols_pi <- unique(unlist(lapply(fit_pi$basis_list, function(basis) {basis$cols})))
  pi <- predict(fit_pi, new_data = W)
  } else {
    task_A <- sl3_Task$new(data.table(W, A = A), covariates = colnames(W), outcome = "A", outcome_type = "binomial")
    fit_pi <- lrnr_A$train(task_A)
    pi <-fit_pi$predict(task_A)
  }
  # R-learner for CATE
 # vars_cate <- colnames(W)[sort(cols_mu)] # For now just use mu variables
  pseudo_outcome <- ifelse(abs(A-pi)<1e-10, 0, (Y - mu)/(A-pi))
  pseudo_weights <- (A - pi)^2 * weights
  if(verbose) print(fit_pi$formula)

  if(verbose) print("Fitting tau = E[Y|A=1,W] - E[Y|A=0,W]")


  # Observations with near zero weights are dropped.
  keep <- which( abs(A-pi) > 1e-10)
  # mu1 <- mu + (1 - pi) * tau
  # mu0 <- mu + (0 - pi) * tau
  # mu1 - mu0 + (A - pi) * (Y-mu-(A-pi)*tau)
  fit_cate <- fit_hal(W[keep,,drop = F], pseudo_outcome[keep], formula = formula_cate, weights = pseudo_weights[keep], max_degree = max_degree_cate, num_knots = num_knots_cate, smoothness_orders = smoothness_orders_cate, family = "gaussian", screen_variables = FALSE, screen_control = list(screen_interactions = FALSE),... )
  tau <- predict(fit_cate, new_data = W)
  basis_list_reduced_tau <- fit_cate$basis_list[fit_cate$coefs[-1] != 0]
  x_basis_tau <- cbind(1,as.matrix(hal9001::make_design_matrix(W, basis_list_reduced_tau)))
  print(mean(tau))
  if(verbose) print(fit_cate$formula)

  fit_cate_relaxed <- glm.fit(x_basis_tau[keep,,drop = F], pseudo_outcome[keep], family = gaussian(), weights = pseudo_weights[keep], intercept = FALSE)
  beta <- coef(fit_cate_relaxed)
  beta[is.na(beta)] <- 0
  tau_relaxed <- x_basis_tau %*% beta




  # IF_conditional <- (A - pi) * (Y - mu - (A-pi)*tau)

  # Compute IF map for projection of CATE
  IF_Y_map <- function(x_proj) {
    n <- length(Y)
    scale <- solve(t(x_proj) %*% diag(pseudo_weights) %*% x_proj / n)
    IF_cate <- (x_proj) *  pseudo_weights* as.vector(pseudo_outcome - tau_relaxed)
    IF_cate <- IF_cate %*% scale + tau_relaxed - mean(tau_relaxed)
  }

  fit_cate$internal <- list(IF_Y_map = IF_Y_map, fit_EAW = fit_pi, fit_EYW = fit_mu, data  = list(tau_relaxed = tau_relaxed, W = W, A = A, Y = Y, pi = pi, mu = mu, tau = tau, pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))
  class(fit_cate) <- c("hal9001", "hal_cate")
  return(fit_cate)
}
