

#' Doubly-robust nonparametric superefficient estimation of the conditional average treatment effect
#' using the lasso-based R-learner
#'
#' This method estimates the conditional average treatment effect function `w - > tau(w) := E[Y | A=1, W =w] - E[Y | A=0, W = w]`
#' under the regression model `E[Y | A, W] = E[Y | A=0, W] + A * tau(W)`.
#' @param W A \code{numeric} \code{matrix} of covariate values.
#' @param A A \code{numeric} vector of treatment values. May be binary or continuous.
#' @param Y A \code{numeric} vector of outcome values.
#' @param pi.hat A \code{numeric} vector containing estimates of the propensity score `pi(W) := P(A=1 | W)`.
#' @param m.hat A \code{numeric} vector containing estimates of the treatment-marginalized outcome regression `m(W) := E[Y | W]`.
#' @param ... Other arguments to be passed to \code{\link[glmnet]{cv.glmnet}}.
#' @import glmnet
#' @export
#'

fit_lasso_cate_partially_linear <- function(W, A, Y, pi.hat = NULL, m.hat = NULL,    verbose = TRUE,...) {
  if(!is.matrix(W)) W <- as.matrix(W)
  weights <- rep(1, length(Y))

  X <- cbind(W,A)


  ###############################################
  ### Fit CATE tau.hat(W) := E[Y_1 - Y_0 | W] using HAL-based R-learner
  ###############################################
  # R learner is implemented as least-squares regression with the following pseudo-outcomes and pseudo-weights:
  pseudo_outcome <- ifelse(abs(A-pi.hat)<1e-10, 0, (Y - m.hat)/(A-pi.hat))
  pseudo_weights <- (A - pi.hat)^2 * weights


  if(verbose) print("Fitting tau = E[Y|A=1,W] - E[Y|A=0,W]")


  # Observations with near zero weights are dropped for stability.
  keep <- which( abs(A-pi.hat) > 1e-10)
  library(glmnet)
  fit_cate <- cv.glmnet(W[keep,,drop = F], pseudo_outcome[keep],   weights = pseudo_weights[keep], family = "gaussian",   ,... )
  tau.hat <- predict(fit_cate, new_data = W)
  # fit relaxed HAL to allow for plug-in inference
  W_post <- cbind(1, W)[,coef(fit_cate, s = "lambda.min")!=0]
  fit_cate_relaxed <- glm.fit(W_post[keep,,drop = F], pseudo_outcome[keep], family = gaussian(), weights = pseudo_weights[keep], intercept = FALSE)
  beta <- coef(fit_cate_relaxed)
  beta[is.na(beta)] <- 0
  tau.hat_relaxed <- x_basis_tau.hat %*% beta


  # Compute IF map for projection of CATE
  IF_Y_map <- function(x_proj) {
    n <- length(Y)
    gamma <- x_basis_tau.hat %*% solve(t(x_basis_tau.hat) %*% diag((A-pi.hat)^2) %*% x_basis_tau.hat / n) %*% colMeans(x_basis_tau.hat)
    IF_cate <- gamma *  (A-pi.hat) * (Y - m.hat - (A-pi.hat)*tau.hat)
    IF_cate <- IF_cate   + tau.hat - mean(tau.hat)
  }

  fit_cate$internal <- list(IF_Y_map = IF_Y_map, fit_EAW = fit_pi, fit_m.hat = fit_m, data  = list(tau.hat_relaxed = tau.hat_relaxed, W = W, A = A, Y = Y, pi.hat = pi.hat, m.hat = m.hat, tau.hat = tau.hat, pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))


  return(fit_cate)
}
