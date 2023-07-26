

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
  tau.hat <- predict(fit_cate, newx = W)
  # fit relaxed HAL to allow for plug-in inference
  nonzero <- as.vector(coef(fit_cate, s = "lambda.min"))!=0
  W_post <- cbind(1, W)[,nonzero]
  fit_cate_relaxed <- glm.fit(W_post[keep,,drop = F], pseudo_outcome[keep], family = gaussian(), weights = pseudo_weights[keep], intercept = FALSE)
  beta <- coef(fit_cate_relaxed)
  beta[is.na(beta)] <- 0
  tau.hat_relaxed <- W_post %*% beta


  # Compute IF map for projection of CATE
  IF_Y_map <- function(x_proj) {
    n <- length(Y)
    gamma <- W_post %*% solve(t(W_post) %*% diag((A-pi.hat)^2) %*% W_post / n) %*% colMeans(W_post)
    IF_cate <- gamma *  (A-pi.hat) * (Y - m.hat - (A-pi.hat)*tau.hat)
    IF_cate <- IF_cate   + tau.hat - mean(tau.hat)
  }

  fit_cate$internal <- list(IF_Y_map = IF_Y_map,    data  = list(tau.hat_relaxed = tau.hat_relaxed, W = W, A = A, Y = Y, pi.hat = pi.hat, m.hat = m.hat, tau.hat = tau.hat, pseudo_outcome = pseudo_outcome, pseudo_weights = pseudo_weights))


  return(fit_cate)
}
