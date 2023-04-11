
#' Estimates and confidence intervals for the ATE.
#' @param fit_cate A \code{hal_cate} object obtained from the function \code{fit_hal_cate}.
#' @param alpha Significant level for confidence intervals
#' @param return_cov_mat A \code{logical} for whether to return the asymptotic covariance matrix of the coefficient estimates.
#' @export
inference_ate <- function(fit_cate,  alpha = 0.05, return_cov_mat = FALSE) {
  return(inference_cate(fit_cate, formula =  ~ 1, alpha, return_cov_mat))
}

#' Estimates and confidence intervals for the projection of the CATE onto a user-specified parametric working model.
#' @param fit_cate A \code{hal_cate} object obtained from the function \code{fit_hal_cate}.
#' @param formula A \code{formula} object specifying a working parametric model for the conditional average treatment effect.
#' For instance, `formula = ~ 1` specifies the marginal average treatment effect `E[CATE(W)]`.
#' More complex formula like `formula = ~ W1` specifies the best `W1`-linear approximation of the true CATE.
#' @param alpha Significant level for confidence intervals
#' @param return_cov_mat A \code{logical} for whether to return the asymptotic covariance matrix of the coefficient estimates.
#' @export
inference_cate <- function(fit_cate, formula =  ~ 1, alpha = 0.05, return_cov_mat = FALSE) {
  internal <- fit_cate$internal
  data <- internal$data
  pseudo_weights <- data$pseudo_weights
  pseudo_outcome <- data$pseudo_outcome
  A <- data$A
  sandwich_weights <- data$sandwich_weights
  if(is.null(sandwich_weights)) sandwich_weights <- pseudo_weights

  coefs <- fit_cate$coefs
  basis_list <- fit_cate$basis_list[coefs[-1]!=0]
  coefs <- coefs[coefs!=0]
  x_basis <- cbind(1,as.matrix(hal9001::make_design_matrix(as.matrix(data$W), basis_list)))
  tau <- x_basis %*% coef(glm.fit(x_basis, pseudo_outcome, weights = pseudo_weights))

  x_basis_proj <- model.matrix(formula, data = as.data.frame(data$W ))
  coef_proj <-  coef(glm.fit(x_basis_proj, tau))
  tau_proj <- x_basis_proj %*% coef_proj

  n <- nrow(x_basis)

  scale <- solve(t(x_basis_proj) %*% x_basis_proj / n) # 1 if intercept model
  IF_proj <-   x_basis_proj * as.vector(tau - tau_proj)  # residual if intercept model
  IF_proj <- IF_proj %*% scale

  IF_cate <- internal$IF_Y_map(x_basis_proj)
  # x_basis_proj = x_basis then we recover nonprojection case
  #scale2 <- solve(t(x_basis_proj) %*% diag(pseudo_weights) %*% x_basis_proj / n)
  #IF_cate <- (x_basis_proj) *  pseudo_weights* as.vector(pseudo_outcome - tau)
  #IF_cate <- IF_cate %*% scale2


  IF_full <- IF_cate #+ IF_proj

  print(sd(IF_cate)/ sqrt(n))
  print(sd(IF_proj)/ sqrt(n) )
  var_mat <- var(IF_full)
  se <- sqrt(diag(var_mat)) / sqrt(n)
  print(se)
  CI <- matrix(coef_proj + abs(qnorm(alpha/2)) * c(-1,1) * se, ncol =2)

  summary <- data.table(variable = colnames(x_basis_proj), coef = coef_proj, se = se, CI_left = CI[,1],  CI_right = CI[,2])
  if(!return_cov_mat) return(summary)
  return(list(summary = summary, cov_mat = var_mat))

}




