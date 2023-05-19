


fit_hal_crate <- function(W, A, Y, weights = NULL, family_Y = c("gaussian", "binomial", "poisson"), formula_crate = NULL, max_degree_crate = 3, num_knots_crate =  c(sqrt(length(Y)), length(Y)^(1/3), length(Y)^(1/5)), smoothness_orders_crate = 1, screen_variable_crate = TRUE,   params_EYAW =  list(max_degree = 3, num_knots =  c(sqrt(length(Y)), length(Y)^(1/3), length(Y)^(1/5)), smoothness_orders = 1, screen_variables = TRUE),  verbose = TRUE,...) {
  family_Y <- match.arg(family_Y)
  params_EYAW$X <- cbind(W,A)
  params_EYAW$Y <- Y
  params_EYAW$weights <- weights
  params_EYAW$family <- family_Y
  mu_fit <- sl3:::call_with_args(fit_hal, params_EYAW)
  mu1 <- predict(mu_fit, cbind(W,1))
  mu0 <- predict(mu_fit, cbind(W,0))
  mu1 <- pmax(mu1, 1e-8)
  mu0 <- pmax(mu0, 1e-8)
  pseudo_weights <- mu1 + mu0
  pseudo_outcome <- mu1/(mu0 + mu1)
  fit_crate <- fit_hal(W, pseudo_outcome, weights= pseudo_weights, family = binomial(), formula = formula_crate, max_degree = max_degree_crate,  num_knots = num_knots_crate, smoothness_orders = smoothness_orders_crate, screen_variable_crate = screen_variable_crate, ...  )
  class(fit_crate) <- c("hal9001", "hal_crate")
}
